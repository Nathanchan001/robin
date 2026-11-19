"""
固定轴（fixed-tilt）光伏发电模型
================================

代码化自《南北半球光伏固定轴发电模型_v1.0.xlsx》，逐时 8760 口径，南北半球通用
（纬度带符号，南纬为负；方位角 0°=正北，90°=正东，180°=正南，与 Excel 一致）。

计算链条（与 Excel「主计算表(北)」「1.1温度计算表(固定轴)」逐列对应）：

1. 太阳位置：Spencer 赤纬角 + EoT 时差 + 真太阳时角；
   cos(θz) 按 Excel 口径钳制为 ≥0.0104（夜间取 +0.0104，不影响发电量）
2. 大气层外辐照（Spencer 公式，1367 W/m² 基准）
3. POA 转译：
   - Beam：IBeam = |GHI − DHI|（水平面 beam），Rb 投影比按 Excel 规则——
     入射角 >90° 时上午（小时≤12）取 |cosθ|/cosθz、下午取 0；否则 cosθ/cosθz
   - 散射四模型可选：PVsyst（默认，即 AS 列所用）、Hay&Davis、Reindl、
     简化环日+各向同性（AX 列）；地面反射反照率 0.2
4. IAM：ASHRAE beam IAM（AOI>87° 为 0，1−0.05·(1/cosθ−1)），只修正 beam 分量
5. 温度（Sandia）：Tm = Epoa·e^(a+b·风速) + 环境温度；
   Tcell = Tm + Epoa/1000·ΔTcnd；功率温度系数 −0.0035/℃（25℃ 基准），
   安装方式参数见 MOUNTING_COEFFICIENTS（默认 Glass/cell/polymer sheet + Open rack）
6. 首年衰减：逐小时线性（默认 2%/8759）；第 n 年再乘 (1 − 0.5%·(n−1)) 长期衰减
7. 综合损失默认 −8%（线损1.5 失配2 污秽2 阴影0.5 连接0.5 停用0.5 铭牌1，%）
8. 逆变器效率默认 0.98

气象输入为当地时间（政治时区）逐时序列，8760 小时（非闰年）。
"""
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

# Sandia 阵列温度系数（a, b, ΔTcnd），键 = (电池板放置方式, 安装方式)，
# 与 Excel「基本参数表」H1:L7 六行一致
MOUNTING_COEFFICIENTS = {
    ('glass_cell_glass', 'open_rack'): (-3.47, -0.0594, 3.0),
    ('glass_cell_glass', 'close_roof_mount'): (-2.98, -0.0471, 1.0),
    ('glass_cell_polymer', 'open_rack'): (-3.56, -0.075, 3.0),
    ('glass_cell_polymer', 'insulated_back'): (-2.81, -0.0455, 0.0),
    ('polymer_thinfilm_steel', 'open_rack'): (-3.58, -0.113, 3.0),
}

DEFAULT_LOSSES = {  # 综述表 C15:C24（阴影遮挡为固定经验值，与几何自遮挡无关）
    '线损': -0.015, '失配': -0.02, '污秽': -0.02, '阴影遮挡': -0.005,
    '连接': -0.005, '停用': -0.005, '铭牌': -0.01,
}

# 散射模型与 Excel 列对应：AS=Epoa 用 pvsyst；AV/AW/AX 为备选输出
DIFFUSE_MODELS = ('pvsyst', 'haydavies', 'reindl', 'isotropic_circumsolar')

ALBEDO = 0.2


@dataclass
class FixedAxisResult:
    """逐时结果与年度汇总。功率单位 kWh/h，辐照单位 W/m²。"""
    hours: int
    # 逐时
    cos_zenith: np.ndarray = field(repr=False, default=None)
    aoi: np.ndarray = field(repr=False, default=None)          # 入射角（度）
    rb: np.ndarray = field(repr=False, default=None)           # Rb 投影比
    beam: np.ndarray = field(repr=False, default=None)         # Eb 倾斜面直射
    diffuse: np.ndarray = field(repr=False, default=None)      # Ed 倾斜面散射（所选模型）
    ground: np.ndarray = field(repr=False, default=None)       # Eg 地面反射
    epoa: np.ndarray = field(repr=False, default=None)         # AS 倾斜面总辐照
    iam_poa: np.ndarray = field(repr=False, default=None)      # AL IAM 修正后
    cell_temp: np.ndarray = field(repr=False, default=None)    # Tcell（℃）
    power: np.ndarray = field(repr=False, default=None)        # BA 最终逐时发电量
    # 年度汇总（对应主计算表 F6:G22）
    total_ghi: float = None
    total_epoa: float = None
    total_beam: float = None
    total_diffuse: float = None
    total_ground: float = None
    energy_raw: float = None        # G17 = ΣEpoa·kWp/1000
    energy: float = None            # G18 = ΣBA（最终发电量）
    ratio_iam_temp: float = None    # G19 = G18/G17
    effective_hours: float = None   # 有效日照时数 = energy/kWp/365

    def annual_by_year(self, year):
        """第 year 年发电量（首年=1），长期衰减 0.5%/年，同 Excel I17 口径。"""
        return self.energy * (1 - 0.005 * (year - 1))


def _solar_position(dayofyear, hour, latitude, longitude, tz_meridian):
    """返回 cosθz(钳制)、θz(度)、方位角(度, 北为0顺时针)、EoT 用的时角等。"""
    n = np.asarray(dayofyear, float)
    h = np.asarray(hour, float)
    frac = (2.0 * np.pi / 365.0) * (n - 1 + (h - 12.0) / 24.0)
    decl = np.degrees(0.006918 - 0.399912 * np.cos(frac) + 0.070257 * np.sin(frac)
                      - 0.006758 * np.cos(2 * frac) + 0.000907 * np.sin(2 * frac)
                      - 0.002697 * np.cos(3 * frac) + 0.00148 * np.sin(3 * frac))
    eot = 229.18 * (0.000075 + 0.001868 * np.cos(frac) - 0.032077 * np.sin(frac)
                    - 0.014615 * np.cos(2 * frac) - 0.040849 * np.sin(2 * frac))
    offset_min = 4.0 * (longitude - tz_meridian) + eot          # S 列
    true_solar_time = h + offset_min / 60.0                     # U 列
    hour_angle = 15.0 * (true_solar_time - 12.0)                # V 列（度）
    d = np.radians(decl)
    phi = np.radians(latitude)
    w = np.radians(hour_angle)
    cos_z = np.sin(d) * np.sin(phi) + np.cos(d) * np.cos(phi) * np.cos(w)  # W 列
    # Excel 口径：|cosθz|≤0.0104 时取 +0.0104（夜间为正，仅影响中间量，不影响电量）
    cos_z = np.where(np.abs(cos_z) <= 0.0104, 0.0104, cos_z)
    zenith = np.degrees(np.arccos(cos_z))                       # X 列
    sin_z = np.sin(np.radians(zenith))
    # 方位角（Z/AA 列，北为 0，顺时针）
    z_arg = -(np.sin(phi) * cos_z - np.sin(d)) / (np.cos(phi) * sin_z)
    z_arg = np.clip(z_arg, -1.0, 1.0)
    az = np.where(hour_angle > 0, 360.0 - np.degrees(np.arccos(z_arg)),
                  np.degrees(np.arccos(z_arg)))
    return cos_z, zenith, az


def simulate_fixed_axis(ghi, dhi, dni, wind_speed, temp_air, dayofyear, hour,
                        latitude, longitude, tilt, surface_azimuth, capacity_kwp,
                        tz_meridian=120.0, mounting=('glass_cell_polymer', 'open_rack'),
                        temp_coeff=-0.0035, temp_ref=25.0,
                        first_year_degradation=0.02, inverter_efficiency=0.98,
                        losses=None, diffuse_model='pvsyst',
                        apply_iam=True, apply_temperature=True, year=1):
    """固定轴光伏逐时发电仿真（南北半球通用，口径同 Excel 模型 v1.0）。

    Parameters
    ----------
    ghi, dhi, dni : array-like
        水平总辐照、散射辐照、法向直射辐照（W/m²，当地时间逐时）。
    wind_speed, temp_air : array-like
        风速（m/s）与环境温度（℃）。
    dayofyear, hour : array-like
        年积日（1—365）与小时（0—23），与气象数组等长，共 8760 条（非闰年）。
    latitude, longitude : float
        纬度（北纬为正）、经度（东经为正），度。
    tilt, surface_azimuth : float
        倾角（度）与方位角（0=北、90=东、180=南，度）。
    capacity_kwp : float
        装机容量（kWp）。
    tz_meridian : float, default 120
        政治时区中央经线（中国 UTC+8 = 120；Excel 由时区对照表查得）。
    mounting : tuple, default ('glass_cell_polymer', 'open_rack')
        MOUNTING_COEFFICIENTS 的键，决定 Sandia 温度系数 (a, b, ΔTcnd)。
    temp_coeff, temp_ref : float
        功率温度系数（%/℃，负值）与基准温度（℃）。
    first_year_degradation : float, default 0.02
        首年线性衰减率（逐小时 L=rate/8759）。
    inverter_efficiency : float, default 0.98
        逆变器效率。
    losses : dict, optional
        综合损失分项（负值），默认 DEFAULT_LOSSES 合计 −0.08；
        传 {} 表示不加综合损失。
    diffuse_model : str, default 'pvsyst'
        'pvsyst' | 'haydavies' | 'reindl' | 'isotropic_circumsolar'。
    apply_iam, apply_temperature : bool, default True
        是否应用 IAM 修正 / 温度修正（对应 Excel BA/BB/BC/BD 四种组合）。
    year : int, default 1
        第几年（首年=1），>1 时乘 (1−0.005·(year−1)) 长期衰减。
    """
    ghi = np.asarray(ghi, float); dhi = np.asarray(dhi, float)
    dni = np.asarray(dni, float); wind = np.asarray(wind_speed, float)
    tamb = np.asarray(temp_air, float)
    n = len(ghi)
    if not (len(dhi) == len(dni) == len(wind) == len(tamb) == n):
        raise ValueError('气象数组长度不一致')

    cos_z, zenith, az = _solar_position(dayofyear, hour, latitude, longitude, tz_meridian)
    beta = np.radians(tilt)
    # 大气层外辐照（AE 列，Spencer）
    ad = 2.0 * np.pi * np.asarray(dayofyear, float) / 365.0
    ext = 1367.0 * (1.00011 + 0.034221 * np.cos(ad) + 0.00128 * np.sin(ad)
                    + 0.000719 * np.cos(2 * ad) + 0.000077 * np.sin(2 * ad))
    # 入射角（AM 列）
    cos_aoi = cos_z * np.cos(beta) + np.sin(np.radians(zenith)) * np.sin(beta) \
        * np.cos(np.radians(az - surface_azimuth))
    cos_aoi = np.clip(cos_aoi, -1.0, 1.0)
    aoi = np.degrees(np.arccos(cos_aoi))

    # Rb（AQ 列，数组公式口径）：AOI>90° 时上午取 |cosθ|/cosθz、下午取 0
    rb_day = cos_aoi / cos_z
    rb_backside_am = np.abs(cos_aoi) / cos_z
    is_backside = aoi > 90.0
    rb = np.where(is_backside & (np.asarray(hour) <= 12), rb_backside_am,
                  np.where(is_backside, 0.0, rb_day))

    beam_h = np.abs(ghi - dhi)                     # AI 列 IBeam（水平 beam）
    beam = np.abs(rb * beam_h)                     # AT 列 Eb
    anisotropy = dni / ext                         # AR 列 Ai（Hay&Davis/Reindl 用）
    safe_cos_z = np.where(np.abs(cos_z) < 1e-6, 1e-6, cos_z)

    if diffuse_model == 'pvsyst':                  # AU 列（AS/Epoa 默认组成）
        f_bh = beam_h / (1368.0 * safe_cos_z)
        diffuse = dhi * ((1 - f_bh) * (1 + np.cos(beta)) * 0.5 + f_bh * rb)
    elif diffuse_model == 'haydavies':             # AV 列
        diffuse = dhi * (rb * anisotropy + (1 - anisotropy) * (1 + np.cos(beta)) / 2)
    elif diffuse_model == 'reindl':                # AW 列
        with np.errstate(divide='ignore', invalid='ignore'):
            circ = np.sqrt(dni * cos_z / np.where(ghi > 0, ghi, np.nan)) \
                * np.sin(beta / 2) ** 3
        circ = np.nan_to_num(circ)
        diffuse = dhi * (cos_aoi * anisotropy
                         + (1 - anisotropy) * (1 + np.cos(beta)) / 2 * (1 + circ))
    elif diffuse_model == 'isotropic_circumsolar':  # AX 列
        diffuse = np.abs(dhi * (1 + np.cos(beta)) / 2
                         + ghi * 0.5 * (0.012 * zenith - 0.04) * (1 - np.cos(beta)))
    else:
        raise ValueError(f'未知散射模型：{diffuse_model}（可选 {DIFFUSE_MODELS}）')

    ground = ghi * ALBEDO * (1 - np.cos(beta)) * 0.5   # AY 列 Eg
    epoa = beam + diffuse + ground                     # AS 列 Epoa

    # IAM（AK/AL 列）：ASHRAE beam IAM，AOI>87° 为 0；只修正 beam
    iam_beam = np.where(aoi > 87.0, 0.0, 1 - 0.05 * (1 / np.where(
        np.abs(cos_aoi) < 1e-6, 1e-6, cos_aoi) - 1))
    iam_poa = iam_beam * beam + diffuse + ground
    irradiance_use = iam_poa if apply_iam else epoa    # AL 或 AS

    # 温度（1.1温度计算表，Sandia；Epoa 用未 IAM 修正值，与 Excel 一致）
    a_m, b_m, dt_cnd = MOUNTING_COEFFICIENTS[mounting]
    t_m = epoa * np.exp(a_m + b_m * wind) + tamb     # G 列
    t_cell = t_m + (epoa / 1000.0) * dt_cnd          # H 列
    temp_factor = 1 + temp_coeff * (t_cell - temp_ref)

    # 首年逐小时线性衰减（N 列，L=rate/8759）
    t_idx = np.arange(n, dtype=float)
    degradation = 1 - t_idx * (first_year_degradation / 8759.0)

    loss_total = sum((losses if losses is not None else DEFAULT_LOSSES).values())
    power = (1 + loss_total) * (irradiance_use * capacity_kwp / 1000.0) * degradation
    if apply_temperature:
        power = power * temp_factor
    power = power * inverter_efficiency
    if year > 1:  # 长期衰减 0.5%/年（I17 口径）
        power = power * (1 - 0.005 * (year - 1))

    energy_raw = float((epoa * capacity_kwp / 1000.0).sum())      # G17
    energy = float(power.sum())                                   # G18
    return FixedAxisResult(
        hours=n, cos_zenith=cos_z, aoi=aoi, rb=rb, beam=beam, diffuse=diffuse,
        ground=ground, epoa=epoa, iam_poa=iam_poa, cell_temp=t_cell, power=power,
        total_ghi=float(ghi.sum()), total_epoa=float(epoa.sum()),
        total_beam=float(beam.sum()), total_diffuse=float(diffuse.sum()),
        total_ground=float(ground.sum()),
        energy_raw=energy_raw, energy=energy,
        ratio_iam_temp=energy / energy_raw if energy_raw else None,
        effective_hours=energy / capacity_kwp / 365.0,
    )


def effective_hours(energy_kwh, capacity_kwp):
    """有效日照时数（h/day）= 年发电量 / 装机容量 / 365，同综述表 J3。"""
    return energy_kwh / capacity_kwp / 365.0


_CHINA_LOOKUP_PATH = Path(__file__).resolve().parent / 'data' / 'china_pv_lookup.csv'


def load_china_lookup(path=_CHINA_LOOKUP_PATH):
    """读取《中国各地区光伏最佳安装倾角及峰值日照速查表》。

    Returns
    -------
    list of dict
        每行含 region/province/city/tilt_deg/peak_sun_hours/
        yield_kwh_per_w/annual_effective_hours/daily_effective_hours。
    """
    import csv
    with open(path, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def lookup_city(city, path=_CHINA_LOOKUP_PATH):
    """按城市名查速查表，返回该城市条目（dict）；未找到返回 None。"""
    for row in load_china_lookup(path):
        if row['city'] == city:
            return row
    return None
