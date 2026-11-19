"""Validation of the fixed-tilt generation model against
《南北半球光伏固定轴发电模型_v1.0.xlsx》 (cached Excel values).

The Excel workbook ships an NSRDB 2021 hourly dataset for 24.54°N, 112.11°E
(exported to ``tests/data/nsrdb_2021_24.54N_112.11E.csv``).  With tilt=10°,
azimuth=260°, 146.76 kWp, the Excel's cached results are:

- Total GHI      = 1,485,256.8        W/m²
- Total Epoa     = 1,591,106.0433     W/m²
- Total Eb       = 825,653.2676       W/m²
- Total Ed(PVsyst) = 763,196.3369     W/m²
- Total Eg       = 2,256.4388         W/m²
- ΣAZ (raw)      = 233,510.7229       kWh
- ΣBA (final, 综述表「最终发电量」) = 183,985.3958  kWh
- G19 ratio      = 0.78790984
"""
import csv
import datetime
from pathlib import Path

import numpy as np
import pytest

from twoaxistracking.fixed_axis import (
    simulate_fixed_axis, effective_hours, lookup_city,
)

DATA = Path(__file__).parent / 'data' / 'nsrdb_2021_24.54N_112.11E.csv'


@pytest.fixture(scope='module')
def weather():
    cols = [[] for _ in range(8)]
    with open(DATA, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            for i, key in enumerate(['month', 'day', 'hour', 'dhi', 'dni',
                                     'ghi', 'wind_speed', 'temp_air']):
                cols[i].append(float(row[key]))
    mo, da, hr, dhi, dni, ghi, ws, ta = cols
    doy = [datetime.date(2021, int(m), int(d)).timetuple().tm_yday
           for m, d in zip(mo, da)]
    return (np.array(ghi), np.array(dhi), np.array(dni), np.array(ws),
            np.array(ta), np.array(doy), np.array(hr))


@pytest.fixture(scope='module')
def result(weather):
    ghi, dhi, dni, ws, ta, doy, hr = weather
    return simulate_fixed_axis(
        ghi, dhi, dni, ws, ta, doy, hr,
        latitude=24.54, longitude=112.11, tilt=10.0, surface_azimuth=260.0,
        capacity_kwp=146.76, tz_meridian=120.0)


REL_TOL = 1e-6


def test_irradiance_totals(result):
    assert result.total_ghi == pytest.approx(1485256.8, rel=REL_TOL)
    assert result.total_epoa == pytest.approx(1591106.0433336997, rel=REL_TOL)
    assert result.total_beam == pytest.approx(825653.2676303632, rel=REL_TOL)
    assert result.total_diffuse == pytest.approx(763196.3368887474, rel=REL_TOL)
    assert result.total_ground == pytest.approx(2256.438814589762, rel=REL_TOL)


def test_final_energy(result):
    # 综述表「最终发电量」= ΣBA
    assert result.energy == pytest.approx(183985.39580011467, rel=REL_TOL)
    # G17 = ΣEpoa·kWp/1000（无 IAM/温度/损失/衰减/逆变器修正）
    assert result.energy_raw == pytest.approx(233510.72291965422, rel=REL_TOL)
    assert result.ratio_iam_temp == pytest.approx(0.7879098377140475, rel=REL_TOL)


def test_effective_hours_and_year15(result):
    assert effective_hours(result.energy, 146.76) == \
        pytest.approx(3.4346523407915015, rel=REL_TOL)
    # I17 = G18·(1 − 0.5%·14)
    assert result.annual_by_year(15) == \
        pytest.approx(171106.41809410663, rel=REL_TOL)


def test_toggle_combinations(weather):
    """BA/BB/BC/BD 四列：IAM 与温度修正的四种开关组合。"""
    ghi, dhi, dni, ws, ta, doy, hr = weather
    common = dict(latitude=24.54, longitude=112.11, tilt=10.0,
                  surface_azimuth=260.0, capacity_kwp=146.76, tz_meridian=120.0)
    bb = simulate_fixed_axis(ghi, dhi, dni, ws, ta, doy, hr,
                             apply_iam=False, apply_temperature=False, **common)
    bc = simulate_fixed_axis(ghi, dhi, dni, ws, ta, doy, hr,
                             apply_iam=False, apply_temperature=True, **common)
    bd = simulate_fixed_axis(ghi, dhi, dni, ws, ta, doy, hr,
                             apply_iam=True, apply_temperature=False, **common)
    assert bb.energy == pytest.approx(208522.97880077065, rel=REL_TOL)   # G20
    assert bc.energy == pytest.approx(192815.54429931645, rel=REL_TOL)   # G22
    assert bd.energy == pytest.approx(198172.12675706734, rel=REL_TOL)   # G24


def test_diffuse_model_variants(weather):
    """AV/AW/AX 备选散射模型年度总量。"""
    ghi, dhi, dni, ws, ta, doy, hr = weather
    common = dict(latitude=24.54, longitude=112.11, tilt=10.0,
                  surface_azimuth=260.0, capacity_kwp=146.76, tz_meridian=120.0)
    hd = simulate_fixed_axis(ghi, dhi, dni, ws, ta, doy, hr,
                             diffuse_model='haydavies', **common)
    rl = simulate_fixed_axis(ghi, dhi, dni, ws, ta, doy, hr,
                             diffuse_model='reindl', **common)
    assert hd.total_diffuse == pytest.approx(633248.1500599388, rel=REL_TOL)  # K8
    assert rl.total_diffuse == pytest.approx(633248.1500599388, rel=REL_TOL)  # K9


def test_china_lookup():
    bj = lookup_city('北京')
    assert bj is not None
    assert float(bj['tilt_deg']) == 35
    assert float(bj['peak_sun_hours']) == pytest.approx(4.21)
    assert lookup_city('不存在的城市') is None
