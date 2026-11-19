# Open source code for calculating self-shading of two-axis tracking solar collectors
twoaxistracking is a python package for simulating two-axis tracking solar collectors, particularly self-shading.

## Fixed-tilt generation model

`twoaxistracking.fixed_axis` implements a fixed-tilt (固定轴) PV generation model
for both hemispheres, ported from the Excel workbook
《南北半球光伏固定轴发电模型_v1.0》 (8760-hourly, NSRDB weather input):

- Solar position (Spencer declination, equation of time, true solar time)
- Plane-of-array transposition: beam with the workbook's Rb backside rule;
  diffuse selectable between PVsyst (default), Hay&Davis, Reindl, and
  isotropic+circumsolar; ground reflection with albedo 0.2
- ASHRAE beam IAM (AOI > 87° → 0)
- Sandia array temperature model with the six mounting coefficient sets from
  the workbook
- First-year linear degradation (2 %/8759 per hour), 0.5 %/yr long-term,
  combined losses (−8 % default), 0.98 inverter efficiency

```python
from twoaxistracking.fixed_axis import simulate_fixed_axis, lookup_city

result = simulate_fixed_axis(
    ghi, dhi, dni, wind_speed, temp_air, dayofyear, hour,
    latitude=24.54, longitude=112.11, tilt=10., surface_azimuth=260.,
    capacity_kwp=146.76, tz_meridian=120.)
result.energy        # annual energy (kWh)
result.effective_hours  # equivalent full-power hours (h/day)
```

`twoaxistracking/data/china_pv_lookup.csv`（`lookup_city()` 读取）is the
《中国各地区光伏最佳安装倾角及峰值日照速查表》: recommended tilt, peak sun
hours, and first-year yield per watt for ~360 Chinese cities.

The model is validated cell-by-cell against the workbook's cached results
(see `tests/test_fixed_axis.py`, agreement within 1e-6).

## Documentation
The documentation can be found at [readthedocs](https://twoaxistracking.readthedocs.io/).

## Installation and dependencies
The package can be installed using pip:

    pip install twoaxistracking

## Citing
If you use the package in published work, please cite:
> Adam R. Jensen, Ioannis Sifnaios, and Kevin Anderson. "twoaxistracking – a python package for simulating
> self-shading of two-axis tracking solar collectors." MethodsX, 9, 101876, (2022).
> [https://doi.org/10.1016/j.mex.2022.101876](https://doi.org/10.1016/j.mex.2022.101876)

and

> Adam R. Jensen, Ioannis Sifnaios, Simon Furbo, and Janne Dragsted. "Self-shading of two-axis
> tracking solar collectors: Impact of field layout, latitude, and aperture shape." Solar
> Energy, 236, 215–224, (2022). [https://doi.org/10.1016/j.solener.2022.02.023](https://doi.org/10.1016/j.solener.2022.02.023)

## Contributing
Contributions to the repository, e.g., bug fixes and improvements to speed up the code are more than welcome.

## License
BSD 3-clause.
