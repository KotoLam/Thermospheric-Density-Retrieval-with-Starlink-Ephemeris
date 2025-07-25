# 基于星链卫星预报星历数据反映热层大气密度变化
# Thermosphere-Density-Variations-Obtained-from-Starlink-Ephemeris
从 **SpaceX** 在 [space-track.org](https://www.space-track.org/#publicFiles)（具体下载API可查阅网站教程）上发布的星链卫星预报星历中获取超过6000颗星链卫星一分钟分辨率的位置和速度信息，并从中计算每颗卫星的机械能衰减量作为热层大气密度变化状况的代理指标。
# 本项目在以下项目的基础上实现
1. 由 Meysam Mahooti 提供的在matlab语言中实现的高精度轨道外推模型 High Precision Orbit Propagator
https://www.mathworks.com/matlabcentral/fileexchange/55167-high-precision-orbit-propagator
2. IAU 提供的标准天文运算库 Standards of Fundamental Astronomy
https://www.iausofa.org/2023_1011_C.html


# TODO
- 将所有matlab代码逐步变为C/C++以脱离MATLAB运行时环境（MCR）

