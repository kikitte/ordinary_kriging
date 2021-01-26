## 普通克里金插值

一个使用ArcGIS临近样本点搜索策略的普通克里金插值C语言实现。

### 1. 亮点

- 使用ArcGIS中对临近样本点的搜索策略，在半变异模型及其参数相同的情况下可得到相同结果
- 多线程支持（默认使用全部可用的CPU），在得到相同插值结果的情况下速度为ArcGIS的4倍以上
- 参数通过JSON文件进行配置

### 2. 如何使用

在进行后续步骤之前请先确保安装GDAL依赖，由于不同操作系统的差异无法一一概括，所以此处不做介绍。

```bash
# 克隆源码
git clone https://www.github.com/kikitte/ordinary_kriging
# 编译
cd ordinary_kriging
make
# 运行测试数据
build/ordinary_kriging test/config.json
```

配置文件介绍（以默认提供的配置文件为例）：

```json
{
  // 输入的样本点JSON文件路径，该JSON内容基本格式为：Array<[x, y, z]>,可参见提供的事例数据
  "INPUT_POINTS": "/home/kikitte/geoanalyst/kriging/chap15/stations.json",
  // 输出的栅格文件路径
  "OUTPUT_RASTER": "/home/kikitte/geoanalyst/kriging/debug/out2.tif",
  // 半变异模型参数
  "VariogramModel": {
    // 半变异模型类型："Spherical", "Exponential", "Gaussian", "Wave", "RationalQuadratic", "Circular"
    "TYPE": "Spherical",
    // 块金
    "NUGGET": 7.179056886219,
    // 基台值
    "SILL": 77.604636335401,
    // 搜索范围
    "RANGE": 480000
  },
  // 临近样本点搜索策略：完全参考ArcGIS
  "NeighborhoodOption": {
    // 扇区类型："OneSector", "FourSector", "FourSector45Offset", "EightSector"
    "SECTOR_TYPE": "FourSector45Offset",
    // 每个扇区内最少样本点数目
    "MIN_NEIGHBORDS": 2,
    // 扇区内最多样本点数目
    "MAX_HEIGHBORDS": 5,
    // 扇区搜索最大距离
    "MAX_DISTANCE": 480000
  },
  // 输出栅格数据选项
  "RasterInfo": {
    // 栅格左上角点X坐标
    "TOPLEFT_X": 227942.9990586524945684,
    // 栅格左上角点Y坐标
    "TOPLEFT_Y": 881299.6237735910108313,
    // 栅格像元分辨率
    "RESOLUTION": 100,
    // 栅格列数
    "COLS": 5388,
    // 栅格行数
    "ROWS": 8128,
    // 栅格空值
    "NODATA_VALUE": 0,
    // 输出的栅格的格式，参考：https://gdal.org/drivers/raster/index.html
    "GDAL_DRIVER": "GTiff",
    // 栅格像元数值类型
    "GDAL_GDT": "Float32"
  }
}
```

### 3. 版权声明

无。该程序为业余时间所作，纯属个人行为。你可以任意修改使用，如果这个工作对你有帮助话那就说明它是有价值的。

### 4. 致谢

程序的完成参考或使用了很多其他人的代码，在此一并感谢，如下：

- [cJSON](https://github.com/DaveGamble/cJSON) 

  用于参数解析和插值样本数据读取。

- AX=B线性方程组求解

  [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode)

  [LU decompossion](https://en.wikipedia.org/wiki/LU_decomposition#C_code_example)

  [Don’t invert that matrix](https://www.johndcook.com/blog/2010/01/19/dont-invert-that-matrix/)

  [Triangular_matrix Forward_and_back_substitution](https://en.wikipedia.org/wiki/Triangular_matrix#Forward_and_back_substitution)

  事实上普通克里金插值实现中求解线性方程组没必要对矩阵完全求逆进而求解权重。使用高斯消元法然后使用反向就可以求得权重。

- 半变异模型

  参考了该仓库的半变异模型: https://github.com/cnca/Kriging

- Esri or ArcGIS

  我只是模仿者。