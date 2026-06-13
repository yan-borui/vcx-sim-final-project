# Batty 变分流固耦合仿真

本项目实现并演示 Batty、Bertails 与 Bridson 在论文
*A Fast Variational Framework for Accurate Solid-Fluid Coupling*
中提出的变分压力投影方法。程序基于三维 MAC 网格和粒子流体，重点展示：

- 流体与刚体的双向压力耦合；
- 粗网格上的变分压力投影；
- 小于一个网格宽度的子网格流动；
- 自由液面从固体壁面自然分离；
- 不规则网格固体和多种密度刚体。

## 构建与运行

环境要求：

- 支持 C++20 的编译器；
- [xmake](https://xmake.io/) 2.6.9 或更高版本；
- Windows、macOS 或 Linux 的 OpenGL 开发环境。

首次构建时 xmake 会自动下载 Eigen、GLM、GLFW、ImGui、tinyobjloader
等依赖。

```powershell
xmake f -m release
xmake build lab4
xmake run lab4
```

## 演示场景

程序左侧列表包含以下场景：

| 场景 | 内容 |
| --- | --- |
| `Coupled Simulation` | 统一变分求解器。支持 Box、Sphere 和 `suzanne.obj` 不规则刚体，显示压力残差与壁面分离 KKT 残差 |
| `Variational Simulation` | 变分耦合核心的独立粗网格演示，可调时间步长、FLIP 比例和压力迭代参数 |
| `Sub-grid Accurate Flow` | 窄通道流动，可切换体积分数权重和二值体素权重进行对照 |
| `Free-Surface Wall Separation` | 水团撞击左壁，可切换传统无穿透边界与自然分离条件 |
| `Fluid-Rigid Coupling 3D` | 多个不同密度球形刚体的三维流固耦合演示 |

所有场景均提供开始、暂停或重置控制。使用鼠标拖动和滚轮控制相机。
在 `Coupled Simulation` 中按住 `Ctrl` 并拖动鼠标左键可向刚体施加交互速度。

## 论文算法对应

### 变分压力投影

流体速度存储在 MAC 面中心。压力通过最小化压力更新后的系统动能求得，
离散矩阵对应论文公式 (6)：

```text
G^T M_F G p = G^T M_F u*
```

`M_F` 使用每个 MAC 速度控制体中的流体质量。自由液面采用 ghost-fluid
距离缩放。P2G 后，对有流体支撑但质量为零的 MAC 样本执行分层速度外推，
避免 G2P 插值将自由液面和固壁附近的切向速度错误地吸向零。

### 子网格精度

`SubgridSimulator` 和 `VariationalCoupledSimulator` 对 MAC 控制体进行 SDF
子采样，分别估计开放体积分数和流体质量分数。压力矩阵使用连续质量权重，
而不是把部分占据网格二值化，因此窄于一个网格的通道仍可传递流量。

### 双向刚体耦合

统一耦合求解器按论文公式 (8)-(13) 累积刚体边界的力和力矩。动态刚体项
保持 `J^T M_S^-1 J` 的 rank-6 结构，并通过 Woodbury 公式加入压力求解。
压力更新后的流体速度和刚体速度来自同一个变分系统。

### 自由液面自然分离

壁面压力变量满足论文公式 (15)：

```text
p >= 0
u.n - v_solid.n >= 0
p * (u.n - v_solid.n) = 0
```

接触状态要求非负压力并满足无穿透；分离状态令边界压力为零。代码使用
active-set 或约化非负二次规划求解，并报告压力残差和 KKT 残差。

### 不规则固体

`MeshSDF` 从 `assets/models/suzanne.obj` 构建缓存的有符号距离场。
同一 SDF 用于子网格体积分数、粒子碰撞、压力边界和刚体反馈。

## 代码结构

核心代码位于 `src/VCX/Labs/FinalProject`：

| 文件 | 作用 |
| --- | --- |
| `FluidSimulator.*` | 粒子积分、FLIP/PIC 传输、MAC 速度外推和基础压力投影 |
| `VariationalCoupledSimulator.*` | 统一变分流固耦合、rank-6 刚体项和壁面分离 QP |
| `SubgridSimulator.*` | 子网格质量权重独立演示 |
| `FreeSurfaceSeparationSimulator.*` | 自由液面壁面分离独立演示 |
| `RigidBody.h` | 刚体状态、惯量、解析形状与网格 SDF 接口 |
| `MeshSDF.*` | OBJ 网格有符号距离场 |
| `Case*.cpp` | 场景 UI、交互和渲染 |

更详细的文件说明见
[FinalProject README](src/VCX/Labs/FinalProject/README.md)。

## 参考文献

Christopher Batty, Florence Bertails, and Robert Bridson.
*A Fast Variational Framework for Accurate Solid-Fluid Coupling*.
ACM Transactions on Graphics, 26(3), 2007.
