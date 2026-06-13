# 项目文件结构说明

## 核心类层次
```
Simulator (基类，在 FluidSimulator.h/cpp)
├── APICSimulator     (在 APICSimulator.h/cpp)
├── CGSimulator       (在 CGSimulator.h/cpp)
├── SubgridSimulator  (在 SubgridSimulator.h/cpp)
├── FreeSurfaceSeparationSimulator (在 FreeSurfaceSeparationSimulator.h/cpp)
├── VariationalCoupledSimulator (在 VariationalCoupledSimulator.h/cpp)
└── VariationalSimulator (在 VariationalSimulator.h/cpp)

RigidBody (独立结构体，在 RigidBody.h)
MeshSDF   (不规则三角网格 SDF，在 MeshSDF.h/cpp)

CaseCoupled                    (双向流固耦合)
CaseSubgrid                    (子网格精度流动)
CaseFreeSurfaceSeparation      (自由液面自然分离)
CaseFluidRigid_zly             (多密度刚体流固耦合)
CaseVariational                (旧版变分求解器对照)

assets/models/suzanne.obj      (不规则固体演示模型)
assets/shaders/fluid.vert/frag (流体粒子渲染与可视化)
assets/shaders/flat.vert/frag  (刚体与网格模型渲染)
```

---

## 各文件作用
| 文件 | 作用 |
|------|------|
| **FluidSimulator.h/cpp** | 流体模拟基类，实现 FLIP/PIC 基础算法（P2G/G2P传输、压力求解、粒子碰撞），同时提供刚体 SDF 边界处理、粒子颜色模式和固体网格重建接口 |
| **APICSimulator.h/cpp** | 继承 Simulator，实现 APIC 方法（比普通 PIC/FLIP 更好地保留局部仿射速度与旋转细节） |
| **CGSimulator.h/cpp** | 继承 Simulator，用共轭梯度法替代 SOR 迭代求解压力，并加入不规则固体边界下的稳定性保护 |
| **SubgridSimulator.h/cpp** | Batty 论文的 MAC 面流体体积分数加权压力投影，保留小于单个网格的通道流动 |
| **FreeSurfaceSeparationSimulator.h/cpp** | Batty 论文第 4 节壁面自然分离条件，使用 active-set 求解 KKT 互补条件 |
| **VariationalCoupledSimulator.h/cpp** | 将切割面体积分数、自然分离约束和刚体反馈接入双向耦合 |
| **VariationalSimulator.h/cpp** | 变分流固耦合：最小化系统动能，满足不可压与接触约束 |
| **MeshSDF.h/cpp** | 不规则三角网格 SDF。使用 tinyobjloader 读取 OBJ，计算点到三角形距离和内外符号，并构建 32³ 距离场缓存以加速查询 |
| **RigidBody.h** | 刚体定义（位置、速度、旋转、惯量、SDF 距离查询），支持 box、sphere、bunny-like 解析形状和 MeshSDF 不规则网格形状 |
| **CaseCoupled.h/cpp** | 双向流固耦合的 UI、渲染和主循环，支持 Box/Sphere/suzanne 形状切换、APIC/CG 切换、压强可视化和鼠标外力交互 |
| **CaseSubgrid.h/cpp** | 子网格窄通道独立演示，可切换体积分数权重与二值权重 |
| **CaseFreeSurfaceSeparation.h/cpp** | 自由液面离墙独立演示，可查看 KKT 残差和分离状态 |
| **CaseFluidRigid_zly.h/cpp** | 多种密度刚体在流体中运动的独立演示，支持密度、半径等参数调节 |
| **assets/models/suzanne.obj** | 不规则固体耦合使用的三角网格模型，相比高精度 bunny 更适合当前实时 SDF 与低分辨率网格演示 |
| **assets/shaders/fluid.vert/frag** | 流体粒子渲染 shader，支持独立相机矩阵和速度、密度、压强等颜色模式 |
| **assets/shaders/flat.vert/frag** | 刚体、边界和 suzanne 网格渲染 shader，显式使用 model/view/projection 矩阵 |
| **tests/SimulatorTests.cpp** | 子网格、自然分离、MeshSDF 和双向耦合的数值回归测试 |

---

## 子网格精度流动
`Sub-grid Flow` 演示使用论文第 2.1 节公式 (4)-(7)。每个 MAC
速度自由度不再使用 0/1 体素权重，而是通过固体 SDF 子采样估计其控制体内
的流体体积分数，并以该分数构造加权压力矩阵 `G^T M_F G`。
演示中可通过 `Sub-grid weights` 开关比较体积分数权重和二值体素权重。半网格窄通道的定量对照保留在 `SimulatorTests.cpp` 中。

## 自由液面自然分离

`Free-surface Separation` 演示的 `Wall separation` 对应论文第 4 节。标准固壁条件会把壁面
法向速度固定为零，使自由液面无法自然离墙。开启 `Wall Separation`
后，对所有有液体支撑的流体-固体边界施加非负边界压力约束。active-set 对每个
候选壁面在 `p = 0` 的分离状态和无穿透接触状态之间切换，并检查
`u dot n >= v_solid dot n` 以及接触冲量方向，允许流体离墙但禁止穿墙。单元测试覆盖负压释放、非穿透、非负接触压力和动态 active-set 收敛。

耦合求解器对部分占据的 MAC 控制体同时保留开放流体通量和刚体边界通量，并根据刚体 SDF 重建切割面到压力点的距离。候选固壁上的压力作为
显式变量加入凸 QP，并施加论文公式 (15) 的 `p >= 0` 约束；约化后的 LCP 使用 active-set 求解，并以完整 KKT 残差验收。刚体矩阵项保持论文
公式 (13) 的 rank-6 形式，通过 Woodbury 恒等式加入，避免生成稠密的 `J^T M_S^-1 J` 块。速度、压力和刚体反馈均来自同一个约束解。

## 基于网格 SDF 的不规则刚体耦合

不规则固体演示对应 `CaseCoupled` 中的 `suzanne` 模式。该模式使用 `assets/models/suzanne.obj` 作为三角网格刚体，由 `MeshSDF` 将 OBJ 模型转换为可快速查询的有符号距离场，再通过 `RigidBody` 的统一 SDF 接口接入流体边界、粒子碰撞、压力投影和刚体受力反馈。

`MeshSDF::LoadOBJ` 首先读取 suzanne 的三角面片，并将模型中心化、归一化到局部坐标系。符号距离的大小由查询点到三角形的最近距离给出，符号由射线相交奇偶性判断。为了避免实时模拟中每个粒子或网格单元都遍历全部三角形，初始化阶段会调用 `BuildDistanceGrid` 构建 32³ 缓存距离场；运行时 `SignedDistance` 主要通过三线性插值返回近似 SDF，从而显著降低不规则边界查询开销。

在耦合过程中，流体网格单元通过 `RigidBody::GetSDF` 判断是否位于 suzanne 内部或接触面附近，粒子进入网格内部时会沿 SDF 法线推出，并修正法向速度。压力求解结束后，接触区域压力沿 SDF 法线累计为刚体受力，再通过力臂计算力矩，从而使 suzanne 能受到流体压力与浮力影响。APIC 模式复用同一 SDF 边界以改进粒子-网格速度转移；CG 模式则在复杂网格边界下加入孤立分量钉扎、预条件求解、非有限值回退和速度限制，以提高 suzanne 模式的稳定性。
