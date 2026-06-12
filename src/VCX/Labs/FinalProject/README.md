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

CaseCoupled                    (双向流固耦合)
CaseSubgrid                    (子网格精度流动)
CaseFreeSurfaceSeparation      (自由液面自然分离)
CaseVariational               (旧版变分求解器对照)
```

---

## 各文件作用

| 文件 | 作用 |
|------|------|
| **FluidSimulator.h/cpp** | 流体模拟基类，实现 FLIP/PIC 基础算法（P2G/G2P传输、压力求解、粒子碰撞） |
| **APICSimulator.h/cpp** | 继承 Simulator，实现 APIC 方法（比普通 PIC 更精确） |
| **CGSimulator.h/cpp** | 继承 Simulator，用共轭梯度法替代 SOR 迭代求解压力 |
| **SubgridSimulator.h/cpp** | Batty 论文的 MAC 面流体体积分数加权压力投影，保留小于单个网格的通道流动 |
| **FreeSurfaceSeparationSimulator.h/cpp** | Batty 论文第 4 节壁面自然分离条件，使用 active-set 求解 KKT 互补条件 |
| **VariationalCoupledSimulator.h/cpp** | 将切割面体积分数、自然分离约束和刚体反馈接入双向耦合 |
| **VariationalSimulator.h/cpp** | 变分流固耦合：最小化系统动能，满足不可压与接触约束 |
| **RigidBody.h** | 刚体定义（位置、速度、旋转、SDF 距离查询） |
| **CaseCoupled.h/cpp** | 双向流固耦合的 UI、渲染和主循环 |
| **CaseSubgrid.h/cpp** | 子网格窄通道独立演示，可切换体积分数权重与二值权重 |
| **CaseFreeSurfaceSeparation.h/cpp** | 自由液面离墙独立演示，可查看 KKT 残差 |
| **tests/SimulatorTests.cpp** | 子网格、自然分离和双向耦合的数值回归测试 |

---

## 子网格精度流动

`Sub-grid Flow` 演示使用论文第 2.1 节公式 (4)-(7)。每个 MAC
速度自由度不再使用 0/1 体素权重，而是通过固体 SDF 子采样估计其控制体内
的流体体积分数，并以该分数构造加权压力矩阵 `G^T M_F G`。

演示中可通过 `Sub-grid weights` 开关比较体积分数权重和二值体素
权重。半网格窄通道的定量对照保留在 `SimulatorTests.cpp` 中。

## 自由液面自然分离

`Free-surface Separation` 演示的 `Wall separation` 对应论文第 4 节。标准固壁条件会把壁面
法向速度固定为零，使自由液面无法自然离墙。开启 `Wall Separation`
后，仅把同时邻接固体和空气的流体格视为候选区域。active-set 对每个
候选壁面在 `p = 0` 的分离状态和无穿透接触状态之间切换，并检查
`u dot n >= v_solid dot n` 以及接触冲量方向，允许流体离墙但禁止穿墙。
单元测试覆盖负压释放、非穿透、非负接触压力和动态 active-set 收敛。

耦合求解器对部分占据的 MAC 控制体同时保留开放流体通量和刚体边界
通量，并根据刚体 SDF 重建切割面到压力点的距离。候选固壁上的压力作为
显式变量加入凸 QP，并施加论文公式 (15) 的 `p >= 0` 约束；约化后的
LCP 使用 active-set 求解，并以完整 KKT 残差验收。刚体矩阵项保持论文
公式 (13) 的 rank-6 形式，通过 Woodbury 恒等式加入，避免生成稠密的
`J^T M_S^-1 J` 块。速度、压力和刚体反馈均来自同一个约束解。
