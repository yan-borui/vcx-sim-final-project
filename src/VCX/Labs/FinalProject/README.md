# 项目文件结构说明

## 核心类层次

```
Simulator (基类，在 FluidSimulator.h/cpp)
├── APICSimulator     (在 APICSimulator.h/cpp)
├── CGSimulator       (在 CGSimulator.h/cpp)
├── SubgridSimulator  (在 SubgridSimulator.h/cpp)
├── FreeSurfaceSeparationSimulator (在 FreeSurfaceSeparationSimulator.h/cpp)
└── VariationalCoupledSimulator (在 VariationalCoupledSimulator.h/cpp)

RigidBody (独立结构体，在 RigidBody.h)

CaseFluid (UI/主循环，在 CaseCoupled.h/cpp)
```

---

## 各文件作用

| 文件 | 作用 |
|------|------|
| **FluidSimulator.h/cpp** | 流体模拟基类，实现 FLIP/PIC 基础算法（P2G/G2P传输、压力求解、粒子碰撞） |
| **APICSimulator.h/cpp** | 继承 Simulator，实现 APIC 方法（比普通 PIC 更精确） |
| **CGSimulator.h/cpp** | 继承 Simulator，用共轭梯度法替代 SOR 迭代求解压力 |
| **SubgridSimulator.h/cpp** | Batty 论文的 MAC 面流体体积分数加权压力投影，保留小于单个网格的通道流动 |
| **FreeSurfaceSeparationSimulator.h/cpp** | Batty 论文第 4 节壁面自然分离条件，使用逐壁面 active-set 近似 KKT 互补条件 |
| **VariationalCoupledSimulator.h/cpp** | 将子网格质量权重和逐壁面分离条件接入流固耦合 case |
| **RigidBody.h** | 刚体定义（位置、速度、旋转、SDF 距离查询） |
| **CaseCoupled.h/cpp** | 主程序入口、UI界面、渲染循环、流体-刚体耦合的主逻辑 |
| **tests/SimulatorTests.cpp** | 子网格、自然分离和双向耦合的数值回归测试 |

---


## 核心算法流程
```
void SimulateTimestep(float const dt) {
    float sdt = dt / numSubSteps;

    for (int step = 0; step < numSubSteps; step++) {
        integrateParticles(sdt);
        handleParticleCollisions();
        if (separateParticles)
            pushParticlesApart(numParticleIters);
        handleParticleCollisions();
        transferVelocities(true, m_fRatio);
        updateParticleDensity();
        solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
        transferVelocities(false, m_fRatio);
    }
    updateParticleColors();
}
```

## 子网格精度流动

Case 1 使用论文第 2.1 节公式 (4)-(7)。每个 MAC
速度自由度不再使用 0/1 体素权重，而是通过固体 SDF 子采样估计其控制体内
的流体体积分数，并以该分数构造加权压力矩阵 `G^T M_F G`。

Case 1 中可通过 `Sub-grid weights` 开关比较体积分数权重和二值体素
权重。半网格窄通道的定量对照保留在 `SimulatorTests.cpp` 中。

## 自由液面自然分离

Case 1 的 `Wall separation` 对应论文第 4 节。标准固壁条件会把壁面
法向速度固定为零，使自由液面无法自然离墙。开启 `Wall Separation`
后，仅把同时邻接固体和空气的流体格视为候选区域。active-set 对每个
候选壁面在 `p = 0` 的分离状态和无穿透接触状态之间切换，并检查
`u dot n >= v_solid dot n` 以及接触冲量方向，允许流体离墙但禁止穿墙。
单元测试覆盖负压释放、非穿透、非负接触压力和动态 active-set 收敛。
