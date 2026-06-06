# 项目文件结构说明

## 核心类层次

```
Simulator (基类，在 FluidSimulator.h/cpp)
├── APICSimulator     (在 APICSimulator.h/cpp)
├── CGSimulator       (在 CGSimulator.h/cpp)
├── SubgridSimulator  (在 SubgridSimulator.h/cpp)
└── FreeSurfaceSeparationSimulator (在 FreeSurfaceSeparationSimulator.h/cpp)

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
| **FreeSurfaceSeparationSimulator.h/cpp** | Batty 论文第 4 节壁面自然分离条件，使用非负压力 active-set 近似 KKT 约束 |
| **RigidBody.h** | 刚体定义（位置、速度、旋转、SDF 距离查询） |
| **CaseCoupled.h/cpp** | 主程序入口、UI界面、渲染循环、流体-刚体耦合的主逻辑 |
| **CaseSubgrid.h/cpp** | 子网格精度流动 case 的独立 UI、场景和渲染逻辑 |
| **CaseFreeSurfaceSeparation.h/cpp** | 水团撞击左壁并自然剥离的独立演示 case，可切换标准壁面条件做对照 |

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

## 子网格精度流动 Case

`Sub-grid Accurate Flow` 使用论文第 2.1 节公式 (4)-(7)。每个 MAC
速度自由度不再使用 0/1 体素权重，而是通过固体 SDF 子采样估计其控制体内
的流体体积分数，并以该分数构造加权压力矩阵 `G^T M_F G`。

场景中的箱体与通道上下壁各只留 `0.5h` 的间隙。默认的子网格权重会保留
这些部分开放的 MAC 面；关闭 `Sub-grid weights` 后，权重退化为
二值体素，窄通道会被封闭，便于直接比较。

## 自由液面自然分离 Case

`Free-Surface Wall Separation` 对应论文第 4 节。标准固壁条件会把壁面
法向速度固定为零，可能通过负压把自由液面吸附在墙上。开启
`Wall Separation` 后，仅对同时邻接固体和空气的流体格施加
`0 <= p` 约束；active-set 将需要负压的格子改为 `p = 0` 自由表面，
并执行 `u dot n >= v_solid dot n`，允许流体离墙但禁止穿墙。演示采用
接近论文 Figure 6 的薄片水团，默认将靠近左墙的粒子标成橙色，并显示
贴墙粒子数和平均离墙距离，便于与关闭分离条件后的吸附结果对照。
