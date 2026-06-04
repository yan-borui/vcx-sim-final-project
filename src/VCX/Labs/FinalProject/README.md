# 项目文件结构说明

## 核心类层次

```
Simulator (基类，在 FluidSimulator.h/cpp)
├── APICSimulator     (在 APICSimulator.h/cpp)
└── CGSimulator       (在 CGSimulator.h/cpp)

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
| **RigidBody.h** | 刚体定义（位置、速度、旋转、SDF 距离查询） |
| **CaseCoupled.h/cpp** | 主程序入口、UI界面、渲染循环、流体-刚体耦合的主逻辑 |

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
