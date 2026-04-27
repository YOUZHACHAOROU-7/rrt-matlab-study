# GitHub shortlist for paper-oriented manipulator trajectory planning

## 目标论文方向
建议优先走：
- 基于改进粒子群算法的六自由度机械臂时间最优轨迹规划
- 备选：基于改进灰狼算法的六自由度机械臂时间最优轨迹规划

原因：
1. 与现有中文核心常见写法一致：D-H建模 + 3-5-3多项式 + 时间优化。
2. GitHub上更容易找到干净的基线代码。
3. 后续只需轻量改动即可形成“改进xx算法”的论文结构。

## 已筛到的高价值公开仓库
### 1) 机械臂轨迹骨架
- mathworks-robotics/trajectory-planning-robot-manipulators
  - 价值：提供 joint space / task space 轨迹脚本、绘图工具、startupExample 入口。
  - 用途：作为机械臂轨迹生成与可视化骨架参考。

### 2) PSO基线
- MatthewPeterKelly/ParticleSwarmOptimization
  - 价值：MATLAB版PSO，结构清晰，包含测试脚本和迭代历史绘图。
  - 用途：作为标准PSO优化器骨架，后续可改惯性权重、学习因子、扰动策略。

### 3) GWO基线
- mzychlewicz/GWO
  - 价值：MATLAB版灰狼优化器，适合做改进GWO基线。
  - 用途：若最终选择“改进灰狼算法”方向，可直接作为起点。

## 真正需要从GitHub上搬的不是整篇论文代码，而是三块零件
1. 机械臂建模与轨迹绘图骨架
2. 标准优化器（PSO / GWO）
3. 适应度函数与约束处理

## 推荐拼接方式
### 路线A：PSO论文主线（首选）
- 机械臂：使用现有MATLAB建模/轨迹脚本
- 插值：自己写3-5-3分段多项式
- 优化器：基于PSO仓库改成IPSO
- 论文创新点建议：
  1. 非线性惯性权重
  2. 动态学习因子
  3. 柯西/高斯/量子扰动三选一
  4. 贪婪保优策略

### 路线B：GWO论文主线（备选）
- 机械臂：同上
- 插值：同上
- 优化器：基于GWO仓库改成IGWO
- 论文创新点建议：
  1. 余弦收敛因子
  2. 融合PSO个体最优信息
  3. 柯西变异
  4. 贪婪策略

## 当前不建议直接做的方向
- 多目标蛾群：多目标复杂度高，工作量大
- 视觉三维重建 + 焊接轨迹：过重，不适合快速出稿
- 完整避障RRT*主线：更偏路径规划，不是最省力的时间优化路线

## 下一步建议
1. 继续在本仓库单独分支上搭 `paper_psO_gwo_baseline` 结构
2. 写 `main_compare.m`：对比标准PSO / 改进PSO / 标准GWO
3. 写 `poly353_traj.m`：3-5-3插值
4. 写 `fitness_time_opt.m`：时间最优目标 + 速度/加速度约束
5. 生成收敛曲线、位移/速度/加速度图、结果表格
