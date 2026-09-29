# 跨领域同类思想调研（防走偏审查，2026-09-29）

用户指令："调研一下其他领域中我们这种类似的思想的应用，确保我们目前没有走偏。"
本笔记把本项目的方法论内核提炼成四条，逐条拿到成熟领域对表，最后给出走偏检查结论。
纯调研，不动任何契约与产物。

## 0. 方法内核提炼

- **K1 已知真值控制例**：不依赖"现场真值未知"的数据评价算法，而是自己施加闭式已知的
  损伤/衰减（损伤阶梯、增益控制例），以构造参考（原始未损伤事件）为质量标签。
- **K2 结果导向、算法无关的指标**：只测"对已知参考的恢复结果"（保真 D_e、幅度 a、
  负控 N_b），不测算法内部性质；亮度/能量/PSNR 类指标明确排除。
- **K3 硬门先于标量分**：可行性门（冻结 τ/ε）不过则整例 infeasible，通过才谈奖励标量。
- **K4 评价先于训练 + 防过拟合纪律**：奖励协议冻结在前、训练在后（G4 未解除）；
  阈值不从候选表现倒推（Cawley–Talbot）；预注册预期（F1–F6）；开发族/测试族 {C5,C8}
  分离；r1/r2 字节一致。

## 1. 医学影像（最接近的同构领域）

- **虚拟成像试验（VIT）/ 数字参考物体（DRO）**：AAPM TrueCT 重建大奖赛用计算假人
  （数字真值精确已知、解剖学真实）评价 CT 重建算法，指标分 task-generic（RMSE/SSIM）
  与 **task-specific（可检测指数 d'、病灶体积精度）**——task-specific 对应我们的
  "事件恢复对已知参考"，task-generic 的 RMSE/SSIM 恰是他们放在次要层的指标
  （对应我们排除能量/亮度类）。AAPM TG387（2024）已有 VIT 可靠开发与使用的共识
  建议；VIC-TRE 证明计算 phantom 可进入**监管评价**。
- **混合注入（hybrid insertion）**：把合成病灶插入真实临床影像再跑同一处理链，
  是 lesion-insertion 类评价的标准做法——**与我们的做法同构**：衰减施加在实测锚点
  签名图（官方 SFCW 链产物）上，而非重新仿真，恰是 hybrid 而非纯仿真。
- **幻影生成机制的算法偏置教训（重要镜像）**：DIR（可变形配准）评审文献明确警告：
  用 B-spline 生成的形变 phantom 会系统性偏置 B-spline 类配准算法（Fatyga 2015 等）。
  镜像到我们：增益效果表 R2 读数——单调上升曲线与施加衰减结构（阶跃/凸起）形状
  错配、只有契约导出 oracle 能精确适配——**说明我们的控制例施加机制（闭式乘性）与
  受测算法不同族，这是优点**；但警戒线已立：未来若把衰减结构设计成某个算法的可逆
  形状，就落入同一偏置。
- **仿真不外推临床（VIT 文献的一致告诫）**：VCT 结论须向 predicate 技术标定、
  ASME V&V 40 按用途评定仿真可信度。对应我们已有的防线：constructed_reference
  标注、色散实测调研先行、设备天线传递函数外置、截幅结论不外推设备。

## 2. 引力波（端到端验证思想）

- **硬件注入（hardware injections, LIGO）**：用执行器物理位移测试质量、注入已知
  波形，对**分析管线做端到端测试**；每次探测事件后注入相似波形复现以验证搜索与
  参数估计；还用于校准符号/时延、检查"安全通道"（注入信号不得出现在无关通道——
  **与我们的 NC 负控（恢复时不得放大负控窗）同构**）；早期还有 blind injections
  （仅小团队知情的盲注）。
- 对应关系：我们的恒等锚/oracle 锚（a≡k、a=1 精确自洽）= 注入-恢复一致性检查；
  冻结契约 + SHA 断言 = 注入波形的版本锁定。

## 3. 天文学巡天（注入-恢复完备性）

- **合成源注入（SSI/Balrog, DES Y6）**：把已知属性的源注入真实巡天图像、用
  **与真实数据完全相同的软件版本**处理注入图像，测 transfer function 与完备性——
  "同一加载链"纪律的独立先例。
- **FRB/MeerKAT 等的 injection-recovery**：注入合成脉冲定 S/N 阈值与完备率曲线
  （如 HEIMDALL 阈值 10 的设定依据）——**阶梯刻度→阈值提案的同款路线**，且阈值
  从注入-恢复数据定、不从候选算法表现定，与我们一致。

## 4. ML 评测方法论（防奖励攻击）

- **CheckList / 变形测试（metamorphic testing）**：解决"无 oracle"评测。三类测试：
  MFT（已知正确输出的小测试集）= 我们的损伤/增益阶梯；INV（不该改变输出的扰动）=
  我们的 nc_untouched 负控；DIR（有已知期望方向的扰动）= 我们的预注册恢复预期
  （rec_dB、nc_ratio 差 k²）。
- **奖励攻击/Goodhart（RL/RLHF 文献）**：代理指标一旦被优化就会失真（Skalse 2022
  形式化；Krakovna 2020 目录；ROUGE 摘要游戏=proxy metric gaming 的范例）。
  **我们的 M7/M8/M9 行（rec 52–60 dB 但 a 门拒）就是一次实测的 proxy gaming
  拦截**。2026 年 AG-RLHF 的 auditor-gated rewards 证明"硬门当审计员"可降低
  过优化——为 K3 提供了独立背书。SpecBench 用验证/隐藏测试分离测 hacking
  （Δ=s_val−s_test）——**与我们的开发族/测试族 {C5,C8} 分离同构**。

## 5. 走偏检查清单（逐条判定）

| 内核 | 领域先例 | 我们当前做法 | 判定 |
|---|---|---|---|
| K1 已知真值控制例 | VIT/TrueCT、hybrid insertion、LIGO 硬件注入、SSI | 损伤/增益阶梯 + 契约导出 oracle，constructed_reference 全程标注 | **没走偏**，且属主流做法 |
| K2 结果导向指标 | task-based IQA（d'）、metamorphic DIR | D_e/a/N_b 对原始参考；能量/亮度排除 | **没走偏**；RMSE/SSIM 在 TrueCT 也居次席，与我们排序一致 |
| K3 硬门先于标量分 | accreditation pass/fail、auditor-gated rewards | 冻结 τ/ε 硬门 → R_bg/R_gain 标量 | **没走偏**；效果表证明门会拦"视觉好"的算子 |
| K4 评价先于训练+防过拟合 | SpecBench val/test 分离、预注册、VCT 标定 | G4 未解除、阈值不倒推、{C5,C8} 分离、F1–F6 预注册 | **没走偏**；但注意 Goodhart 是动态威胁（见 §6 R3） |

## 6. 风险登记（其他领域的教训映射回我们）

- **R1 仿真-真实差距**（VIT/V&V40 的最大告诫）：我们的防线已部分存在（色散实测调研、
  天线传递函数外置、不外推设备声明），但**"仿真可信度按用途评定"还没有成文**——
  建议在 START_HERE 级别补一条：本仓库证据链的可信度边界随场景族扩展逐步评定。
- **R2 控制例生成机制的算法偏置**（B-spline phantom 镜像）：当前施加机制与算法族
  无关，安全；未来设计新结构时须自问"这是否恰好在为某类算法量身定制"。
- **R3 Goodhart 是动态威胁**：奖励冻结在前只解决"设计时"过拟合；G4 解除开始训练后，
  被优化的指标会暴露新攻击面。领域做法是**训练后复测**（SpecBench 的 Δ 审计）。
  建议：解除 G4 后的首次训练完成时，重跑阶梯+效果表+敏感性全套，确认排名未对冻结
  奖励过拟合，再进测试族。
- **R4 同一加载链纪律**：Balrog 先例强化了我们"官方 SFCW 链 + Hann 窗"的选择正确性；
  任何换链（如换窗函数）必须视为新证据链，不能与旧记录混比。

## 7. 来源

- AAPM TrueCT Reconstruction Grand Challenge（PMC11973969）与 aapm.org 赛事页；
  AAPM TG387 共识（2024，经 TrueCT 文引）；VIC-TRE（Badano 2017/2018，经 OpenReview
  综述引）；virtual clinical trials 综述（White Rose eprints）；
  Nenoff 2023 DIR 不确定性评审（B-spline 偏置，boris.unibe.ch）。
- Abbott et al., Advanced LIGO hardware injection system（arXiv:1612.07864）；
  LIGO S2 burst upper limits（arXiv:gr-qc/0505029，注入用于 veto 安全性）。
- DES Y6 Balrog（arXiv:2501.05683）；EFTE injection-recovery（arXiv:2011.02495）；
  Nanshan FRB 完备性（RAA 2022）；MeerKAT cube 注入（ApJ 2025）。
- Ribeiro et al. CheckList（经 arXiv:2504.18827 转述）；metamorphic testing 综述
  （HKU TSE；LLMorph arXiv:2603.23611）；Skalse 2022/Krakovna 2020/Goodhart
  （经 arXiv:2605.21384、2607.18064 转述）；AG-RLHF auditor-gated rewards
  （arXiv:2602.01750）；SpecBench（arXiv:2605.21384）。
