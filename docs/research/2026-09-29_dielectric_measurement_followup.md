# 介电色散实测文献跟进（2026-09-29）

状态：调研记录存档；**当前批次已是色散材料（t3 批量 99/99 例使用 Debye 色散覆盖层），本调研结论暂不改动任何冻结配置**。后续若做场地标定或低频锚点复核时再启用。

## 背景与目的

当前仿真配置（`configs/research/dispersion_materials_v0.1.json`，2026-09-28 冻结）的 7 个频点锚点全部来自文献综合（West 2003 / Ishida 2000 / González-Teruel 2022），非场地实测；其中 ε'(20 MHz)=28 为最不确锚点，且 gprMax 无 Cole-Cole 导致 Debye 拟合在低频端有已知残差（拟合 22.8 vs 锚点 28）。本轮检索 2023–2026 年土壤/岩石复介电常数与电导率**实测**文献，评估是否有更贴 20–170 MHz 频段的实测数据可替换或校验现有锚点。

检索原始数据（Scholar 引擎 CSV，存于工作区根目录，不入库）：
- `E:\automation_djh\scholar_dielectric_soil.csv`（soil complex permittivity conductivity broadband measurement，7 条）
- `E:\automation_djh\scholar_dielectric_rock.csv`（dielectric properties rock GPR，2 条）
- `E:\automation_djh\scholar_em_landslide.csv`（electromagnetic parameter landslide soil water content GPR，12 条）

## 有用结果（按对本项目的用途排序）

1. **Zhang et al. 2025, Measurement**——微带线测量系统实测 20–400 MHz 土壤介电常数与电导率。频段完整覆盖本项目 20–170 MHz，比 West 2003（75 MHz 起测）更贴低频端；若复核 20 MHz 锚点，此为主要候选实测依据。<https://www.sciencedirect.com/science/article/pii/S0263224125010358>
2. **Hakiki & Lin 2025, ACS Measurement Science Au**（开放获取）——宽带复介电谱的 Cole–Cole vs 电路模型对比；可量化"Debye 无展宽"残差对 ε'(f)/σ(f) 的实际影响，为 gprMax 缺 Cole–Cole 的已知限制提供误差标尺。<https://pmc.ncbi.nlm.nih.gov/articles/PMC12532061/>
3. **Bobrov et al. 2023, IEEE Trans.**（被引 19）——土壤复介电谱 10 kHz–8 GHz 宽带弛豫分解（Part II），可检验现有 24.65 MHz 弛豫峰取值合理性。<https://ieeexplore.ieee.org/abstract/document/10354377/>
4. **Stellini et al. 2023, Sensors**（开放获取，被引 13）——土壤宽带复介电常数测量方法汇总，含损耗因子与离子电导讨论。<https://www.mdpi.com/1424-8220/23/11/5357>
5. **Budzeń et al. 2024, Measurement**（被引 12）——压实度对复介电谱测量的影响；对应滑坡堆积体松散-压实差异，是过渡带参数（当前为半强度极点的设计假设）未来找实测依据时的方向。<https://www.sciencedirect.com/science/article/pii/S0263224124011928>
6. **Hasar et al. 2023, IEEE Trans. Microwave Theory Tech.**（被引 24）——同轴夹具土壤介电常数提取的简化校准；对"作者夹具谱→材料谱反演"（R13-D 剩余缺口）有方法参考价值。<https://ieeexplore.ieee.org/abstract/document/10005082/>
7. **Lin, Hakiki & Lin 2026, SSSAJ**——高含水率颗粒介质单探针宽带复介电谱（本项目覆盖层为含水粉质粘土，相关性高）。<https://acsess.onlinelibrary.wiley.com/doi/abs/10.1002/saj2.70298>

## 与仓库既有资料的关系

- 不重复：R11（González-Teruel 2020, Sensors）与 R12（González-Teruel 2022）已在 `docs/research/2026-09-24_cycle02_materials_and_source.md` 登记；本轮结果为其补充频段更近、更新（2023–2026）的实测文献。
- R13-D（Schmidt et al. Zenodo 15 个 S2P）仍未反演为材料参数，本轮文献不替代该工作；若未来反演，Hasar et al. 2023（条目 6）的夹具校准思路可参考。
- Peplinski 类混合模型本轮同样未见 20–170 MHz 带内实测依据，维持 AGENTS.md 长期记录"不使用"的决定不变。

## 决定

- 当前 t3 批量（`configs/research/batch2d_slope_t3_co/`，99/99 例）已使用 Debye 色散覆盖层（ε_∞=18.017, Δε=7.878, τ=6.4567 ns, σ_dc=0.003）+ 非色散基岩（9/0.001），色散对照批已量化其影响（界面回波 +2.5 ns、-3.3~-6.6 dB）。
- 色散参数 v0.1 的已知限制（20 MHz 锚点不确定、Debye 残差、过渡带无实测）在本轮调研中没有发现需要立即修改配置的实测证据；**维持冻结，继续主线工作**。
- 若后续出现场地钻孔/实验室标定，或需要低频端复核，优先启用条目 1（Zhang 2025）与 R13-D 反演。
