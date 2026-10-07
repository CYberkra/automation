# 后处理续查：剩余条带的原始时间来源

日期：2026-10-08；源码基线d42a097。用户要求继续本地排查。只复用428道native与已审计复谱缓存，CPU重建，没有FDTD、训练、远程任务、材料/几何/原包改动。前序：[后处理错误与基线](2026-10-08_line9_postprocessing_audit.md)。

## 新结论

1. **原Gaussian深部背景与规范Hann下的剩余晚时响应，必须分开解释。** 原Gaussian中前120ns的带限泄漏占据深部参考窗；Hann显著压低它后，剩余信号主要来自实际原始记录的240–400ns，不能再全部称为早时旁瓣。
2. **1m包的低频“深部增强”，还有浅时段响应展宽的影响。** 即使去掉前120ns，底砂参考窗内的20–39.8MHz剩余形态，仍与只取120–240ns原始场重建的形态高度相关，约0.95–0.97；不能据此认定深部砂岩更强。
3. **修正处理链能够保留人为注入的倾斜回波。** 在实际三包复谱上加入−70/−80/−90dB已知曲线，正/反相配对差分均保留正确形态，最大峰位误差≤0.411ns。没有发现将所有倾斜界面系统拉平或大幅移位的逆变换错误。

上述没有认证真实砂岩检出。实际240–400ns场可能包括目标、其他地下/侧向路径、多次反射和数值误差；本地后处理不能把这些物理身份自动分离。

## 一、将原始时间轴分为五段

原始段为0–120、120–240、240–400、400–600、600ns至记录末。交界用10ns与20ns两种平滑重叠，各段权重逐样点相加严格为1。所有分量分别做实际源归一化及精确501点复变换，随后用相同Gaussian或Hann重建；复谱及复时间结果重新相加验证闭合。

这不是地质消融：分段会改变频谱，段边缘会产生新的带限结构，所以不把段当作某一材料的纯回波。两档交界一致性用于检查结果是否严重依赖边缘选择。复谱/复时间闭合相对L2≤7.0e−15。

在实际体素底部砂岩的95MHz近似到时±10ns窗内，20ns交界结果：

|包|Gaussian全场与仅0–120ns分量相关|Hann全场与仅240–400ns分量相关|Hann该窗全场RMS，相对固定Gaussian峰值|
|---|---:|---:|---:|
|pkg1|0.99984|0.97472|−103.04dB|
|pkg2|0.98072|0.99960|−81.13dB|
|pkg3|0.99189|0.99829|−83.54dB|

10ns交界结论相近。**复数相关绝对值不是能量占比。** 该表说明不同频窗下深部显示的主要时间来源改变，不证明240–400ns分量等于砂岩反射。

Hann下，120–240ns分量在底砂参考窗的RMS约−147.29/−125.86/−118.92dB，分别远低于240–400ns分量约−103.22/−81.13/−83.53dB。这支持全频Hann的剩余晚时信号不是被浅段的宽波包完全解释。

对“深部都是早时泄漏”的表述作明确限制：只适用于已核验的原Gaussian背景及指定窗，不能不加区分地用于修正后所有频窗、所有剩余成分。

## 二、低频图中还有浅时段展宽

保持原扫频网格，实际选20–39.8MHz共67点。分成0–120、120–240、240ns之后三组，按源归一化，频点数归一化正确，不做AGC/叠道。每组复谱重新相加与审计缓存一致，相对L2≤2.5e−14。

在每道底砂近似到时±10ns窗，先去掉0–120ns分量，再将剩余全场与仅120–240ns分量比较：

|包|Gaussian相关：10ns/20ns交界|Hann相关：10ns/20ns交界|
|---|---:|---:|
|pkg1|0.255/0.265|0.362/0.362|
|pkg2|0.396/0.406|0.414/0.420|
|pkg3|**0.954/0.969**|**0.948/0.965**|

该现象主要出现在1m包，不能推广到另两包全部低频响应。1m实际底砂近似到时为258–290ns，较早120–240ns原始响应包含覆盖层底、浅互层及其他路径；窄带重建将它们展宽到后续时间。前序Gaussian理想单回波半幅全宽约96.8ns，与这些事件间隔同量级，混叠风险明显。

因此“去掉早时场后仍有低频深部亮带”也不是充分的目标证据；还需确认其来自哪段native时间。分段的谱变化仍是限制，10/20ns复核提升稳定性证据，不消除此限制。绝不把0.969写成“96.9%能量来自覆盖层”。

## 三、在实际数据上注入已知倾斜信号

为区分“处理链不能保形”和“原数据中目标身份/强度不明”，在实际审计复谱上增加已知复响应：

`H_inject(f,x)=a·exp(−i2πfτ(x))`

τ(x)取实际H5的底砂近似参考曲线，仅用于构造已知形状；a取原始全频Gaussian峰值的−70/−80/−90dB，各自生成正、反两种注入。该人为信号不是真实FDTD目标，不加入训练/验证数据，也不展示成恢复后的砂岩。

对每个包、幅度和Gaussian/Hann窗口：`process(H±H_inject)−process(H)`与单独处理的已知信号比较，相对L2最大≤1.8e−11；正反相符号保持，配对差分恢复的人为回波最大峰位误差≤0.411ns，低于0.832ns采样的一半。

这排除了修正线性链对上述已知曲线的系统性拉平/大延时错误。它**没有**证明混合图中每道都能盲检出注入信号，也没有证明真实场存在同强度同形状的底砂响应。尤其不能将配对差分的已知真值优势当成实际无背景处理的能力。

## 四、本地处理的当前可用结论

- 使用全频501点、保相位、正确幅度与时间轴的规范总场，Gaussian/Hann并列。不要再以低频图的亮带或AGC结果单独声称深部恢复。
- 保留Hann下的真实晚时记录响应，它们不是应被无差别删除的“早时噪声”。目前也不因其在深部窗就命名为砂岩。
- 不将时间分割设为默认生产流程：它需要完整native记录且可能伤及真实事件；本轮只用于查来源。
- 没有新发现一个修正后就能让三包恢复完整模型轮廓的后处理错误。已经证实的错误已由独立基线纠正；剩余是可见性和响应身份问题，不能继续把责任全推给IFFT或全推给航高。

## 产物与复现

全频原始来源图，仓库相对路径：

- [pkg1](../../artifacts/research_checks/2026-10-08_line9_native_time_origin_r1/line9_pkg1_v5_d02m_194st_native_time_origin.png)
- [pkg2](../../artifacts/research_checks/2026-10-08_line9_native_time_origin_r1/line9_pkg2_x160-110_agl35_142st_native_time_origin.png)
- [pkg3](../../artifacts/research_checks/2026-10-08_line9_native_time_origin_r1/line9_pkg3_x160-120_agl1_70st_native_time_origin.png)
- [真实1m数据上的人为回波注入图，须按行区分](../../artifacts/research_checks/2026-10-08_line9_native_time_origin_r1/pkg3_known_echo_injection_on_actual_data.png)

低频原始来源图：

- [pkg1](../../artifacts/research_checks/2026-10-08_line9_low_band_origin_r1/line9_pkg1_v5_d02m_194st_low_band_origin.png)
- [pkg2](../../artifacts/research_checks/2026-10-08_line9_low_band_origin_r1/line9_pkg2_x160-110_agl35_142st_low_band_origin.png)
- [pkg3](../../artifacts/research_checks/2026-10-08_line9_low_band_origin_r1/line9_pkg3_x160-120_agl1_70st_low_band_origin.png)

两批JSON保留输入审计/缓存/脚本哈希、交界宽度、闭合、复数相关及注入检验；所有原始资料保持私有。图件共享固定全频Gaussian参考，未逐面板归一化。

```powershell
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/trace_line9_time_origin.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --cache artifacts/local_checks/2026-10-08_line9_postprocessing_cache_r1 --out artifacts/local_checks/native_origin_reproduce
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/trace_line9_low_band_origin.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --cache artifacts/local_checks/2026-10-08_line9_postprocessing_cache_r1 --out artifacts/local_checks/low_band_origin_reproduce
```

输出需新目录，依赖V4 SFCW工具箱、NumPy、h5py、Matplotlib和中文字体。该分析复用前轮缓存，脚本对原始native身份及重算复谱再次核验，未调用求解器。阶段完成仅代表本地后处理诊断，不代表物理验收、独立实测验证或模型训练资格。
