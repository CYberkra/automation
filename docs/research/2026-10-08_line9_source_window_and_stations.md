# 源、时窗与另外两站：检验夹层归因能否复现

2026-10-08。上轮[远处夹层及同域航高对照](2026-10-08_line9_localized_height_controls.md)是一个经模型提示选择的站位，不能因此断言整线原因已唯一找到。本批继续用户“研究明白”和SSH自主研究授权，求解前冻结三项源/时窗与六项站位/材料控制，共9项，不扩成整线扫描。

## 冻结的科学比较

域110×40m，2.5cm，TMz线源，40A源幅度，1.3m收发间隔，固定y35m，HORIPML80格，原材料数据库及界面平均，V4.0.0原生FP64，全部无快照。源/时窗组复用原c0008/H0的冲激800ns，新增Ricker100MHz/800ns、冲激1600ns、Ricker100MHz/1600ns，形成2×2受控比较。Ricker不是原设备发射波形；这是激励一致性检验，最终都按实际源的复谱相除后比较同一精确20–170MHz/0.3MHz/501点，不把原始脉冲直接叫SFCW。

另外两站选原包已有c0058/c0108，各相隔10m，剖面里程分别169.25/159.25m，local中点69.25/59.25m；每站新算完整H1、去底连通砂岩H0、H0中仅x45–60m上覆砂岩→泥岩。原c0008高航高H1/H0复用，不重算。这个ROI对原站在远处，对c0108已接近收发点下方，因此不能把三站的替换都统一称“远处去杂波”。各站底砂评价窗由实际体素层序与完整材料谱的一次反射近似在求解前生成，峰±12ns，不按结果挑选。

坐标核对：原包裁剪模型x80–190m，local中点+80为裁剪前模型坐标；原几何`profile_x_offset_m=20`，因此剖面里程为local+100。原c0008的模型中点159.25m与包内剖面里程179.25m不是同一个坐标。新两站同时核对原native收发坐标和package_manifest的chainage，防止只凭文件名解释区段。既有报告中的剖面X指后者。

后处理不加尾渐消、AGC、移时或复幅相拟合；Hann和Blackman并列，差分先保留复数再取包络。还将直接以长时窗原生冲激场和实际Ricker源作离散因果卷积，检验独立Ricker正演，避免只证明两套频域代码一致。

## 执行契约

独立目录`E:/line9_basal_pair_20261008_r1/source_station_r1`。部署ZIP两端SHA256为`afd670175d82e757abfebb20b9b44d20d4bf40f67c784aa5755ba39232e97a9c`；最终执行合同SHA256为`5dab26e39413d8621925d35e42c40dc41507d88b92716f78c491d9d4a2bc9552`。完整源/时步/原生样本数、代码/运行时/输入哈希与容量门禁在合同中。最大9次、不得原地重试，单项≤30min/批次≤60min、600s会话心跳、USER_STOP和共享GPU锁沿用。无修改ROG原仓库、原始资料或旧冻结模型。

## 完成结果：9/9，另复用历史配对

9次启动、9次完成，无失败或重试；原生V4.0.0/float64/实际源、位置、时步及样本数通过审核。800ns为13569样本，1600ns为27136样本。接收序列的短窗在长窗前缀逐位一致（冲激、Ricker两者均为零差异），排除了更改记录长度反过来改变早期FDTD。

|模型|求解耗时 / s|
|---|---:|
|ricker_800|77.719|
|impulse_1600|150.844|
|ricker_1600|150.375|
|c0058_H1|78.406|
|c0058_H0|79.266|
|c0058_far|81.110|
|c0108_H1|80.796|
|c0108_H0|81.297|
|c0108_far|81.313|

累计单项耗时861.126s；含启动和审核的实际批次约14.5min，峰值所属进程RSS≤2.035GiB。完成审核SHA256 `37024d532cd9e73dae869215f463edc876ec36aec73cce842831a1c339ef2601`；回传ZIP SHA256 `c1f8d7d69210e8e5aa7c92ba41dcc132b87a52ed2a647e7ae57762c7e5eb40c0`。

### 源/记录长度并不是335ns强峰的主因，但晚窗不可照用

|源与记录长度|目标窗复差 / 冲激800ns（Hann）|同量（Blackman）|
|---|---:|---:|
|ricker_800|0.010732%|0.004435%|
|impulse_1600|0.010644%|0.004399%|
|ricker_1600|0.010732%|0.004435%|

三种新结果的目标窗主峰均335.994677ns；Hann相关均>0.99999999。独立平滑源正演也保留该峰，不支持它是冲激单独引发的非线性/高频污染。长冲激原生场与实际Ricker源作离散因果卷积（除40A，不拟合时延/系数），全序列相对误差约1e-14，250–450ns弱场窗≤2.79e-12；另用逐点直接求和而非FFT独立复算≤3.79e-14。

必须保留反例：600–800ns晚窗的新平滑源/延窗结果范数仅原冲激800ns的约8.43e-5（Hann），复差约100%；说明冲激有限记录截断确实能显著污染很弱的晚窗。它不足以解释335ns固定目标窗，但不能据此声称全部深晚时窗都可靠。源×时窗交互亦有记录，不把这些变化归为独立可加噪声。

### 夹层遮挡依位置变化，不能泛化为“整线都看不到底砂”

|剖面里程 / m|H1/模板相关 H/B|底砂差场/模板相关 H/B|H0/底砂差场范数 H/B|差场峰 / ns（Hann）|
|---|---:|---:|---:|---:|
|179.25|0.2749/0.4526|0.9616/0.9757|3.7007/3.6276|342.648|
|169.25|0.8996/0.9403|0.9061/0.9469|0.0566/0.0415|324.351|
|159.25|0.5823/0.6492|0.6370/0.7293|0.2334/0.2245|295.243|

H/B为Hann/Blackman。179.25m站位复现原结论：H0强于底砂差场约3.7倍且复内积接近反相，原场的局部模板相关很差。169.25m站位H0只有底砂差场约4–6%，H1与差场接近，底砂已在原场中占主导；这否定“高航高各站一律被夹层淹没”。159.25m站位差场占原场约97–99%，但差场本身与局部一次反射模板相关仅0.637/0.729，峰295.243ns比局部预测306.886ns早约11.64ns。此处不能把失配全部归咎于剩余背景；尚需核查非局部底砂路径/几何和模板的适用性。

两新站的ROI45–60m替换并非稳定“去杂波”：169.25m站位H0窗内范数反增13.82/18.00倍，159.25m约1.059/1.003倍、复差25.5%/20.2%。ROI相对新站已变近，删除引入新边缘并改变相消，不能视作通用改善方案。原站ROI敏感仍成立，但整线解释必须分区域，不按结果人为删除地层。

### 独立审核

501点DTFT与直接求和最大相对误差3.12e-12；逆变换独立复算最大1.64e-14。原生输出/输入/材料/几何哈希、source metadata和解析Ricker样本、短长前缀、事件9/9、两频窗复范数/相位及每张图哈希全部通过。只认证这些数值和模型对照，不认证全线或实测。

数值及可读审核：`artifacts/research_checks/2026-10-08_line9_source_station_results_r1/{analysis,delivery_audit,execution_contract,completed_verification}.json`；原生H5和完整复数SFCW缓存留私有`artifacts/local_checks/2026-10-08_line9_source_station_results_r1/`。


## 图件与可复现入口

本批新配对只覆盖三站中的两站，原站复用；图件已明确标记其他站位未新增配对，不把10m空隙画成新算的稠密数据。已有142道总场仅作站位定位和对照，不能冒称本批新B-scan。

```powershell
python scripts/analyze_line9_source_station_controls.py --source <retrieved-private-batch> --prior <previous-pair-results> --parent <basal-pair-prepared> --out <fresh-public-results> --numerical <fresh-private.h5> --cache <original-pkg2-cache.npz> --review <original-pkg2-audit.json>
```

- [source_window_grayscale.png](../../artifacts/research_checks/2026-10-08_line9_source_station_results_r1/source_window_grayscale.png)
- [source_window_envelopes.png](../../artifacts/research_checks/2026-10-08_line9_source_station_results_r1/source_window_envelopes.png)
- [native_convolution_check.png](../../artifacts/research_checks/2026-10-08_line9_source_station_results_r1/native_convolution_check.png)
- [three_station_grayscale.png](../../artifacts/research_checks/2026-10-08_line9_source_station_results_r1/three_station_grayscale.png)
- [three_station_envelopes.png](../../artifacts/research_checks/2026-10-08_line9_source_station_results_r1/three_station_envelopes.png)
- [existing_bscan_three_stations.png](../../artifacts/research_checks/2026-10-08_line9_source_station_results_r1/existing_bscan_three_stations.png)

## 证据边界

源型一致不等于现场天线/端口S21标定，延长到1600ns也不等于所有时窗/网格都收敛。三站控制是对原站解释的外延检验，不能认证整个195m测线、有限三维或实测可见性。ROI替换会引入边缘与改变交互；H1−H0包含材料改变的全部效应和下部PML材料延续，不是纯一次反射或实测clean，不用作训练标签。

这项检验也回应gprMax开发者关于线源冲激尾部不严格为零、有限记录截断需量化、可用独立平滑脉冲核查的建议。[作者相关正文说明，2023年10月31日及11月1日](https://groups.google.com/g/gprmax/c/hZQrvyrMzh8)。本次只采用“控制源型和时窗并实际量化”的方法依据；是否满足本模型弱回波精度由新数据决定。
