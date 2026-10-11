# 六项精细波场诊断

先打开[可点击报告](wavefield_report.html)，再看[覆盖层工作带往返证据](x170_cover_roundtrip.png)。[完整研究报告](../../../docs/research/2026-10-11_fine_wavefield_results.md)说明物理机制、原生快照的高频限制及复现。

- 6项gprMax4.0.1原生FP64，8122帧。源、主接收记录同既有基线逐字节相同；模型、材料、激励、正式SFCW不变。
- 当前证据支持覆盖层重复返回与侧偏传播；底砂相关响应约399→356ns已存在。不支持把170m232ns强带主要归因右侧PML。
- 原始冲激动画含工作频带外格点纹理；工作带结论用完整探针的准确501点official direct。动画不是工作带场，不把亮纹当唯一射线路径。
- 32PNG/10GIF；局部细图使用2.5cm原生数组。预览约缩图3倍，完整局部动画206–351MB/段及约90.6GiB原快照留ROG。全域原生动画体量允许，随本目录入库。
- R2只修正190m诊断门窗误取首次反射旁瓣及显示轴、增加晚期关键帧共同尺度细图；未新增求解或改变正式数值结果。

验收：[delivery_verification.json](delivery_verification.json)。原生求解身份：[execution_contract.json](execution_contract.json)、[completed_verification.json](completed_verification.json)。当前文件哈希：[final_public_manifest.json](final_public_manifest.json)。原R1取回身份另存`public_retrieval_manifest.json`，修订图不冒充原图哈希。

私有原生探针、46个关键帧及数值数组：`D:\自动处理\artifacts\local_checks\2026-10-11_line9_fine_retrieval_r1\data`。

ROG原生全集：`E:\djh\line9_fine_wavefield_20261011_r1\execution\prepared\cases`；完整动画：`E:\djh\line9_fine_wavefield_20261011_r1\delivery_r1`。ROG本批owned进程已退出。

这不是3D/真实天线、现场复电性、网格/PML/记录长度全面收敛或实测一致性验收；未授权任何后继新求解。
