# HS4 输入模型的第三方浏览器查看

用户希望独立查看当前模型，核对 B-scan 与展示几何的差异。本轮将归档 HS4 / COL6 输入箱体转换为第三方网站可读文件，不运行求解器，不读取实测资料，不改变处理链。

## 可直接使用的交付

目录：`artifacts/research_checks/2026-10-03_hs4_browser_geometry/`，操作说明见其中 `README.md`。

| 网站 | 本轮实际导入 | 用途 |
| --- | --- | --- |
| [Online 3D Viewer](https://3dviewer.net/) | `hs4_interface.obj`；`hs4_scene.obj` + `hs4_scene.mtl` | 旋转界面、隐藏覆盖层、查看实际收发位置 |
| [Kitware Glance](https://kitware.github.io/glance/app/) | `hs4_interface.vtp` | 用 `z_m` 为界面高度着色 |

原作者的 [格式清单](https://github.com/kovacsv/Online3DViewer/blob/master/README.md) 与 [Glance 加载说明](https://kitware.github.io/glance/doc/loading_files.html) 已读取相关正文。两个网站均用本轮模型实际打开，未使用 Share 发布功能。3D Viewer 每次导入后需设置 Z 轴朝上；完整地层需隐藏 `cover` 才能看到内部界面。Glance 导入 VTP 的初始视角接近侧视，旋转后观察。没有验证体数据切片、VTKHDF 或波场加载。

## 核验结果

`scripts/export_hs4_browser_geometry.py` 解析 25 份 COL6 `.in`，全部 2,304 个覆盖箱体相同，恢复高度与 `hs4_interface_binned_table.npz` 逐项一致。导出前后核验 51 个输入 SHA-256；H5 仅读 Position 属性，未读取接收时序。

`scripts/check_hs4_browser_geometry.py` 独立重读导出文本，验收 2,304 个水平箱面的坐标/面积/朝向、21.05 m² 竖直台阶面积、OBJ/VTP 坐标与面索引一致、覆盖层/基岩体积、50 个收发标记中心及 CSV 剖面。结果 PASS，详细数值见 `verification.json`。另写新目录重复导出，核心六份文件 SHA-256 一致。

- 界面 OBJ：4,107 顶点、7,888 三角面；X/Y 为 0–12 m，Z 为 8.6–9.4 m。无平滑、无纵向夸大，保留台阶。
- 地层与位置 OBJ：8,898 顶点、16,948 三角面；四个对象名称与浏览器实际显示一致。
- 基岩体积 1296.096875 m³，覆盖层体积 431.903125 m³，总和 1728 m³；空气域不绘制。
- H5 实际 COL6 X=1.60 m（输入标称 1.625 m），收发 Z=27 m，地表 Z=12 m；位置标记为直径 0.12 m 的符号八面体。
- COL6 中点沿 Y=3.25–9.25 m；中点覆盖深度 2.75–3.05 m，不能用整张三维界面的总起伏代替这一测线剖面。

## 结论边界

这是归档输入箱体几何复原，尚未独立核验求解器实际离散材料网格。网站可帮助用户排查几何展示的轴向、比例与测线选择；不会把传播时间轴变成深度轴，也不能单独说明回波形态为何不同。当前地下形态与边界/侧向散射的物理归因仍需同条件参考及数值对照，不因本次导出升级签认。

复现命令见交付 README；原始胶囊、冻结契约/G4、旧图及旧哈希均保持原样。新增脚本不进入默认数组回归链，无求解/训练/实测性能结论。本研究单元独立分支提交推送。
