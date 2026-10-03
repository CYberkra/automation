# HS4 / COL6：浏览器查看输入模型

这些文件复原当前归档 HS4 的输入箱体几何。单位为米，Z 向上，原始比例，无平滑、无纵向夸大。源文件和分箱表没有改动；没有调用 gprMax。

## 推荐先看界面

打开 [Online 3D Viewer](https://3dviewer.net/)，点左上角 **Open from your device**，选择本目录的 `hs4_interface.obj`。加载后点工具栏 **Set Z axis as up vector**，可选 **Orthographic camera**，用鼠标旋转、缩放。

这个文件只显示基岩顶面：X/Y 各 12 m，Z 为 8.6–9.4 m，总起伏 0.8 m。48×48 个 0.25 m 箱体的水平台阶和相邻竖直跳变均保留。网页没有单位元数据，数值按米解释。

## 看地层和收发位置

在同一网站中，**同时选择** `hs4_scene.obj` 和 `hs4_scene.mtl`。每次重新导入后设 Z 轴朝上。左侧 Meshes 有四个对象：

- `rock`：基岩，从 Z=0 到起伏界面。
- `cover`：覆盖层，从起伏界面到地表 Z=12。
- `Tx_COL6_markers`：红色发射位置标记。
- `Rx_COL6_markers`：绿色接收位置标记。

点 `cover` 旁眼睛隐藏覆盖层，露出基岩顶面；点 `rock` 旁适配按钮放大地层。看收发位置时再适配整个模型。空气域隐含，不画遮挡模型的空气箱体。H5 的实际位置 X=1.60、Z=27，距地表 15 m；0.12 m 八面体只是位置标记，不代表天线形状/尺寸。

## 第二个独立网站

打开 [Kitware Glance](https://kitware.github.io/glance/app/)，点 **Open → Local → browse**，选择 `hs4_interface.vtp`，再点 **Load**。展开数据设置，Representation 用 Surface，Color by 用 `z_m`，颜色范围 8.6–9.4。鼠标旋转查看；初始视角可能接近侧视。

VTP 与 OBJ 的坐标和三角面逐项相同，`z_m` 是 Z 坐标，不是覆盖深度。该 VTP 是表面几何，不能用它查看体内材料切片或电磁波场。

## 与 B-scan 对照时

`col6_positions_and_profile.csv` 给出 25 道的实际收发位置、中点和中点处界面。中点沿 Y=3.25–9.25 m，覆盖深度为 2.75–3.05 m。整张三维曲面与这一条测线覆盖的几何不同；中点剖面也只是几何参照，不是双站雷达回波的预测真值。B-scan 的纵轴若为传播时间，不可直接当作模型 Z 轴。

导出来自归档 `.in` 与分箱表，并读 H5 位置属性，没有重新导出求解器实际离散材料网格。网页导入成功和几何核验通过不构成地下回波形态或其物理成因验收。

## 复现与证据

在仓库根目录、具有 NumPy/h5py 的 Python 环境中：

```powershell
python scripts/export_hs4_browser_geometry.py --out artifacts/local_checks/hs4_browser_reexport
python scripts/check_hs4_browser_geometry.py --directory artifacts/local_checks/hs4_browser_reexport
python scripts/check_hs4_browser_geometry.py --directory artifacts/research_checks/2026-10-03_hs4_browser_geometry
```

导出目录须不存在；复跑写新目录。`export_provenance.json` 记录 51 个输入身份及源码哈希；`manifest.json` 只覆盖导出器生成的六份核心文件。`verification.json` 是独立重读 OBJ/VTP/CSV 的验收记录；`package_manifest.json` 另覆盖交付目录文件（不含自身与 ZIP）。截图记录第三方网站的实际导入，截图为展示证据，不替代数值验收。

原作者资料：[Online3DViewer 格式清单](https://github.com/kovacsv/Online3DViewer/blob/master/README.md)、[Glance 文件加载](https://kitware.github.io/glance/doc/loading_files.html)。2026-10-03 查阅相应正文并实测上述本地导入流程。
