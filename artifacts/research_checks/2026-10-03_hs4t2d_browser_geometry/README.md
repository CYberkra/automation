# 已运行的 HS4T2D 二维剖面

`hs4t2d_arch.obj` 对应 121 道起伏模型；`hs4t2d_step.obj` 对应 121 道阶梯对照。两者共用 `hs4_scene.mtl`，蓝灰为基岩、棕色为覆盖层。Online 3D Viewer 可同时选择一个 OBJ 和 MTL 导入；设置 Z 轴朝上。

本目录模型是已运行 `.in` 箱体的 XZ 平面剖面，截面 Y=0.025 m；输入的 Y 单胞厚度为 0.05 m，二维物理模型沿 Y 不变，不代表有限厚度的三维目标。为聚焦地层，仅绘制 Z=0–12 m，原空气域 Z=12–33 m 未绘制；源/接收 Z=27 m 没有投影到地层内。原比例、不平滑、不纵向夸大。

起伏模型 48 箱高度与原 `hs4t2d_transect_table.npz` 一致，Z=8.85–9.25 m，覆盖深度 2.75–3.15 m。阶梯在 X=6 m 左右分为 Z=9.4/8.6 m，覆盖深度 2.6/3.4 m。两类模型的全部 121 道输入箱体分别一致，485 个输入身份核验通过。H5 仅读 Position 元数据（实际 Y=0），无接收数组或求解器调用。

`verification.json` 记录独立重读 OBJ 面顶点对照归档箱体的结果，表面积合计各为 144 m²；`export_provenance.json` 给出源码和输入 SHA-256。`manifest.json` 覆盖六份导出核心文件。该导出不替代求解器材料网格或 B-scan 物理验收。

复跑须使用新目录（NumPy/h5py 环境）：

```powershell
python scripts/export_hs4t2d_browser_geometry.py --out artifacts/local_checks/hs4t2d_browser_reexport
```

独立网站：[Online 3D Viewer](https://3dviewer.net/)。原作者[用户手册](https://3dviewer.net/info/index.html)支持本地 OBJ/MTL 及 GitHub 文件加载；[嵌入查看源码](https://github.com/kovacsv/Online3DViewer/blob/master/source/website/embed.js)支持相机/正交参数，用于正面显示剖面。完整站点每次导入默认 Y 向上，需手动改 Z；嵌入链接直接指定 Z 向上正面相机。
