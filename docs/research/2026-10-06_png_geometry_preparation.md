# PNG 到 gprMax V4 几何：准备阶段

状态：已完成官方工具与私有图片分区准备；场地颜色含义、尺寸仍待确认，尚未生成该图片的可运行 gprMax 几何。没有求解、训练、修改原图或旧 H5。

## 官方依据与本机核验

读取日期 2026-10-06。阅读范围：官方 Utilities 的 PNG 转换节、转换模块源码、本机同名模块及相关格式测试源码；不是整个手册或求解器复现。

- [官方 Utilities](https://docs.gprmax.com/en/latest/inc_Utilities.html) 的 `convert_png2h5` 可交互点选颜色，或以 Python `convert_png` 指定颜色。
- [官方转换源码](https://github.com/gprMax/gprMax/blob/master/gprMax/toolboxes/Utilities/convert_png2h5.py) 按选色顺序生成编号，将图片行列转换为 XYZ 几何数组，未选像素为 -1（导入时保留底材，不能简单说它自动成为空气）。
- [官方材料数据库](https://docs.gprmax.com/en/latest/material_databases.html)：当前 PNG 转换同时生成 `material_keys` 与 JSON 材料表，颜色不能推断本构参数；旧 TXT 不是当前 V4 几何导入的直接配套数据库。
- 本机 `D:\gprmax_v4_gpu_env\Scripts\python.exe` 的版本实读为 4.0.0。文档 latest 当前显示 4.0.1；以本机实际格式试件为准，不从网页版本推断本机版本。
- 官方转换 API 小试件通过：XYZ 数组、两个不变方向切片、图片方向映射、-1 透明规则、材料键和 JSON 空本构模板。没有调用 `run` 或任何求解器。该环境没有 pytest，尝试指定单项格式测试未能启动，随后使用直接 API 断言验证；不能写成 pytest 通过。

格式试件、版本和转换模块 SHA-256 保存在忽略目录 `artifacts/local_checks/2026-10-06_line9_png_labels_r1/official_converter_fixture/summary.json`。

## 已准备的私有图片候选

用户指出原图下方为地表/空气；图片不是自动按顶部为空气解释。源图片没有米制刻度和图例。配套材料 TXT 给出空气、土层、泥岩、砂岩及其本构参数，但不编码颜色。

使用 [prepare_colour_geometry_image.py](../../scripts/prepare_colour_geometry_image.py) 显式指定主色。原图保持只读；蓝色最后一行/列仅作为待核验边框候选另存，过渡色保留原始未知标签、最近颜色建议与平局掩码，不当作已确认材料边界。

私有产物均在 `artifacts/local_checks/2026-10-06_line9_png_labels_r1/`，不入 Git：

- `image_labels_not_solver_geometry.h5`：像素分区草稿，缺米制几何属性，**不是 gprMax 输入**；原方向保留。
- `palette_clean_candidate.png`：纯色分区建议，不是确认后的地质模型。
- `summary.json`：源及产物哈希、裁边候选、未知/更改像素与歧义记录。H5/PNG 回读一致，源 SHA-256 不变。
- `image_legacy_alignment.json`：按规范化像素位置抽取旧 H5，比较四种上下/左右方向及颜色编号置换的诊断；不证明旧 H5 是由这张 PNG 直接生成。

旧 H5 的位置诊断支持用户的上下方向提示，但不是四色到四种材料的一一映射：两个不同主色主要对应同一旧材料编号，不能只根据 TXT 行序将图片从上到下命名。因此暂不将图例或物理尺寸冻结为事实；没有修改用户目录中已有 cleaned 候选。

## 下一步

收到颜色与物理范围后，记录新几何契约，明确图片方向、XYZ 存储、裁边/过渡色策略、网格和薄方向单元数，再生成新的 H5 + V4 JSON + ParaView 预览。旧 H5 的 35×45 m 是其属性，不自动成为 PNG 的场地尺度。

如保留现场剖面形态，该模型属于九号线场地适配，不能进入独立九号线盲测或仿真开发集。按此前材料设计使用电性时另列为设计假设，不从图片获得现场电性。求解需要独立完整运行契约；本阶段只是几何准备。
