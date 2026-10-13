# 当前模型：V3.1连续测线严格同位置对照

用户已批准150–180m、51道×2网格，共102解。工作网格精确1/30m，2.5cm为对照。20–110MHz/301、50×50m、平滑y、1100ns、3mPML、GPU single、官方impulse/direct/Gaussian/补零4。

`profile.in`代表150m粗网格；`profiles/`含两档全部102份输入；`geometries/`含对应体素与常数材料JSON；`acquisition.json`冻结实际坐标和全部参数；`source.json`说明源及SFCW。

102解已在ROG的E:/djh/line9_matched_bscan_20261013_r1完成，不重复启动。此目录是模型实体与哈希入口，不是完整运行胶囊。完整私有包的data目录另含input_manifest、执行回执、原始输出与已执行study.py；包内code/analyze_line9_continuous_bscan.py可只重放分析。run_study.py是冻结worker副本，会启动正演，不是验收入口。

源高度采用两网格共同0.1m节点，最多改原高度5cm；域原点为站位−25.05m，实际收发里程为站位±0.65m。源格点0.1m不是计算网格。

旧模型完整保存在model_archive/2026-10-13_before_matched_bscan/。新几何/源为已批准版本；旧结果保持历史，不改标。

新V3.1体素H5只在本机及完整私有数据包中交付，不上传Git。纯克隆需取回本轮完整包的data/geometries目录，并逐项核对manifest哈希后才能运行。

结果：执行和变换核验通过，原网格筛查0/51通过；不能标记已收敛。见docs/research/2026-10-13_line9_matched_bscan_results.md。
