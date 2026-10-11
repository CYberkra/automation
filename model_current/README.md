# model_current —— 现行 gprMax 模型（固定文件夹）

**这个文件夹永远放当前最新推进的 gprMax 模型，文件夹名和位置永远不变。**
找最新模型：打开本文件夹即可，不用去别处找。

## 内容（即当前模型实体）

- `cases/x190_H0_local/profile.in`：无滑坡（H0）标准算例的 gprMax 文本输入。
- `cases/x190_H1_local/profile.in`：含滑坡（H1）标准算例的 gprMax 文本输入。
- `geometries/H0.h5` / `geometries/H1.h5`：两个场景的几何体素文件，profile.in 以相对路径 `../../geometries/` 直接引用，可原样运行。
- `geometries/line9_materials_v0_2.json`：材料契约（覆盖层 Debye、泥岩 εr18/σ0.006、砂岩 Debye）。
- `source.json`：源契约（官方 impulse，起始时间 0）。
- `acquisition.json`：采集设计契约（站位/道距等）。

## 规则

1. 模型有任何经批准的修改（材料/源/几何/采集），**先**把本文件夹旧内容整体移到 `model_archive/<日期>_<说明>/`，**再**放入新内容。旧版本永不覆盖、永不删除。
2. 本文件夹只放"模型定义"（输入），不放求解结果。每批求解的结果、图件和验收见 `configs/research/current_model.json` 指向的报告与证据目录。
3. 修改同步义务：改动本文件夹的同一工作单元内，更新 `configs/research/current_model.json` 并在 `docs/research/decision_log.md` 登记。

当前版本对应最新冻结批次：`line9_fine_wavefield_20261011_r1`（2026-10-11）。
