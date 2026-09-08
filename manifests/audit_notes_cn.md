# 投稿代码审计说明

## 已采取的保护措施

原始工作区代码未被修改或删除。投稿副本只做了可追踪的路径配置和论文主参数补丁；主分析脚本仍保留论文未展示的扩展分析代码块，避免破坏依赖和运行顺序。

## 已确认的论文分析

- 论文主 benchmark 使用四种 pipeline、all/promoter、strict/filled 和 hop0，共 16 种 Looplook 模式。
- ChIPseeker 在当前 master 中统一使用全部 nearest-gene 注释，而不是对 promoter 模式切换成 promoter-only ChIPseeker 集合。
- Signal-density background 定义为表达分析 universe 中既不属于 Looplook、也不属于 ChIPseeker 的基因。
- 论文 GSEA 的主分析设定为每次抽样 200 个基因、共 300 次迭代；density 和 cumulative moving average 面板均应来自 size200。
- 最终选图目录中使用了 `GSEA_Unique_Density_*_size200_newplot.pdf`，但没有发现 `GSEA_Unique_EffectSize_Forest_*_newplot.pdf` 或 `GSEA_Unique_ThreeSet_ES_Forest_*_newplot.pdf` 被纳入最终选图。
- 已补入并由作者确认 Looplook 靶基因生成入口 `generate_multi_tss_data_NO_bedpe.R`：脚本从处理后的 loops、factor peaks、TPM、ATAC 与组蛋白修饰轨迹开始，生成八类 assignment 对象和四个 TSS-window RData。
- FOSL2 的四个 TSS 窗口均有日志记录完成。历史运行日志记录的内部版本为 0.99.15；公开可复现实现统一发布为 Looplook 0.99.19。

## 投稿前必须确认的问题

1. `analysis/looplook_benchmark_v12_plot_data_master.R` 的投稿副本已固定为 hop0、size200 和 300 次迭代；旧的 100/500 主抽样默认值不再保留在活动配置中。
2. `analysis/figure_scripts/1.plot_unique_gsea_cached.R` 已限定读取 size200/hop0。论文实际使用的 density 部分读取 master 生成的 `tmp_unique`；未使用的 effect-size forest 分支不应作为投稿结果。
3. 同一脚本的 three-set 新计算分支含有 `!vapply(es3_raw, is.data.frame, ...)` 过滤方向问题。因为对应 forest 未用于论文，本快照没有擅自修复。若将来公开该输出，需单独修复、版本化并验证。
4. 上游 mapping 脚本不纳入投稿代码目录；相关软件、参考基因组、过滤阈值和 peak-calling 参数以论文 Methods 为准。
5. 本目录的 package snapshot 已确认来自 `/Users/zhangying/Documents/3.software/looplook/maniscript/maniscript/v5/脚本/looplook`，其 DESCRIPTION 版本为 0.99.19。
6. 当前入口脚本保存的是一次具体运行的 `cfg`，其中 `project_name` 仍带有 `brd4_arv825`，而 `target_file` 指向 FOSL2 peaks。公开运行前应按数据集修改配置，但不要改动核心分析逻辑。
7. `overwrite_consensus = FALSE` 表示复用已有 consensus BEDPE；若公开数据包不提供该文件，必须另外提供其生成步骤或将配置改为重新生成并完成验证。

## 不应纳入公开代码的内容

- 文件名含 `old`、`old1`、`old2`、`old3`、`副本` 或 `已恢复` 的历史版本。
- `.RData`、临时 `.rds`、iteration checkpoints、plot queue、日志和重复 PDF。
- 未进入论文、且没有在 Methods 或补充结果中报告的探索性 forest、ranking 和 sensitivity 输出。
- 包目录中的 `.git`、`.Rproj.user`、R CMD check 临时目录、生成的网站和本地历史文件。
- 上游 mapping 脚本；其分析参数由论文 Methods 描述，不随投稿分析代码发布。

## 建议发布策略

首次投稿可以提交此保守快照加数据与代码可用性说明。正式公开归档前，再建立一个经过四套数据完整重跑验证的 portable release；不要直接在本快照上进行大规模重构。
