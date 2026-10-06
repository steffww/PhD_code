# ONE-IBD 一站式分析数据（2026-10-06）

这个文件夹只放**自己重跑分析与绘图必需的输入**。所有数据均为**年龄≥18 岁，并剔除已明确记录的 2 型糖尿病、癌症史或长期抗生素使用阳性者**；原始项目文件未移动或覆盖。`manifest.csv` 记录每个副本的来源、大小和 SHA-256。

## 推荐先打开

| 文件 | 用途 |
|---|---|
| `metadata/baseline_core_450.csv` | 一人一行的精简临床协变量、心理评分、IBD 表型和资格状态；**450 人＝283 IBD＋167 非 IBD**。 |
| `metadata/baseline_full_450.csv` | 更完整的原始临床字段、疾病表型来源值与冲突标记。 |
| `metadata/followup_long_323.csv` | 湘雅 IBD 基线后 **323 次访视／170 人**；含访视日期、PHQ-9/GAD-7、IBDQ、部分活动度和原始住院/手术字段。一人可以有多行。 |
| `microbiome/mt_ibd_283.rds` | microeco `microtable`：**283 人、223 个物种**，适合 IBD 内分析。 |
| `microbiome/mt_all_449.rds` | microeco `microtable`：**449 人、268 个物种**，适合 IBD/非 IBD 比较；其中 1 名临床保留者经全人群物种过滤后全零而被移除。 |
| `metabolomics/full_panel/` | **1,017 个特征×299 人**（250 IBD＋49 非 IBD）的原始浓度、重新标准化矩阵、特征注释和 RDS。 |
| `metabolomics/HMDB_only/` | 同一 299 人的 **797 特征**独立处理分支；不要与 full panel 的标准化规则混用。 |
| `paired_ibd_250.csv` | 菌群与粪便代谢组均可用的 **250 名 IBD** 样本 ID。 |

`metadata/excluded_known_positive_44.csv` 保留被剔除者的 ID 与原因，供核对分母。**本包没有复制旧模型系数、p/q 或旧答辩图**，因为那些结果仍来自仅排除未成年人的旧分析。

## 代谢组各文件

- `raw_*.csv`：已按本版资格筛选的 assay 浓度矩阵，**特征为行、样本 ID 为列**；这不是质谱仪 RAW 谱图。
- `normalized_*.csv`：在这 **299 名保留者**上重新计算的分析矩阵，不是从旧成人 z-score 中简单删除列。
- `feature_annotation.csv`：检测特征的现有注释。
- `input_eligible.rds`：R list，包括 `raw`、`mat`、`meta`、`meta_enriched`、`lambda` 和 `sample_medians`。full panel 为 MedianNorm→广义 Log2Norm→AutoNorm；HMDB-only 为 NULL→广义 Log2Norm→AutoNorm。
- `metadata/fecal_model_metadata_299.csv`：两分支共用的原模型字段和样本次序。

## R 读取示例

```r
metadata <- read.csv("metadata/baseline_core_450.csv", check.names = FALSE)
followup <- read.csv("metadata/followup_long_323.csv", check.names = FALSE)
mt <- readRDS("microbiome/mt_ibd_283.rds")
otu <- as.matrix(mt$otu_table)       # 物种 × 样本
micro_meta <- mt$sample_table
mbx <- readRDS("metabolomics/full_panel/input_eligible.rds")
metabolites <- mbx$mat              # 1,017 特征 × 299 样本
stopifnot(all(colnames(otu) %in% metadata$sample_id))
stopifnot(all(colnames(metabolites) %in% metadata$sample_id))
```

若用 CSV，请用首列作为 feature ID，保持样本列名原样，并按 `sample_id` 而不是按行顺序连接临床表。

## 必须保留的解释边界

1. **IBD 的癌症史未被系统记录。** 这 283 人是 `cancer_history=unknown`，不是 `no`；本包只能称为“剔除已知阳性后的暂保留数据”，不能称为已证实全员无癌症史。
2. 其他零散合并症记录不是系统的是/否筛查，未擅自追加排除。重复疾病字段有冲突时，完整 metadata 保留冲突标记，协调值留空。
3. 随访仅纳入有可靠表头的湘雅 Sheet1，且访视日期晚于基线；Sheet2 的结果字段未静默合并。住院/手术字段是原始记录，尚不能直接作为“IBD 复发或进展”结局。
4. 两个 microeco 对象在各自保留人群中**独立重新过滤物种**；两条代谢组分支也各自重新标准化。不能直接套用旧成人版的候选菌、代谢物或显著性结果。

## 更新：IBDQ-32 四个领域分数（2026-10-06）

已在 `metadata/baseline_core_450.csv`、`baseline_full_450.csv` 及 `followup_long_323.csv` 中加入 `IBDQ_bowel`、`IBDQ_systemic`、`IBDQ_emotional`、`IBDQ_social`；两个 microeco RDS 的 `sample_table` 与两个代谢组 `input_eligible.rds` 的 `meta_enriched` 也同步加入**基线**领域分数。原 `IBDQ_total` 和 `IBDQ_sub` 保留不改。

| 领域 | 题号 | 可能分数范围 |
|---|---|---:|
| 肠道症状 | 1、5、9、13、17、20、22、24、26、29 | 10–70 |
| 全身症状 | 2、6、10、14、18 | 5–35 |
| 情绪功能 | 3、7、11、15、19、21、23、25、27、30、31、32 | 12–84 |
| 社会功能 | 4、8、12、16、28 | 5–35 |

题号映射依据公开的 [IBDQ-32 临床试验统计分析计划，附录 3](https://cdn.clinicaltrials.gov/large-docs/26/NCT02891226/SAP_001.pdf)。各项按 1–7 分计，**越高表示健康相关生活质量越好**。本包采用保守的**领域内所有题目有效才计算领域总分，不插补**；完整 32 题时四领域之和应等于按题目重算的总分。不同领域题数不同，不应直接比较未经标准化的领域总分。

在保留的 283 名 IBD 中，277 人的 32 题全部有效；肠道/全身/情绪/社会领域分别可计算 **281/283/281/281** 人。323 次基线后访视中，256 次 32 题全部有效；四领域分别可计算 **259/264/263/265** 次。完整记录的四领域和与按题目重算总分逐条一致。**2 名基线患者**的原锁定 `IBDQ_total` 与 32 题重算值相差 3 和 4 分，已在 `IBDQ_total_locked_discordant` 中标记；没有静默修改原总分。

`IBDQ_sub` 是旧工程中“排除与 PHQ-9/GAD-7 重叠项的总分”，**不是**这里的四领域之一。IBDQ 情绪领域本身含情绪题目；分析 PHQ-9/GAD-7 时，把它当作普通疾病活动度协变量可能带来题项重叠。评分脚本为 `score_ibdq_domains.py`，RDS 同步脚本为 `sync_ibdq_rds.R`。
