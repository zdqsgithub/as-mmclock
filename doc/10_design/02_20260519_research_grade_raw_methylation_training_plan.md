# Research-Grade Training Plan: Raw Mouse Methylation Data to Interpretable Age Clock

Date: 2026-05-19

## 1. Project Goal

本计划的目标是从小鼠甲基化原始实验数据出发，建立一个科研级、可复现、可解释的年龄预测模型。

当前项目边界沿用 v13/v14 结论：

- 当前有效 headline 只限于 `support-covered` chronological-age prediction。
- `GSE121141 old104+ brain_cortex/heart/lung` 是 stress-test/blocker metric，不是当前 headline pass/fail。
- 新的 full-lifespan old target-tissue claim 需要 Route A 生成/合作数据，或另行审批的 Route B minimal FASTQ/Bismark pilot。
- 任何 biological-age 或 CR/rapamycin 结论都必须来自真实 held-out predictions，不能来自 metadata、占位 AUC 或训练集残差。

最终模型要同时满足两个目标：

1. 年龄预测：从 RRBS/WGBS/compatible bisulfite methylation 数据预测小鼠 chronological age。
2. 生物解读：把被模型使用的 CpG/region/latent factor 转化为可审计的生物学解释，包括年龄相关方向、组织/批次混杂、区域稳定性、邻近基因与染色质背景。

## 2. 原始实验信号到甲基化定量值

### 2.1 分子层：甲基化状态如何变成可测序信号

DNA methylation clock 的底层输入不是“年龄”，而是 DNA 上胞嘧啶修饰状态的群体读数。

标准 bisulfite sequencing 的核心化学逻辑是：

- 未甲基化 cytosine 在 bisulfite 处理后转化为 uracil，PCR/测序中表现为 T。
- 5-methylcytosine 和 5-hydroxymethylcytosine 通常抗转化，在测序中仍表现为 C。
- 因此常规 bisulfite/RRBS/WGBS 的读数更准确地说是 `5modC = 5mC + 5hmC`，除非实验设计使用 oxBS/TAB/其他 5hmC 分辨方法。

这意味着模型的 beta value 不是单个细胞的二值状态，而是组织样本中多个细胞、多个等位基因、多个 DNA 分子在某 CpG 或 region 上的 modified-cytosine fraction。

### 2.2 文库层：RRBS/WGBS/EM-seq 的差异

本项目主线继续优先支持 RRBS/WGBS 和兼容的 bisulfite methylation processed schema。

- RRBS：MspI 等 restriction enzyme 富集 CpG-rich fragments，成本低，适合多样本时钟，但跨数据集 CpG 覆盖差异很大。
- WGBS：覆盖更全面，但测序成本和数据量更大；和 RRBS 混合时不能假设所有 RRBS CpG 都能被 WGBS 足够覆盖。
- EM-seq：酶法替代 bisulfite，减少 DNA 损伤，但进入项目 benchmark 前必须单独验证 schema、assembly、coverage 和 region overlap。

项目规则：

- RRBS 不默认去重。RRBS 的 restriction fragments 会天然产生相同起止位点，盲目 dedup 可能删除真实分子信号。
- WGBS 是否 dedup 必须由实验协议、UMI、library complexity 和 duplicate metrics 决定，并写入 manifest；不能混入未记录的 preprocessing 差异。

### 2.3 仪器层：Illumina 信号如何成为 FASTQ

Illumina SBS 的核心是 sequencing-by-synthesis：

- flow cell 上的 DNA cluster 进行逐轮碱基合成；
- 每轮加入带荧光标记的 reversible terminator dNTP；
- 仪器成像记录每个 cluster 的荧光颜色和强度；
- base caller 把图像强度、相位/预相位和背景噪声转为 A/C/G/T 与 Phred quality。

对甲基化项目而言，FASTQ 中的 `C` 和 `T` 不是普通 SNP 信号，而是经过 bisulfite/EM-seq 转换后的 epigenetic-state proxy。因此下游比对必须使用 bisulfite-aware aligner，而不是普通 DNA aligner 直接解释 C/T mismatch。

### 2.4 计算层：FASTQ 到 beta value

标准定量路径固定如下：

1. `FASTQ` QC：read count、read length、Q30、adapter content、per-base quality、GC/AT skew。
2. trimming：adapter、低质量尾部、RRBS MspI 相关端部偏差。
3. Bismark/Bowtie2 alignment：对 bisulfite-converted reference 进行方向性或非方向性比对。
4. methylation extraction：对 CpG cytosine 统计 methylated 和 unmethylated calls。
5. coverage file：形成 `chrom, start, end, methylation_percent, count_methylated, count_unmethylated`。
6. beta value：`beta = count_methylated / (count_methylated + count_unmethylated)`。
7. region aggregation：按 5kb 或预注册 region 聚合 CpG beta，默认使用 `mean(beta)` 或 coverage-aware sensitivity analysis。

同时保留以下不直接作为 headline age predictor 的技术 QC 特征：

- total reads、mapping rate、dedup/library complexity、conversion efficiency；
- M-bias、coverage distribution、CpG count per region；
- mean beta、global methylation proxy、missingness；
- reference assembly、liftover status、sex/MT exclusion count。

这些技术特征主要用于判断批次/文库/仪器混杂，不能在未审计情况下作为生物年龄解释。

## 3. 生物物理与生物学解释框架

### 3.1 从分子状态到组织 beta

beta value 的物理含义是采样到的 DNA molecule fraction：

```text
beta(CpG_i) = methylated_reads_i / (methylated_reads_i + unmethylated_reads_i)
```

它受三类因素共同影响：

- 分子修饰：DNMT/TET 活性、maintenance methylation、被动/主动去甲基化。
- 细胞组成：血液免疫细胞比例、组织细胞类型比例、衰老相关细胞群变化。
- 技术采样：coverage、PCR amplification、bisulfite conversion、restriction fragment capture。

因此，模型解释必须避免把所有重要 region 都直接称为“衰老机制”。更严格的说法是：这些 region 是与年龄相关且通过技术/组织混杂审计的 methylation signals。

### 3.2 年龄相关信号的生物学假设

可解释分析优先验证下列假设：

- 跨组织共性年龄信号：多组织中方向一致的 CpG/region，可能反映 conserved aging program。
- 组织特异年龄信号：只在 cortex/heart/lung/liver/blood 等组织中稳定，可能反映 tissue-specific maintenance 或细胞组成变化。
- region-level 稳定性：相邻 CpG 共同变化比单 CpG 更抗 RRBS coverage dropout。
- biological-age signal：CR/rapamycin 等干预组是否在真实 held-out residual 中表现出 age deceleration。

### 3.3 可解释特征层级

解释必须分层输出，不只看单个模型的重要性：

1. 单 CpG：beta、coverage、age correlation、presence rate。
2. 5kb region：mean beta、n CpGs、presence、age direction、fold stability。
3. region cluster/embedding：SVD/NMF latent factor、cluster_id、top loading regions。
4. annotation：CpG island/shore/shelf、promoter/enhancer/gene body、重复序列、染色质状态、邻近 gene。
5. confounding audit：dataset、tissue、sex、strain、library type、coverage、mean beta。

如果某 region 的 age association 同时被 tissue 或 dataset 强烈解释，则它应标记为 `confounded_or_context_specific`，不能作为泛化年龄机制证据。

## 4. 数据与实验设计

### 4.1 Route A 推荐 wet-lab / 合作数据设计

若目标是解决当前 full-lifespan old target tissue 缺口，Route A 是优先策略。

最低科研级设计：

- 物种/品系：mouse，优先固定 C57BL/6 或明确记录 strain。
- 组织：`brain_cortex/cortex`、`heart`、`lung`，每个组织独立建模/评估。
- 年龄：young/mid/old 覆盖，old 必须包含 `>=104w`；年龄必须 sample-specific exact age。
- 每组织年龄层样本：建议每个 age bin 至少 `n>=6`，更稳妥为 `n>=10`，以满足 support-covered benchmark。
- sex：平衡或预注册单性别；不可混杂在 batch 中。
- batch：sample extraction/library/sequencing 随机化；不能让 age=tissue=batch 完全重合。
- assay：RRBS 或 WGBS；必须记录 protocol、read layout、read length、adapter、conversion kit、reference assembly。
- control：包含 unmethylated spike-in 或其他 conversion efficiency proxy。
- 交付：sample sheet、file manifest、processed methylation、checksum、raw FASTQ 可追溯路径。

Route A 数据进入训练前必须通过：

- metadata gate；
- adapter smoke；
- matrix gate；
- RALPH Learn readiness；
- 显式训练授权。

### 4.2 Route B public/raw pilot 规则

Route B 只用于公开数据 rescue 或最小 raw pilot。

- 先查官方 GEO/SRA/ENA metadata；
- 只允许 2-3 sample pilot；
- 只在 sample-specific age、bulk target tissue、assay、run size、FASTQ URLs、checksum 都通过后下载；
- Bismark pilot 只用于判断是否能构建可用 region matrix；
- common 5kb regions `<50000` 时停止，不能训练。

v14 的 GSE83947 raw pilot 已经证明：工具链可执行，但 common region overlap 只有 `1814`，所以该数据集保持 auxiliary-only。

## 5. 训练数据矩阵设计

### 5.1 首选 feature type

主线 feature 仍为 5kb region beta matrix：

- 行：region_id，例如 `chr:start-end`；
- 列：sample_id；
- 值：region mean beta；
- sidecar：region stats，包括 chrom/start/end/n_cpgs/n_samples_present/mean_beta/std_beta/coverage summary。

选择 region 而非单 CpG 的理由：

- RRBS 单 CpG 跨数据集覆盖不稳定；
- region mean 可以降低 read-depth dropout；
- region 更适合 biological annotation 和 cluster-level interpretation；
- 项目历史 v3-v14 证明 common-region gate 是跨数据集可用性的核心门槛。

### 5.2 训练输入标准列

所有 predictions 必须输出标准 schema：

- `sample_id`
- `dataset_batch`
- `tissue`
- `intervention`
- `age_days_true`
- `age_days_pred`
- `age_weeks_true`
- `age_weeks_pred`
- `residual_weeks`

所有 metadata 必须至少包含：

- `sample_id`
- `age_days`
- `age_weeks`
- `tissue`
- `sex`
- `strain`
- `intervention`
- `dataset_batch`
- `assay`
- `metadata_source`

## 6. 无泄漏训练方案

### 6.1 固定评估层级

按项目规范，评估顺序如下：

1. Matrix gate：common 5kb regions、metadata overlap、age coverage、beta range、sex/MT exclusion。
2. GroupKFold：按 `dataset_batch` 分组，评估跨数据集泛化。
3. LODO：每个 dataset 做 leave-one-dataset-out。
4. Target held-out：`all_except:GSE121141 -> GSE121141`、`all_except:GSE80672 -> GSE80672`。
5. Support-covered headline：只统计同 tissue train_n `>=10` 且 age gap `<=8w` 的 held-out rows。
6. Unsupported stress-test：单独报告 GSE121141 old104+ 等超出支持范围的样本。

### 6.2 fold-internal preprocessing

以下步骤必须只在 train fold/train dataset 拟合：

- feature presence filter；
- imputation；
- scaling / robust scaling / quantile transform；
- feature selection；
- SVD/NMF embedding；
- target transform；
- model fit。

禁止：

- 全数据先筛 feature；
- 全数据先拟合 imputer/scaler/embedding；
- 在 unsupported old104+ stress metric 上做 autoresearch；
- dummy CR-AUC；
- human clock CpG 直接映射到 mouse。

### 6.3 模型路线

科研级默认路线不从 deep learning 开始：

1. Baseline：ridge / elasticnet。
2. Nonlinear baseline：LightGBM。
3. Sensitivity：Random forest 只作为辅助，不作为主搜索。
4. Embedding：SVD/NMF only，且如果参与预测必须 CV-safe；如果用于解释必须明确是否 post-hoc。
5. Deep learning：除非新 RFC 证明样本量、组织覆盖和 matrix gate 足够，否则不启用 MLP/CNN/Transformer。

## 7. 解释分析计划

每个通过 benchmark gate 的最佳模型必须输出：

- top regions by fold-stable importance；
- ridge/elasticnet coefficient direction；
- LGBM permutation/SHAP fallback importance；
- age correlation and monotonicity；
- tissue/dataset/coverage confounding score；
- latent factor / cluster summary；
- annotation enrichment；
- CR/rapamycin residual stratification when true held-out intervention data exist。

解释报告必须回答：

1. 哪些 region cluster 驱动年龄预测？
2. 它们在 fold/dataset/tissue 间是否稳定？
3. 这些 signal 是年龄、组织、批次、coverage 还是 cell composition 驱动？
4. 是否存在干预组 age deceleration？该结论是否有 shuffled sanity？
5. 是否能提出可实验验证的 candidate region/gene，而不是只给模型权重列表？

## 8. QC Gates and Stop Rules

### 8.1 Raw/Bismark gate

每个 raw pilot 或 production ETL 必须记录：

- FASTQ URLs、bytes、md5；
- reference assembly and SHA256；
- Bismark/Bowtie2/Samtools versions；
- mapping rate；
- methylation extraction report；
- coverage distribution；
- CpG/region count；
- beta min/max；
- sex/MT exclusion；
- common 5kb region overlap。

### 8.2 Matrix gate

进入 Learn/Benchmark 前必须满足：

- sample metadata overlap `>=95%`；
- exact age coverage `>=95%`；
- beta range within `[0,1]`；
- common 5kb regions with reference `>=50000`；
- dataset/tissue/intervention counts recorded；
- no duplicate sample_id；
- assembly/liftover traceable。

### 8.3 Benchmark gate

进入 constrained optimization 前必须满足：

- random-label sanity passes：`abs(r)<0.2` and MAE returns to random-like level；
- shuffled CR/intervention sanity passes when intervention metrics are reported；
- support-covered headline and unsupported stress metrics are reported separately；
- prediction schema is standard；
- result manifest records all leakage controls。

## 9. Success Criteria

### 9.1 当前模型可声明成功

如果使用现有 v13 boundary：

- support-covered MAE 维持或优于 `20.286w`；
- unsupported stress-test 单独报告；
- GSE121141 old104+ 不作为 headline 成败；
- CR AUC 只作为研究级 held-out validation。

### 9.2 新 full-lifespan target-tissue clock 成功

如果 Route A 或新 Route B 数据进入 headline：

- old `brain_cortex/heart/lung` 有 same-tissue age support；
- GSE121141 old104+ MAE 相对 v7.5 baseline `75.386w` 改善 `>=10w`，或经预注册新 benchmark 证明旧指标不再适合作为 headline；
- GSE121141 all-age held-out MAE `<=40.033w`；
- GroupKFold MAE `<=25.767w`；
- sanity checks pass；
- biological interpretation report identifies stable, non-confounded region clusters。

### 9.3 失败或转向条件

任一条件成立，应停止模型调参：

- no P1/P3 data after official refresh and pilot gates；
- common 5kb regions `<50000`；
- old target tissue remains unsupported；
- calibration only improves target set but fails non-target LODO；
- sanity check fails；
- model importance dominated by batch/tissue/coverage technical features。

## 10. Execution Plan

### Phase A: Raw Signal SOP Freeze

- 固化 Bismark/FASTQ/processed schema SOP。
- 保留 v14 GSE83947 作为 negative pilot example：技术可跑，但 matrix gate failed。
- 更新 Route A/B RFC，使合作方数据必须直接满足 matrix gate 所需字段。

### Phase B: Route A Data Intake

- 收集 sample sheet、file manifest、processed methylation。
- 跑 `validate_route_a_submission.py`。
- 跑 adapter smoke。
- 跑 Route A matrix builder。
- 只在 `ready_for_ralph_learn_pending_explicit_training_approval` 后准备 fixed benchmark command manifest。

### Phase C: Fixed Benchmark

- 不做 autoresearch。
- 先跑 ridge/elasticnet/lgbm 固定配置。
- 输出 GroupKFold、LODO、GSE121141 held-out、GSE80672 CR held-out、support-covered/stress split。
- 跑 random-label and shuffled intervention sanity。

### Phase D: Feature Interpretation

- 对通过 gate 的最佳模型生成 region/cluster解释。
- 标记年龄主导、组织主导、批次主导、coverage 主导。
- 给出可实验验证的候选 regions，不给过度 biological-age claim。

### Phase E: Promote or Stop

- 若达到 success criteria，进入 constrained autoresearch，最多 30-50 configs。
- 若未达到，写 blocker/RFC，转 Route A 数据生成或新 Route B candidate，不继续在旧矩阵上调参。

## 11. Required Sources

- Illumina SBS technology: https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html
- Bismark project and methylation extraction behavior: https://www.bioinformatics.babraham.ac.uk/projects/bismark/
- Bismark genome preparation: https://felixkrueger.github.io/Bismark/options/genome_preparation/
- Krueger and Andrews 2011, Bismark: https://academic.oup.com/bioinformatics/article/27/11/1571/216956
- Gu et al. 2011, RRBS protocol: https://www.nature.com/articles/nprot.2010.190
- Stubbs et al. 2017, mouse multi-tissue clock: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1203-5
- Meer et al. 2018, whole lifespan mouse multi-tissue clock: https://elifesciences.org/articles/40675
- Simpson et al. 2023, region-based RRBS clock: https://mayoclinic.elsevierpure.com/en/publications/region-based-epigenetic-clock-design-improves-rrbs-based-age-pred/
- 5hmC/bisulfite caveat: https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-017-0123-7

