# CUT&Tag Snakemake 7 + SLURM Pipeline

这是一个用于转录因子基因组结合位置发现的 CUT&Tag 双端测序数据分析 pipeline。项目基于 Snakemake 7 开发，适用于 SLURM 集群批量投递。pipeline 的核心分析流程如下：

```text
paired-end FASTQ
  -> fastp trimming / quality filtering
  -> Bowtie2 alignment
  -> BAM 过滤与分组 pooled
  -> deepTools BigWig / correlation / PCA / fingerprint / TSS profile
  -> MACS3 peak calling
  -> IDR（可选，恰好两个 treatment biological replicates 时运行）
  -> 并行 motif 分支：HOMER / MEME-ChIP（默认关闭内部 STREME）
  -> MultiQC / 可选结果打包
```
本 pipeline 默认合并同组 biological replicates（pooled）后进行主 peak calling，以增加可用于检测的测序深度；这不保证所有样本的灵敏度或特异性都会提高。用户可结合 IDR 和其他 QC 审查重复间一致性，不能用 pooled 结果掩盖单个重复的数据质量问题。

本文档是 pipeline 的用户文档，覆盖项目初始化、环境、配置、metadata、参考基因组、SLURM 投递、进度监控、结果解释和常见故障。


## 阅读导航

- [使用方法](#1-使用方法)
- [软件和环境](#2-软件和环境)
- [样本 metadata](#3-samplestsv样本-metadata)
- [比较关系](#4-comparisonstsv比较关系)
- [完整配置](#5-configyaml完整配置说明)
- [参考基因组准备](#6-参考基因组准备)
- [投递前检查](#7-投递前检查)
- [SLURM 投递](#8-slurm-配置和投递)
- [进度监控](#9-监控任务进度)
- [Pipeline 输出内容](#10-pipeline-输出内容)
- [结果打包](#11-结果打包)
- [故障处理](#12-常见故障和处理方法)
- [运行检查清单](#13-单次运行检查清单)
- [科学解释限制](#14-科学解释上的重要限制)
- [分发维护原则](#15-Pipeline-Release-维护原则)
- [软件和数据库引用](#16-软件和数据库引用)

---

## 1. 使用方法
### 1.1 Pipeline Release 文件结构
Release 目录存放 workflow 和脚本：
```text
cuttag_pipeline_release/
├── README.md
├── VERSION              # release 版本记录
├── run.sh               # 主启动脚本
├── submit.slurm         # SLURM投递脚本
├── init_project.sh      # 创建项目配置模版
├── config/
├── schemas/
├── scripts/
├── workflow/
├── profiles/slurm/
└── tests/               # 小型回归测试，不是分析输入数据
``` 
### 1.2 快速启动示例
```bash
# 1. 创建项目模板
/path/to/cuttag_pipeline_release/init_project.sh /path/to/my_cuttag_project

# 2. 编辑 config.yaml、samples.tsv、comparisons.tsv (见1.3)

# 3. 提交前检查
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/my_cuttag_project/config.yaml \
  --check-only

# 4. 查看 DAG，不投递任务
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/my_cuttag_project/config.yaml \
  --dry-run

# 5. 投递 SLURM driver
# 先激活 Snakemake 7；在此目录提交可使 driver 日志落在项目中
cd /path/to/my_cuttag_project
sbatch \
  --account=YOUR_ACCOUNT \
  --partition=YOUR_PARTITION \
  /path/to/cuttag_pipeline_release/submit.slurm \
  /path/to/my_cuttag_project/config.yaml
```

### 1.3 Project 文件结构及配置
project 文件夹是用户实际分析流程配置文件、结果和日志的存放位置，不要求与pipeline release存在于同一目录。
```text
my_cuttag_project/
├── config.yaml           # pipeline配置文件：pipeline运行参数
├── samples.tsv           # pipeline配置文件：样本 metadata
├── comparisons.tsv       # pipeline配置文件：实验对照组别
├── 01.RawData/           # 可选，按用户目录组织
└── results/              # pipeline 运行结果
```
Pipeline的`init_project.sh`将会复制配置文件到用户指定路径`/path/to/my_cuttag_project/`，用户需对`config.yaml`、`samples.tsv`、`comparisons.tsv`  共 3 个文件进行配置后才可启动分析。（样本、比较关系及分析参数见后续对应章节。）

---

## 2. 软件和环境

### 2.1 运行前提

本 Pipeline 需求的运行环境为：

- Linux 环境；
- Snakemake 7.x；
- Conda（自动环境管理必需）及可选 Mamba solver；
- SLURM 客户端（集群运行时）；
- 可访问项目数据、参考基因组和输出目录。

检查 Snakemake 版本：

```bash
snakemake --version
```

必须确认是 **Snakemake 7.x**。`run.sh` 不负责安装 Snakemake，只会检查当前 `PATH` 中是否存在 `snakemake`。
分析软件依赖见 `workflow/envs/`。使用 `null` 配置时，Snakemake 可创建分析环境，但前提是启动环境、Conda 和软件包源均已准备好；不会自行安装 Snakemake 或 Conda。

典型的环境激活方式：

```bash
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate snakemake7
snakemake --version
```

### 2.2 两种分析环境模式

#### 模式 A：使用 pipeline release 自带 environment YAML

项目配置`config.yaml`中保持：

```yaml
software:
  conda_envs:
    core: null
    homer: null
    meme: null
    idr: null
```

Snakemake 会根据以下文件创建隔离环境：

```text
workflow/envs/core.yaml
workflow/envs/homer.yaml
workflow/envs/meme.yaml
workflow/envs/idr.yaml
```

如果集群不能联网，建议管理员提前创建这些环境，或者使用模式 B。

#### 模式 B：使用服务器已有命名环境

例如：

```yaml
software:
  conda_envs:
    core: cuttag_pipeline_core
    homer: cuttag_homer
    meme: cuttag_meme
    idr: cuttag_idr
```

这些环境必须在投递前已经存在，pipeline 不会自动修复它们。环境名是用户自定义的，不要求与示例一致；`software.conda_envs.meme` 是 MEME Suite 环境名，即使 motif 方法名为 `memechip`，这里仍使用键 `meme`。

例如可以在允许安装软件的节点创建命名环境（只创建环境，不启动分析）：

```bash
PIPELINE=/absolute/path/to/cuttag_pipeline_release
conda env create -n cuttag_pipeline_core -f "$PIPELINE/workflow/envs/core.yaml"
conda env create -n cuttag_homer -f "$PIPELINE/workflow/envs/homer.yaml"
conda env create -n cuttag_meme -f "$PIPELINE/workflow/envs/meme.yaml"
conda env create -n cuttag_idr -f "$PIPELINE/workflow/envs/idr.yaml"
```

不应直接复用一个未经检查的共享环境。`idr_env` 只是过去某个服务器的环境名，并非普遍不可用；判断环境是否可用应检查实际依赖和功能，而不是依据名称。

### 2.3 Snakemake 启动环境与依赖

分析环境不包含 Snakemake 本身。启动环境还需要 `PyYAML`；`jsonschema` 用于配置 schema 校验，缺失时 preflight 会警告并跳过该项检查。Snakemake 的自动建环境还需要 Conda；Mamba 是可选的 solver，不是 Conda 的完全替代品。

若安装的 Snakemake 默认使用 Mamba，而服务器没有 Mamba，可在启动参数中明确选择 Conda：

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml \
  --dry-run -- --conda-frontend conda
```

`submit.slurm` 只接收项目配置，不会转发额外的 `run.sh` 参数。若要对 driver 使用这个设置，可在使用的 `profiles/slurm/config.yaml` 中增加 `conda-frontend: conda`。自动环境创建需要包源和网络；不要把“配置检查通过”理解为“环境已经创建并测试通过”。

---


## 3. samples.tsv：样本 metadata

示例模板：

```tsv
sample_id	group_id	replicate	label	r1	r2	batch	pair_id	description
TF_R1	TF	1	TF rep1	data/TF_R1_1.fq.gz	data/TF_R1_2.fq.gz	B1	P1	transient TF replicate 1
TF_R2	TF	2	TF rep2	data/TF_R2_1.fq.gz	data/TF_R2_2.fq.gz	B1	P2	transient TF replicate 2
IgG_R1	IgG	1	IgG rep1	data/IgG_R1_1.fq.gz	data/IgG_R1_2.fq.gz	B1	P1	control replicate 1
IgG_R2	IgG	2	IgG rep2	data/IgG_R2_1.fq.gz	data/IgG_R2_2.fq.gz	B1	P2	control replicate 2
```

字段说明（文件必须以 **Tab** 分隔，不能用空格代替，也不能把 Excel 的 CSV 直接改扩展名）：

| 字段 | 必需 | 说明 |
|---|---:|---|
| `sample_id` | 是 | 唯一的样本 ID，必须是路径安全字符 |
| `group_id` | 是 | 生物学分组，决定 pooled 分析 |
| `replicate` | 是 | 组内 biological replicate 编号，从 1 开始 |
| `label` | 否 | 图形展示标签 |
| `r1` | 是 | R1 FASTQ 路径 |
| `r2` | 是 | R2 FASTQ 路径 |
| `batch` | 否 | 批次记录，当前不参与批次校正 |
| `pair_id` | 否 | matched treatment/control 记录，当前 pooled 分支不使用 |
| `description` | 否 | 自由文本说明 |

### 3.1 一行代表什么？

一行代表一个 **biological replicate**，不是一个技术 lane。

例如同一个 biological replicate 有多个测序 lane，应在进入 pipeline 前先合并 FASTQ，或者明确整理为一个样本目录。不要把 lane1、lane2 填成两个 biological replicates。

### 3.2 `group_id` 和 `replicate` 的作用

- 同一 `group_id` 的样本会 pooled；
- pooled BAM 用于主 MACS3 peak calling；
- `replicate` 用于识别 biological replicate；
- treatment 组的 replicate 数决定 IDR 是否运行；
- `label` 只影响展示，不影响计算。

### 3.3 FASTQ 路径

相对路径相对于 `config.yaml`，例如：

```text
data/TF_R1_1.fq.gz
```

也可以使用绝对路径：

```text
/public/data/project/TF_R1_1.fq.gz
```

R1 和 R2 必须：

- 都存在；
- 不是同一个文件；
- 属于同一 paired-end 文库；
- read name 能够正确配对。

---

## 4. comparisons.tsv：比较关系

模板：

```tsv
comparison_id	treatment_group	control_group	description
TF_vs_IgG	TF	IgG	pooled TF versus pooled IgG
```

字段说明（文件必须以 **Tab** 分隔，不能用空格代替，也不能把 Excel 的 CSV 直接改扩展名）：

| 字段 | 必需 | 说明 |
|---|---:|---|
| `comparison_id` | 是 | 唯一的比较 ID |
| `treatment_group` | 是 | 必须存在于 `samples.tsv` 的 `group_id` |
| `control_group` | 是 | 必须存在且不能与 treatment 相同 |
| `description` | 否 | 比较说明 |

当前 control 是 pooled control。`pair_id` 不会自动改变比较关系。

---

## 5. config.yaml：完整配置说明

以下 `fastp`、`alignment`、`coverage`、`qc`、`macs3`、`idr`、`annotation`、`motif` 代码块均位于 `analysis:` 下，不能作为顶层配置粘贴。`project`、`software`、`reference` 和 `outputs` 是顶层键。优先编辑 `init_project.sh` 复制的完整模板，以避免缩进或缺少必需参数。


### 5.1 project

```yaml
project:
  name: cuttag_example
  output_dir: results
  samples: samples.tsv
  comparisons: comparisons.tsv
```

- `name`：项目名称，仅用于记录；
- `output_dir`：结果目录，可使用相对路径；
- `samples`：samples.tsv 路径；
- `comparisons`：comparisons.tsv 路径。

### 5.2 software

```yaml
software:
  conda_envs:
    core: null
    homer: null
    meme: null
    idr: null
```

填 `null` 使用 release 的 environment YAML；填写环境名则使用已有命名环境。

### 5.3 reference

```yaml
reference:
  assembly: B73_RefGen_v5
  fasta: /absolute/path/to/genome.fa
  gtf: /absolute/path/to/annotation.gtf
  bowtie2_index: null
  effective_genome_size: 2182075994
  chromosomes: ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
```

#### reference 配置原则

FASTA、GTF、Bowtie2 index、chromosome list 和 effective genome size 必须属于同一个 assembly。

#### bowtie2_index

填 `null` 时，pipeline 会根据本次 FASTA 在`results/00.reference/bowtie2/`建立 index。

如果使用已有 index，填写的是 prefix：

```yaml
bowtie2_index: /data/genome/v5/bowtie2/genome
```

不是：

```text
/data/genome/v5/bowtie2/genome.1.bt2
```

#### chromosomes

请确保染色体名称必须与 FASTA 完全一致。如果 FASTA 使用 `chr1`，就必须配置 `chr1`，不能自行改为 `1`。

#### effective_genome_size

该值同时用于：

- MACS3 `-g`；
- deepTools RPGC normalization。

它不是简单地复制 FASTA 总长度就一定合理。请根据物种、assembly 和分析口径确认来源，并在项目记录中保留依据。

检查 FASTA contig：

```bash
samtools faidx genome.fa
cut -f1,2 genome.fa.fai
```

### 5.4 fastp

```yaml
fastp:
  qualified_quality_phred: 15
  length_required: 36
  cut_front: false
  cut_tail: true
  cut_window_size: 4
  cut_mean_quality: 20
  average_qual: 0
```

用于 adapter/质量过滤。对于大的 FASTQ 可能需要较长 walltime，profile 中默认 fastp 为 4 小时。

### 5.5 alignment

```yaml
alignment:
  mode: end-to-end
  preset: very-sensitive
  no_mixed: true
  no_discordant: true
  insert_min: 10
  insert_max: 700
  mapq_min: 30
```

- 当前默认 paired-end end-to-end 比对；
- `mapq_min` 用于 BAM 过滤；
- `insert_min/max` 用于 Bowtie2 插入片段约束；
- 只有在明确理解后再修改 `no_mixed` 和 `no_discordant`。

### 5.6 coverage

```yaml
coverage:
  bin_size: 25
  smooth_length: 50
  extend_reads: 200
  normalize_using: RPGC
  ignore_duplicates: false
```

`ignore_duplicates` 只影响 BigWig 生成，不修改 BAM，也不等价于 MACS3 去重。

`bin_size` 是覆盖轨迹的采样宽度，`smooth_length` 是平滑窗口，两者不是同一个参数。BigWig 后续被相关性/PCA 和 TSS profile 使用，因此它的 bin size 和平滑会影响下游信号分辨率；不能靠把 TSS bin 调得更小恢复已在 BigWig 中丢失的细节。

当前命令传递 `--extendReads 200`。对于有正常 mate 的 paired-end fragment，deepTools 会优先按两个 mate 的实际位置确定 fragment 长度；用户提供的 200 bp 主要作为 singleton、mate 距离异常或跨染色体时的 fallback，而不是把所有 paired-end fragment 强行设为 200 bp。不要把它解释成 CUT&Tag 文库片段长度的测量值。

### 5.7 qc

```yaml
qc:
  summary_bin_size: 25
  correlation_method: spearman
  fingerprint_bin_size: 25
  fingerprint_samples: 100000
  tss_before: 3000
  tss_after: 3000
  tss_bin_size: 25
```

用于 deepTools 的相关性、fingerprint、PCA 和 TSS-centered signal profile；这些输出是 QC/可视化，不是独立的 TSS enrichment score。

#### deepTools 的 bin size 统一为 25 bp

当前模板将四项 bin size 均设置为 **25 bp**；配置键仍分别保留，方便用户按分析目的调整。虽然数值相同，它们的统计对象不同：

| 配置键 | 默认值 | 实际命令/输入 | 用途 |
|---|---:|---|---|
| `analysis.coverage.bin_size` | 25 | `bamCoverage --binSize`；sample/group BAM | 生成相对细粒度的 BigWig 覆盖轨迹 |
| `analysis.qc.summary_bin_size` | 25 | `multiBigwigSummary bins --binSize`；sample BigWig | 汇总 25 bp 窗口内的信号，用于相关性和 PCA；二者共用同一矩阵 |
| `analysis.qc.fingerprint_bin_size` | 25 | `plotFingerprint --binSize`；sample BAM | 抽样窗口的 read count 分布，用于 fingerprint |
| `analysis.qc.tss_bin_size` | 25 | `computeMatrix reference-point --binSize`；sample BigWig + GTF | 构建 TSS 上下游信号矩阵，用于热图和 profile |

`fingerprint_samples: 100000` 是抽样窗口数量，不是窗口大小。当前 fingerprint 直接读取 BAM，不经过 BigWig；相关性/PCA 和 TSS 则读取 BigWig。

四项分析的 bin size 已统一为 25 bp，仍通过各自的配置键传递。相关性/PCA 的全基因组 25 bp 汇总矩阵可能较大，应关注内存、磁盘和运行时间。fingerprint 的抽样窗口数量仍为 100000。跨样本、跨运行比较时，还应保持平滑、归一化和输入过滤条件一致；bin size 相同不代表这些分析的统计方法或输入相同。

TSS 默认窗口为上游 3000 bp、下游 3000 bp；建议各窗口长度能被 `tss_bin_size` 整除。当前 TSS 输出是聚合信号热图和 profile，**没有计算独立的 TSS enrichment score，也没有自动 QC 合格阈值**。

### 5.8 MACS3

```yaml
macs3:
  format: BAMPE
  qvalue: 0.05
  call_summits: true
  trackline: false
  broad: false
  keep_dup: auto
```

当前要求 paired-end `BAMPE`、narrowPeak、`call_summits: true` 和 `broad: false`。`keep_dup: auto` 为默认设置；schema 也允许正整数，改变它会影响 MACS3 保留重复片段的方式，不会物理修改输入 BAM。

pipeline 不做 Picard 物理去重，也不额外执行自定义 peak filter。

### 5.9 IDR

```yaml
idr:
  mode: auto
  threshold: 0.05
  rank: p.value
```

`auto` 的行为：

| treatment biological replicates | 行为 |
|---:|---|
| 1 | 跳过 IDR |
| 2 | 运行 IDR |
| >2 | 当前版本跳过，不静默截取前两个 |

IDR 是 replicate-level QC/高置信 peak 子集，不替代 pooled MACS3 peaks。treatment group 必须在 metadata 中存在且有样本；没有样本的 group 不是合法的“0 个重复”运行模式。

设置：

```yaml
mode: off
```

可关闭全部 IDR。

### 5.10 annotation

```yaml
annotation:
  enabled: false
  method: homer
  promoter_max_abs_tss_distance: 5000
```

当前 annotation 只支持 HOMER。
用户可自行使用`ChIPseeker`、`clusterProfiler`等开展下游分析。

### 5.11 motif

```yaml
motif:
  enabled: true
  methods: [homer, memechip]
  input_scope: all
  homer:
    size: 200
    lengths: [6, 8, 10, 12]
  memechip:
    flank: 100
    center_cut: 0
    min_width: 6
    max_width: 12
    meme_p: 16
    meme_nmotifs: 3
    streme_nmotifs: 0
    known_motif_db: null
```

#### 输入

motif 使用 MACS3 treatment peak summit 序列。

`flank: 100` 表示对 1 bp summit 区间两侧各扩展 100 bp，一般得到 **201 bp** 序列；靠近染色体边界时会裁剪。HOMER 则通过自身 `-size 200` 获取固定窗口，二者并非逐碱基完全相同的输入。

#### HOMER

HOMER 使用参考基因组匹配背景；没有配置 IgG `-bg`。HOMER 直接读取 summit BED 和参考 FASTA，不依赖 `motif_fasta` 或 MEME-ChIP，二者是独立分支。

#### MEME-ChIP

MEME-ChIP 使用：

```bash
-ccut 0
-meme-p 16
-meme-nmotifs 3
-streme-nmotifs 0
```

`streme_nmotifs: 0` 是显式关闭 MEME-ChIP 内部 STREME 的必要配置。

`center_cut: 0` 对应 `-ccut 0`，表示**不进行额外的中心裁剪**，使用完整的 summit FASTA 序列。它不是显著性 cutoff；`flank: 100` 仍控制 summit 两侧扩展，通常生成 201 bp 序列。与 `-ccut 100` 相比，MEME 搜索输入可能更长、耗时可能增加；本次不改变 motif 显著性或筛选阈值。

当前 rule 固定申请 16 个线程，`meme_p` 控制实际传给 MEME 的并行数；不要把 `meme_p` 配置为大于已申请线程数的值。`streme_nmotifs` 允许非负整数；大于 0 会重新启用内部 STREME，preflight 会警告，不再属于本项目推荐的默认运行方式。

`known_motif_db: null` 时不传递 `-db`。填写时应为有效的 MEME-format motif 数据库路径；当前 preflight 不检查数据库内容，应由用户自行核验格式、物种适用性及访问权限。没有数据库时，不能期待已知 motif 数据库匹配注释。

当前没有 external negative/control FASTA，不使用`-neg`参数。因此当前 motif 结果不是严格意义上的 treatment-vs-IgG differential motif enrichment。

### 5.12 outputs

```yaml
outputs:
  keep_trimmed_reads: true
  package:
    enabled: true
    archive: cuttag_results.tar.gz
```

`package.enabled: true` 时，所有分析目标完成后生成结果压缩包。

---

## 6. 参考基因组准备

建议先准备并检查：

```bash
samtools faidx /path/to/genome.fa
cut -f1,2 /path/to/genome.fa.fai
```

如果使用预建 Bowtie2 index：

```bash
bowtie2-inspect -n /path/to/index/prefix | head
```

确认：

- FASTA 与 GTF seqname 一致；
- FASTA 与 Bowtie2 index contig 一致；
- `chromosomes` 中的名称存在于 FASTA；
- effective genome size 口径明确；
- 参考版本没有混用。

如果 `bowtie2_index: null`，第一次运行会建立 index，可能需要较长时间和较大磁盘空间。

---

## 7. 投递前检查

### 7.1 preflight

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml \
  --check-only
```

推荐状态：

```text
Preflight summary: 0 error(s), 0 warning(s).
```

`error` 会阻止启动；`warning` 不一定阻止启动，但需要用户确认原因，例如缺少 `jsonschema`、无法验证预建 index，或 treatment 重复数大于 2。preflight 不检查 FASTQ 全部内容或 read name 配对，也不验证 motif 数据库内容。

preflight 能检查配置、metadata、路径和部分参考一致性，但不能完全保证：

- 节点上的动态库兼容；
- Conda 环境没有损坏；
- 资源足够；
- 作业一定不会超时。

### 7.2 dry-run

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml \
  --dry-run
```

检查：

- 是否出现预期 rule；
- 是否意外重算大量结果；
- 是否出现 missing input；
- motif 是否出现 `homer_motif` 和 `memechip`；
- 不应出现独立 `meme` 或 `streme` rule；
- 如果结果已存在，是否显示 `Nothing to be done`。

---

## 8. SLURM 配置和投递

### 8.1 配置 account 和 partition

编辑：

```text
profiles/slurm/cluster.yaml
```

修改：

```yaml
__default__:
  account: YOUR_ACCOUNT
  partition: YOUR_PARTITION
```

不要把具体实验室的 account 和 partition 固定写入通用 release。

`cluster.yaml` 还控制各 rule 的 walltime 和总内存，例如：

| rule | CPU | 默认内存 | 默认时间 |
|---|---:|---:|---:|
| `fastp` | 8 | 16000 MiB | 4 h |
| `bowtie2_align` | 16 | 64000 MiB | 6 h |
| `homer_motif` | 16 | 64000 MiB | 8 h |
| `memechip` | 16 | 128000 MiB | 8 h |

profile 的 `mem_mb` 会以 `sbatch --mem=<数值>M` 传递，是总内存，不是每 CPU 内存。实际分配 CPU 数可能受站点策略或内存/CPU 约束影响，不能只根据 rule 的 threads 推断最终 allocation。MEME-ChIP 内部部分步骤可能只用单线程，但官方综合程序仍以一个整体 rule 运行。

### 8.2 方式 A：登录节点运行 Snakemake driver

只有在集群允许登录节点运行长期 driver 时使用：

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml
```

driver 负责提交 SLURM 子作业；本地调度、DAG 构建和环境创建也可能消耗资源。登录节点上不应执行实际分析程序；不允许长期 driver 时使用方式 B。不要在生产启动时加 `--local`。

### 8.3 方式 B：把 driver 也提交到 SLURM

推荐集群环境使用：

```bash
cd /path/to/project
sbatch \
  --account=YOUR_ACCOUNT \
  --partition=YOUR_PARTITION \
  /path/to/cuttag_pipeline_release/submit.slurm \
  /path/to/project/config.yaml
```

注意：`sbatch` 的 account/partition 只保证 driver 的提交；子作业的 account/partition 仍应在 `profiles/slurm/cluster.yaml` 中配置。`submit.slurm` 不会激活环境；提交前需保证 `python3`、PyYAML 和 Snakemake 7 在继承的 PATH 中。driver 默认 1 CPU、2G 内存、48 小时时限；长流程需按站点要求覆盖这些资源。

### 8.4 限制并行提交数

profile 默认：

```yaml
jobs: 20
```

它表示同时允许的 SLURM 子作业数量，不等于 20 个 CPU。

临时限制为 5 个子作业：

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml \
  -- \
  --jobs 5
```

---

## 9. 监控任务进度

假设 driver Job ID 为 `123456`。

### 9.1 查看自己的所有作业

```bash
squeue -u "$USER" \
  -o "%.18i %.45j %.9T %.10M %.10l %R"
```

常见状态：

| 状态 | 含义 |
|---|---|
| `RUNNING` | 正在运行 |
| `PENDING` | 等待资源或依赖 |
| `Dependency`（Reason，不是 State） | 通常表示 `PENDING` 任务的 SLURM 依赖未满足 |
| `COMPLETED` | 成功完成 |
| `FAILED` | 程序失败 |
| `TIMEOUT` | 超过 walltime |
| `OUT_OF_MEMORY` | 内存不足 |
| `CANCELLED` | 被取消 |

### 9.2 查看 driver

```bash
squeue -j 123456
sacct -j 123456 \
  --format=JobID,JobName,State,ExitCode,Elapsed,Start,End,NodeList
```

### 9.3 查看子任务

```bash
sacct -u "$USER" \
  --starttime YYYY-MM-DD \
  --format=JobID,JobName%45,State,ExitCode,Elapsed,Start,End,NodeList
```

### 9.4 日志位置

```text
project/.snakemake/log/
project/slurm_logs/
project/results/logs/
```

用途：

| 目录 | 内容 |
|---|---|
| `.snakemake/log/` | Snakemake driver 总调度日志 |
| `slurm_logs/` | 每个 SLURM 子作业的 stdout/stderr |
| `results/logs/` | fastp、Bowtie2、MACS3、HOMER、MEME-ChIP 等工具日志 |
| `results/benchmarks/` | 每个 rule 的运行时间和资源记录 |

driver 经 `submit.slurm` 提交时，另有 `cuttag_driver.<JOBID>.out/.err`，默认位于 `sbatch` 的提交工作目录；上面的投递示例先切换到项目目录以避免找不到日志。

实时查看：

```bash
PROJECT=/absolute/path/to/project
LATEST_LOG=$(ls -t "$PROJECT"/.snakemake/log/*.snakemake.log | head -n 1)
tail -F "$LATEST_LOG"
# 或直接检查某个工具的日志
# tail -F "$PROJECT/results/logs/motif_memechip/TF_vs_IgG.log"
```

### 9.5 如何判断真正完成

必须同时满足：

1. driver 为 `COMPLETED`；
2. driver `ExitCode` 为 `0:0`；
3. Snakemake 日志显示 `N of N steps (100%) done`；
4. 没有失败的子任务；
5. 再次 dry-run 显示没有待执行任务。

上述条件验证的是调度和 workflow 完成状态，并不保证每个第三方子程序都没有内部报错；应继续检查关键输出和工具日志。有些综合程序可能生成主 HTML、返回 0，却报告某些子步骤失败。

检查：

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml \
  --dry-run
```

---

## 10. Pipeline 输出内容

完整运行后，`results/` 的目录结构如下（以下以 `project.output_dir: results` 为例）。`<sample>`、`<group>`、`<comparison>` 分别对应 metadata 中的 `sample_id`、`group_id`、`comparison_id`。具体文件是否出现，取决于样本数量、比较关系以及是否启用了 IDR、annotation、motif 和 package。

```text
results/
├── 00.reference/
├── 01.trimmed_reads/
├── 02.mapped/
├── 03.merged/
├── 04.deeptools/
├── 05.peaks/
├── 06.motif/
├── 07.qc/
├── logs/
└── benchmarks/
```

### 10.1 `00.reference/`：参考基因组和索引

```text
00.reference/
├── genome.fa -> 用户配置的参考 FASTA
├── genome.fa.fai
├── genome.sizes
├── bowtie2/
│   └── genome.*.bt2[l]
├── bowtie2_index.ready
└── preparsed/
    ├── genome.fa.200.*
    └── genome.fa.200.ready
```

- `genome.fa` 是指向用户参考 FASTA 的软链接，不是复制的完整 FASTA；
- `genome.fa.fai` 和 `genome.sizes` 由 `samtools faidx` 生成；
- `bowtie2/` 仅在 `reference.bowtie2_index: null`、由 pipeline 在结果目录构建 index 时出现；使用外部预建 index 时，结果目录通常只保留验证完成标志；
- `preparsed/` 是 HOMER 参考背景预解析文件；
- 不要删除 `.fai`、`genome.sizes` 或 `bowtie2_index.ready` 后直接判断结果仍完整。

### 10.2 `01.trimmed_reads/`：fastp 输出

每个 sample 一个目录：

```text
01.trimmed_reads/<sample>/
├── <sample>.R1.fastq.gz
├── <sample>.R2.fastq.gz
├── <sample>.fastp.json
└── <sample>.fastp.html
```

如果 `outputs.keep_trimmed_reads: false`，trimmed R1/R2 可能在下游任务完成后被 Snakemake 自动删除；fastp JSON/HTML 仍用于 QC。`unpaired.*` 和 `failed.fastq.gz` 是中间输出，通常标记为 temporary，不应作为最终分析结果依赖。

### 10.3 `02.mapped/`：sample-level BAM

```text
02.mapped/<sample>/
├── <sample>.filtered.bam
└── <sample>.filtered.bam.bai
```

最终保留的是经过：

- Bowtie2 paired-end 比对；
- MAPQ 过滤；
- 未比对 reads 过滤；
- `reference.chromosomes` 染色体过滤；
- 排序和索引；

之后的 sample-level BAM。未排序 BAM 和临时排序 BAM 不属于最终输出，通常会被 Snakemake 清理。

### 10.4 `03.merged/`：group-level pooled BAM

每个 `group_id` 一个目录：

```text
03.merged/<group>/
├── <group>.bam
├── <group>.bam.bai
├── <group>.flagstat.txt
└── <group>.idxstats.txt
```

- `<group>.bam` 是同一 group 内 biological replicates pooled 后的 BAM；
- 主 MACS3 peak calling 使用 treatment group 和 control group 的 pooled BAM；
- `flagstat.txt` 和 `idxstats.txt` 是基础比对/染色体分布统计，也会进入 MultiQC 输入。

### 10.5 `04.deeptools/`：BigWig 和 QC 图

```text
04.deeptools/
├── sample_bigwig/<sample>.bw
├── group_bigwig/<group>.bw
└── <comparison>/
    ├── correlation.npz
    ├── correlation_matrix.tsv
    ├── correlation_heatmap.png
    ├── PCA.png
    ├── fingerprint.png
    ├── fingerprint_counts.tsv
    ├── TSS_matrix.gz
    ├── TSS_heatmap.png
    └── TSS_profile.png
```

- `sample_bigwig/`：sample-level filtered BAM 的覆盖轨迹；
- `group_bigwig/`：pooled group BAM 的覆盖轨迹；
- `correlation.npz`：`multiBigwigSummary bins` 的矩阵；
- `correlation_heatmap.png` 和 `PCA.png`：基于该矩阵的样本关系图；
- `fingerprint_counts.tsv` 和 `fingerprint.png`：直接从 sample BAM 统计的 fingerprint 结果；
- `TSS_matrix.gz`、`TSS_heatmap.png`、`TSS_profile.png`：以 GTF 的 TSS 为 reference point 的信号矩阵和图形；
- 当前模板四项 deepTools bin size 均为 25 bp，但四类分析的输入和统计方式不同，不能把它们解释为同一种 QC 指标。

### 10.6 `05.peaks/`：MACS3 peaks 和 IDR

主 pooled peak 结果：

```text
05.peaks/<comparison>/
├── <comparison>_peaks.narrowPeak
├── <comparison>_peaks.xls
└── <comparison>_summits.bed
```

- `narrowPeak`：主 peak 区间及 MACS3 统计字段；
- `xls`：MACS3 详细 peak 表；
- `summits.bed`：主 peak summit，供 annotation、HOMER 和 MEME-ChIP 使用。

### 10.7 `05.peaks/<comparison>/idr/`：replicate-level IDR

只有在 treatment group 恰好有两个 biological replicates 且 `idr.mode: auto` 时才会生成：

```text
idr/
├── <replicate1>_peaks.narrowPeak
├── <replicate1>_peaks.xls
├── <replicate1>_summits.bed
├── <replicate2>_peaks.narrowPeak
├── <replicate2>_peaks.xls
├── <replicate2>_summits.bed
├── rep1.sorted.narrowPeak
├── rep2.sorted.narrowPeak
├── idr.tsv
└── idr.tsv.png
```

`idr.tsv` 是 IDR 分析结果，`idr.tsv.png` 是 IDR 诊断图。IDR 结果是重复间一致性/高置信 peak 的 QC 分支，不会自动替代 pooled MACS3 peaks。

### 10.8 annotation 输出

当 `analysis.annotation.enabled: true` 时，在 `05.peaks/<comparison>/` 生成：

```text
<comparison>.annotation.tsv
<comparison>.promoter.bed
<comparison>.promoter_summits.bed
```

- `annotation.tsv` 是 HOMER `annotatePeaks.pl` 输出；
- `promoter.bed` 和 `promoter_summits.bed` 是按绝对 TSS 距离筛选的 summit 区间，当前并非完整 MACS3 peak 边界，也不是参考基因组中所有基因的 promoter 区域；
- 这些文件不是 GO/KEGG enrichment 结果。

### 10.9 `06.motif/`：HOMER 和 MEME-ChIP

```text
06.motif/<comparison>/
├── summits.fasta
├── summits.slop.bed
├── HOMER/
└── MEME-ChIP/
```

`summits.fasta` 和 `summits.slop.bed` 由 MEME-ChIP 的上游 rule 生成；只运行 HOMER 时不要求生成这两个文件。HOMER 和 MEME-ChIP 为独立分支，输出取决于 `motif.methods`。

`HOMER/` 常见文件：

```text
knownResults.txt
knownResults.html
homerResults.html
homerMotifs.all.motifs
nonRedundant.motifs
motifFindingParameters.txt
```

`MEME-ChIP/` 常见文件：

```text
.complete
meme-chip.html
summary.tsv
combined.meme
meme_out/
centrimo_out/
fimo_out_1/
spamo_out_1/
motif_alignment.txt
```

实际子目录和文件会随 MEME Suite 版本、输入数据和运行是否成功而变化。`meme-chip.html`、`summary.tsv`、`combined.meme` 和 `meme_out/` 是优先检查的核心结果；`.complete` 只是 Snakemake 完成标志，不是生物学质量判断。

### 10.10 `07.qc/`：MultiQC

```text
07.qc/
├── multiqc_report.html
└── multiqc_report_data.zip
```

MultiQC rule 当前传入 fastp JSON、group-level flagstat 和 idxstats；报告实际展示的模块取决于安装的 MultiQC 版本及其识别能力。deepTools、IDR 和 motif 图不会自动嵌入本报告。它不会替代对单个 sample、peak 或 motif 日志的检查。

### 10.11 `logs/` 和 `benchmarks/`

```text
logs/
├── prepare_reference.log
├── fastp/<sample>.log
├── bowtie2/<sample>.pipeline.log
├── macs3/<comparison>.log
├── idr/<comparison>.log
├── motif_homer/<comparison>.log
├── motif_memechip/<comparison>.log
└── ...

benchmarks/
├── <rule>.<wildcard>.txt
└── ...
```

- `logs/` 保存工具 stdout/stderr，是排查失败的首要位置；Bowtie2 比对率主要记录在 `bowtie2/<sample>.summary.log`，不要只检查可能为空的 `.pipeline.log`；
- `benchmarks/` 保存 Snakemake benchmark 信息，用于比较耗时和资源；
- `slurm_logs/` 位于项目根目录，不在 `results/` 中，保存 SLURM 子作业 stdout/stderr。

### 10.12 哪些结果优先查看？

建议顺序：

1. `07.qc/multiqc_report.html`：检查 fastp 和基础比对统计；
2. `03.merged/<group>/*.flagstat.txt`、`*.idxstats.txt`：检查 pooled BAM；
3. `04.deeptools/<comparison>/fingerprint.png`、`correlation_heatmap.png`、`PCA.png` 和 TSS 图；
4. `05.peaks/<comparison>/<comparison>_peaks.narrowPeak`：主 peak 集合；
5. `05.peaks/<comparison>/idr/idr.tsv`：两个重复时的 IDR QC；
6. `06.motif/<comparison>/HOMER/` 和 `MEME-ChIP/`：motif 结果；
7. `logs/` 和 `benchmarks/`：异常结果的追踪信息。

pipeline 显示 `100% done` 只说明 workflow 目标完成；最终解释仍需结合比对率、重复率、染色体分布、peak 数量、重复间一致性和实验背景。

## 11. 结果打包

配置：

```yaml
outputs:
  package:
    enabled: true
    archive: cuttag_results.tar.gz
```

生成：

```text
cuttag_results.tar.gz
cuttag_results.tar.gz.manifest.txt
cuttag_results.tar.gz.sha256
```

默认包括：

- config 和 metadata；
- FASTA 的 `.fai` 和 genome sizes；
- BigWig/QC；
- peaks/IDR；
- motif；
- MultiQC；
- logs；
- benchmarks。

默认不包括：

- 原始 FASTQ；
- trimmed FASTQ；
- BAM；
- Bowtie2 index；
- 完整参考 FASTA。

校验：

```bash
# 在原始项目路径可直接检查
sha256sum -c cuttag_results.tar.gz.sha256

# 若已下载到新目录：当前 checksum 文件含服务器绝对路径，改为校验本地文件名
printf '%s  %s\n' "$(awk '{print $1}' cuttag_results.tar.gz.sha256)" \
  cuttag_results.tar.gz | sha256sum -c -

# 查看完整文件清单
less cuttag_results.tar.gz.manifest.txt
```

BigWig 可能占据压缩包的大部分空间。当前打包 rule 仍有兼容性限制：

- 实现固定扫描项目下的 `results/04.deeptools` 等目录；启用打包时请保持 `output_dir: results`，自定义输出目录尚未得到正确支持；
- 打包配置文件固定取项目内的 `config.yaml`，建议保持该文件名；
- metadata 和参考索引结果必须在项目目录之内（原始 FASTQ、原始参考 FASTA 可以位于外部）；
- 会收集扫描目录中遗留的旧结果，并不自动按当前启用的分析方法过滤；从独立 MEME/STREME 版本迁移后，旧目录可能也进入压缩包；
- 不包含项目根目录的 driver 日志或 `slurm_logs/`；这些日志如需归档应另行保存；
- 分析结果压缩包不等同于 pipeline release，也不包含全部复现输入。

---

## 12. 常见故障和处理方法

### 12.1 `snakemake is not available in PATH`

原因：没有激活 Snakemake 7 环境。

```bash
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate snakemake7
snakemake --version
```

### 12.2 Conda 环境创建失败

先确认：

```bash
conda info
mamba --version
```

如果计算节点不能联网：

- 在登录节点提前创建环境；
- 使用已有命名环境；
- 使用共享 Conda 环境目录；
- 将 Conda 包缓存配置到可访问路径。

不要在 rule 内临时安装软件。

### 12.3 `Directory cannot be locked`

先确认没有 driver：

```bash
squeue -u "$USER"
ps -fu "$USER" | grep -E 'snakemake|submit.slurm'
```

确认没有正在运行的 Snakemake 后执行：

```bash
snakemake --unlock \
  --snakefile /path/to/cuttag_pipeline_release/workflow/Snakefile \
  --configfile /path/to/project/config.yaml \
  --directory /path/to/project
```

不要在仍有 driver 运行时删除 `.snakemake/`。

### 12.4 `IncompleteFilesException`

通常说明上一次作业中断或输出不完整。先查看对应 rule 日志，然后使用：

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml \
  -- \
  --rerun-incomplete
```

不要直接删除整个 `results/`。

### 12.5 `PENDING (Dependency)`

这通常不是报错，表示该任务等待上游任务完成：

```bash
scontrol show job JOBID
squeue -j JOBID -o "%.18i %.45j %.9T %R"
```

### 12.6 TIMEOUT 或 OOM

检查：

```bash
sacct -j JOBID \
  --format=JobID,State,ExitCode,Elapsed,MaxRSS,AllocTRES
```

根据 rule 修改 `profiles/slurm/cluster.yaml` 的：

```yaml
time: "HH:MM:SS"
mem_mb: 64000
```

不要只增加全局资源，因为不同 rule 的资源需求差异很大。

### 12.7 HOMER 出现 `CXXABI` 或 `Illegal division by zero`

检查：

```bash
grep -E 'CXXABI_|Illegal division by zero|Might have something wrong' \
  results/logs/motif_homer/*.log
```

当前 pipeline 会设置 HOMER 自带动态库路径，并在日志中检测这些错误。若仍出现：

- 检查 HOMER 环境是否损坏；
- 检查节点系统库差异；
- 优先重建 HOMER 环境；
- 不要把错误日志当作成功结果。

### 12.8 MEME 的 Ghostscript 警告

例如：

```text
Cannot convert EPS file to PNG
```

通常只影响 motif logo 图片转换，不一定影响：

```text
meme-chip.html
combined.meme
summary.tsv
meme.txt
```

先检查核心输出是否非空，不要只依据该警告判定整个 rule 失败。

### 12.9 只想测试某个 rule

`--until` 是 DAG 范围限制，不是“仅运行一个样本”的开关；可能包含多个比较及尚缺失的上游任务。仅测试一个 comparison 时，指定 **绝对目标路径**，并先 dry-run：

```bash
/path/to/cuttag_pipeline_release/run.sh \
  --config /path/to/project/config.yaml \
  --dry-run -- \
  /path/to/project/results/06.motif/TF_vs_IgG/MEME-ChIP/.complete
```

实际运行时去掉 `--dry-run`，且不要加 `--local`，由 profile 提交计算任务。此方式会让 driver 常驻当前 shell；如站点不允许，可写一个独立 driver job 在计算节点调用同一命令。默认 `submit.slurm` 不转发单目标参数。

### 12.10 如何停止任务

先停止 driver，再确认子作业：

```bash
scancel DRIVER_JOB_ID
squeue -u "$USER"
```

如仍有不需要的子作业，再单独取消对应 Job ID。确认没有 Snakemake 进程后，必要时执行 `--unlock`。

---

## 13. 单次运行检查清单

### 运行前

- [ ] Snakemake 是 7.x；
- [ ] Conda 可用，所选 solver 为 Conda 或可用的 Mamba；
- [ ] account 和 partition 已配置；
- [ ] FASTQ 路径存在；
- [ ] R1/R2 正确配对；
- [ ] samples.tsv 一行对应一个 biological replicate；
- [ ] comparisons.tsv 中 treatment/control group 存在；
- [ ] FASTA、GTF、index 和 chromosomes 属于同一 assembly；
- [ ] effective genome size 有明确来源；
- [ ] `run.sh --check-only` 通过；
- [ ] dry-run 没有意外的大规模重算。

### 运行中

- [ ] 保存 driver Job ID；
- [ ] 使用 `squeue` 检查子作业；
- [ ] 使用 `sacct` 检查失败、超时和 OOM；
- [ ] 检查 `.snakemake/log/` 和 `slurm_logs/`；
- [ ] 不要同时启动第二个 driver。

### 运行后

- [ ] driver `COMPLETED`；
- [ ] ExitCode 为 `0:0`；
- [ ] Snakemake 显示 `100% done`；
- [ ] 当前启用的 peak、IDR、motif、QC 输出符合预期；
- [ ] 如开启 MEME-ChIP，`.complete` 存在，且子步骤日志没有未处理错误；
- [ ] MultiQC 报告存在；
- [ ] 如开启打包，结果压缩包 checksum 通过；
- [ ] 保存最终 config、metadata、日志和软件版本。

---

## 14. 科学解释上的重要限制

1. pooled peak 是主要分析结果，IDR 是 replicate-level QC/高置信子集；
2. IDR 自动模式只对恰好两个 treatment replicates 运行；
3. MACS3 使用 control BAM，但 motif 默认不使用 IgG negative FASTA；
4. HOMER 使用参考基因组匹配背景；
5. MEME-ChIP 当前显式关闭 STREME；
6. motif 结果不能直接解释为 treatment-vs-IgG differential motif enrichment；
7. pipeline 成功完成只表示计算任务成功，不代表所有实验数据都通过生物学质量审查。

---

## 15. Pipeline Release 维护原则

- 通用 Pipeline Release 不写死具体实验室 account、partition 和绝对 Conda 路径；
- 服务器专用配置放在项目或 profile 中；
- 不在 pipeline rule 内自动修改 Conda 环境；
- 环境损坏应重建环境；
- 结果和日志不能仅凭目录存在判断成功；
- 重新运行前先检查是否已有 driver 或 SLURM 子作业；
- 修改配置后先做 preflight 和 dry-run，再投递完整流程。


---

## 16. 软件和数据库引用

使用本 pipeline 产生论文、报告或公开结果时，请自行查阅并引用**实际使用的软件和算法的原始文献**，例如 Snakemake、fastp、Bowtie2、SAMtools、BEDTools、deepTools、MACS/MACS3、IDR、HOMER、MEME Suite/MEME-ChIP 和 MultiQC。还应按使用情况引用参考基因组、基因注释及 motif 数据库来源。引用本 pipeline 不能代替对这些原始工作的引用；未运行的可选模块不必仅因模板中出现就写为已使用。

请在 Methods 或补充材料记录实际软件版本、参考 assembly/注释版本、关键参数、metadata 和过滤口径。pipeline 参数模板不是文献引用，也不是对所有数据都适用的分析标准。
