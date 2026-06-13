# GNNPhase

GNNPhase 是一个面向 read-based haplotype assembly / genome phasing 的图神经网络流程。项目把 reads 作为图节点，把 read-read overlap 作为边，并将 SNP 覆盖信息投影到节点和边上作为可训练特征；模型输出每条 read 属于 haplotype 1 的概率，再基于 read-level score 对每个 SNP 位点投票，生成 phased VCF，并进一步用 HapCUT2 的统计工具计算 switch error/rate、mismatch 或 Hamming error、N50/AN50 等指标。

---


## 1. Environment

推荐使用 conda/mamba 新建环境：

```bash
mamba create -n gnnphase python=3.10 -y
mamba activate gnnphase
pip install -r requirements.txt
```

如果使用 GPU，请按照本机 CUDA 版本单独安装匹配的 PyTorch 和 DGL。安装完成后确认：

```bash
python - <<'PY'
import torch, dgl, pysam
print("torch:", torch.__version__)
print("cuda:", torch.cuda.is_available())
print("dgl:", dgl.__version__)
print("pysam:", pysam.__version__)
PY
```

---

## 2. External tools

除 Python 包外，还需要下列命令行工具。请确保它们都在 `PATH` 中。

| Tool | Purpose | Check |
|---|---|---|
| minimap2 | reads-to-reference/contig alignment; read-read overlap PAF | `minimap2 --version` |
| samtools | BAM sort/index/stats/downsample | `samtools --version` |
| bcftools | VCF 过滤、标准化、consensus | `bcftools --version` |
| PBSIM3 | 模拟 ONT reads，用于 synthetic data | `pbsim --version` |
| PECAT | ONT read error correction 和第一阶段 contig assembly | `pecat.pl` |
| Clair3 | 无真实 SNP VCF 时进行 long-read SNP calling | `run_clair3.sh --help` |
| HapCUT2 | phasing 评估指标：switch rate、mismatch rate、N50、AN50 | `HAPCUT2 --help` |
| hifiasm | 分 haplotype reads 后进行组装 | `hifiasm --version` |
| seqtk, pigz/bgzip, tabix | 可选；用于 read subset、压缩和索引 | `seqtk`, `bgzip`, `tabix` |

常见安装方式：

```bash
mamba install -c conda-forge -c bioconda minimap2 samtools bcftools seqtk hifiasm -y
```

PECAT、Clair3、HapCUT2 和 PBSIM3 建议按各自官方仓库安装；Clair3 通常单独放在 `clair3` conda 环境或 Docker/Singularity 环境中。

---

## 3. Input files

### 3.1 Required input

真实数据或模拟数据最终都需要准备：

```text
reads.fasta / reads.fastq(.gz)        # 原始或校正后的 reads
reference.fa 或 contigs.fa            # 参考基因组或 PECAT 第一阶段 contigs
sample.vcf.gz                         # SNP VCF；真实 VCF 优先，没有则由 Clair3 生成
```

### 3.2 Real SNP VCF available

如果已有真实的 SNP VCF，例如 GIAB truth VCF、Platinum Genomes VCF 或其他可靠的 heterozygous SNP callset，则直接使用它作为 read 节点 SNP 特征来源：

```bash
VCF=sample.real_het_snp.vcf.gz
```

建议只使用 biallelic heterozygous SNPs，并保持 contig 命名与 BAM/reference 一致，例如 `chr1` 与 `1` 不要混用。

### 3.3 No real SNP VCF: call SNPs with PECAT + minimap2 + Clair3

如果没有真实 SNP 位点，需要先从 reads 中提取 SNP。推荐流程如下：

#### Step 1. PECAT error correction

```bash
pecat.pl correct pecat.cfg
```

PECAT 官方流程中，纠错后的 reads 通常位于：

```text
${PROJECT}/1-correct/corrected_reads.fasta
```

#### Step 2. PECAT first-stage assembly to contigs

```bash
pecat.pl assemble pecat.cfg
```

第一阶段组装得到的 contigs 通常位于：

```text
${PROJECT}/3-assemble/primary.fasta
```

在本流程中，`primary.fasta` 可作为临时 reference/contig，用于后续 Clair3 SNP calling。

#### Step 3. Align corrected reads to PECAT contigs with minimap2

ONT reads 示例：

```bash
THREADS=32
CORR_READS=${PROJECT}/1-correct/corrected_reads.fasta
CONTIGS=${PROJECT}/3-assemble/primary.fasta
BAM=clair3_input.corrected_reads_to_contigs.bam

minimap2 -ax map-ont -t ${THREADS} ${CONTIGS} ${CORR_READS} \
  | samtools sort -@ ${THREADS} -o ${BAM}

samtools index -@ ${THREADS} ${BAM}
```

PacBio/HiFi reads 可将 preset 改为 `map-pb` 或 `map-hifi`，取决于 reads 类型：

```bash
minimap2 -ax map-hifi -t ${THREADS} ${CONTIGS} ${CORR_READS} \
  | samtools sort -@ ${THREADS} -o ${BAM}
samtools index -@ ${THREADS} ${BAM}
```

#### Step 4. Run Clair3 to generate VCF

ONT 模型示例：

```bash
THREADS=32
CLAIR3_MODEL=${CONDA_PREFIX}/bin/models/ont_r10_dorado_sup_5khz
CLAIR3_OUT=clair3_out

run_clair3.sh \
  --bam_fn=${BAM} \
  --ref_fn=${CONTIGS} \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=${CLAIR3_MODEL} \
  --output=${CLAIR3_OUT} \
  --include_all_ctgs
```

Clair3 的最终 VCF 通常为：

```text
${CLAIR3_OUT}/merge_output.vcf.gz
```

将它作为后续 SNP 特征来源：

```bash
VCF=${CLAIR3_OUT}/merge_output.vcf.gz
```

> 注意：如果 Clair3 VCF 是以 PECAT contig 坐标生成，那么用于 `build_read_graph()` 的 BAM 和 VCF 必须在同一 contig/reference 坐标系下。若使用 GRCh38 truth VCF，则 BAM 也应是 reads-to-GRCh38 的 BAM。

---

## 4. Prepare reads and alignments

### 4.1 Rename reads to numeric IDs

`snp.py` 中提供 `rename_fasta_with_ids()`，会生成：

```text
sample_id.fasta
sample.id2name
```

建议将 read header 统一改成数字 ID，便于 DGL 节点编号和 BAM/PAF/VCF 投票时的 read name 对齐。

示例 Python：

```python
from snp import rename_fasta_with_ids

rename_fasta_with_ids(
    input_fasta="sample.reads.fasta",
    output_fasta="sample_id.fasta",
    id2name_file="sample.id2name",
)
```

### 4.2 Read-read overlap PAF

ONT reads：

```bash
minimap2 -x ava-ont -t 32 sample_id.fasta sample_id.fasta > sample.paf
```

PacBio/HiFi reads：

```bash
minimap2 -x ava-pb -t 32 sample_id.fasta sample_id.fasta > sample.paf
```

`snp.py` 中的 `run_read_to_read_alignment()` 默认使用近似 PacBio preset 的 `ava-pb`，并额外设置：

```bash
minimap2 -x ava-pb -t 32 -g 3000 -w 30 -k 19 -m 100 -r 500 sample_id.fasta sample_id.fasta > sample.paf
```

### 4.3 Reads-to-reference or reads-to-contig BAM

ONT reads：

```bash
minimap2 -ax map-ont -t 32 reference_or_contigs.fa sample_id.fasta \
  | samtools sort -@ 32 -o sample.bam
samtools index -@ 32 sample.bam
```

HiFi reads：

```bash
minimap2 -ax map-hifi -t 32 reference_or_contigs.fa sample_id.fasta \
  | samtools sort -@ 32 -o sample.bam
samtools index -@ 32 sample.bam
```

---

## 5. Build SNP features and DGL graph

### 5.1 SNP feature definitions

`build_read_graph()` 从 VCF、BAM 和 PAF 中提取 SNP 相关特征：

Node SNP features:

```text
[n_het_cov, ref_frac, alt_frac, conflict_frac]
```

Edge SNP features:

```text
[n_shared_het, n_same_allele, n_diff_allele, allele_consistency]
```

完整节点特征由 read length 加 SNP 特征组成；完整边特征由 overlap length、alignment similarity 加 SNP edge features 组成。构图后会对节点和边特征做标准化，并保存为 DGL `.bin` 文件。

### 5.2 Example graph construction script

建议新建 `scripts/build_graph.py`，内容示例：

```python
import argparse
import dgl
from snp import (
    build_read_graph,
    collect_read_reference_positions,
    save_feats,
    load_feats,
    build_graph,
    preprocess_and_add_sg_factors,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--id2name", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--prefix", default="sample")
    parser.add_argument("--min_mapq", type=int, default=10)
    parser.add_argument("--min_baseq", type=int, default=13)
    args = parser.parse_args()

    node_dict, edge_dict = build_read_graph(
        truth_vcf_path=args.vcf,
        bam_path=args.bam,
        paf_path=args.paf,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
    )

    node_json = f"{args.prefix}.node_feats.json"
    edge_json = f"{args.prefix}.edge_feats.json"
    save_feats(node_dict, edge_dict, node_json, edge_json)
    node_var_feats, edge_var_feats = load_feats(node_json, edge_json)

    read_ref_pos = collect_read_reference_positions(
        bam_path=args.bam,
        min_mapq=args.min_mapq,
    )

    g, _, _, _ = build_graph(
        id2name_file=args.id2name,
        fasta_file=args.fasta,
        paf_files=[args.paf],
        node_var_feats=node_var_feats,
        edge_var_feats=edge_var_feats,
        hits_dict=None,
        read_ref_pos=read_ref_pos,
        use_positional_prune=True,
        use_backbone_prune=False,
    )

    g = preprocess_and_add_sg_factors(g)
    dgl.save_graphs(args.out, g)
    print(f"Saved graph to {args.out}")


if __name__ == "__main__":
    main()
```

运行：

```bash
python scripts/build_graph.py \
  --vcf ${VCF} \
  --bam sample.bam \
  --paf sample.paf \
  --id2name sample.id2name \
  --fasta sample_id.fasta \
  --out data/train_graphs/sample_graph.bin \
  --prefix sample
```

将训练、验证、测试图分别放入：

```text
data/train_graphs/*.bin
data/val_graphs/*.bin
data/test_graphs/*.bin
```

---

## 6. Training

在 `hyperparameters.py` 中设置关键参数：

```python
save_train_graph_path = "data/train_graphs"
save_val_graph_path   = "data/val_graphs"
save_test_graph_path  = "data/test_graphs"
save_garph_path       = "results/test_graphs"   # 原代码变量名为 save_garph_path

test_mode = False
stage = "pretrain"       # or "finetune"
device = "cuda:0"        # or "cpu"
```

启动训练：

```bash
python train.py
```

模型会读取 DGL `.bin` graph，使用 `GraphGatedGCNModel` 训练 read-level haplotype assignment。训练阶段使用 BCE loss，并结合 edge consistency loss 约束共享 SNP 边上的 haplotype 一致性。

Checkpoint 默认保存到：

```text
checkpoints/best_model_<dataset_name>_pretrain.pth
checkpoints/info_<dataset_name>_pretrain.json
```

如果进行 fine-tuning：

```python
stage = "finetune"
pretrained_path = "checkpoints/best_model_<dataset_name>_pretrain.pth"
finetune_phase = 1   # or 2
```

---

## 7. Testing and read-level score output

在 `hyperparameters.py` 中设置：

```python
test_mode = True
save_test_graph_path = "data/test_graphs"
save_garph_path = "results/test_graphs"
```

运行：

```bash
python train.py
```

测试完成后，保存的 graph 中会包含：

```text
graph.ndata["score"]       # 每条 read 属于 haplotype 1 的概率
graph.ndata["pred_valid"]  # 是否有有效预测
graph.ndata["prob"]        # 概率副本，若代码开启
graph.ndata["hard_label_04"]
```

---

## 8. Export reads by score

将 read-level score 转换成 hap0/hap1/ambiguous read list：

```bash
python dump_reads_by_score.py \
  --graph results/test_graphs/test_0_graph.bin \
  --id2name sample.id2name \
  --threshold 0.6 \
  --prefix results/GNNPhase \
  --use_original_name
```

输出：

```text
results/GNNPhase.hap0.reads.txt
results/GNNPhase.hap1.reads.txt
results/GNNPhase.ambiguous.reads.txt
```

阈值规则：

```text
score >= 0.6  -> hap1
score <= 0.4  -> hap0
0.4 < score < 0.6 -> ambiguous
```

---

## 9. Read-level score to SNP-level phased VCF

`evaluate_saved_graph.py` 会把 read-level score 投影到每个 heterozygous SNP 位点上，并根据 hap0/hap1 read 支持的 REF/ALT 数量投票，输出 phased VCF。

```bash
python evaluate_saved_graph.py \
  --graph results/test_graphs/test_0_graph.bin \
  --bam sample.bam \
  --vcf ${VCF} \
  --out results/GNNPhase.phased.vcf \
  --id2name sample.id2name \
  --threshold 0.6 \
  --min_mapq 10 \
  --min_baseq 13 \
  --min_total_reads 3 \
  --min_reads_per_hap 1
```

如果整体 hap0/hap1 方向相反，可加：

```bash
--flip
```

输出文件：

```text
results/GNNPhase.phased.vcf
```

建议压缩和索引：

```bash
bgzip -f results/GNNPhase.phased.vcf
tabix -f -p vcf results/GNNPhase.phased.vcf.gz
```

---

## 10. Evaluation: switch rate/error, mismatch/Hamming error, N50/AN50

本项目建议与参考论文一致，使用 HapCUT2 的 `calculate_haplotype_statistics.py` 计算 phasing metrics。该工具会比较 predicted phased VCF 与 truth phased VCF，输出 switch rate、mismatch rate、flat error rate、phased count、N50、AN50 等指标。

### 10.1 Prepare per-chromosome VCFs

HapCUT2 统计脚本通常按染色体读取 VCF。示例：

```bash
mkdir -p metrics/vcf_by_chr

PRED_VCF=results/GNNPhase.phased.vcf.gz
TRUTH_VCF=truth.phased.vcf.gz

for chr in chr{1..22}; do
  bcftools view -r ${chr} -Ov -o metrics/vcf_by_chr/pred.${chr}.vcf ${PRED_VCF}
  bcftools view -r ${chr} -Ov -o metrics/vcf_by_chr/truth.${chr}.vcf ${TRUTH_VCF}
done
```

### 10.2 Run HapCUT2 statistics

```bash
HAPCUT2_DIR=/path/to/HapCUT2

python ${HAPCUT2_DIR}/utilities/calculate_haplotype_statistics.py \
  -v1 metrics/vcf_by_chr/pred.chr1.vcf metrics/vcf_by_chr/pred.chr2.vcf metrics/vcf_by_chr/pred.chr3.vcf metrics/vcf_by_chr/pred.chr4.vcf metrics/vcf_by_chr/pred.chr5.vcf metrics/vcf_by_chr/pred.chr6.vcf metrics/vcf_by_chr/pred.chr7.vcf metrics/vcf_by_chr/pred.chr8.vcf metrics/vcf_by_chr/pred.chr9.vcf metrics/vcf_by_chr/pred.chr10.vcf metrics/vcf_by_chr/pred.chr11.vcf metrics/vcf_by_chr/pred.chr12.vcf metrics/vcf_by_chr/pred.chr13.vcf metrics/vcf_by_chr/pred.chr14.vcf metrics/vcf_by_chr/pred.chr15.vcf metrics/vcf_by_chr/pred.chr16.vcf metrics/vcf_by_chr/pred.chr17.vcf metrics/vcf_by_chr/pred.chr18.vcf metrics/vcf_by_chr/pred.chr19.vcf metrics/vcf_by_chr/pred.chr20.vcf metrics/vcf_by_chr/pred.chr21.vcf metrics/vcf_by_chr/pred.chr22.vcf \
  -v2 metrics/vcf_by_chr/truth.chr1.vcf metrics/vcf_by_chr/truth.chr2.vcf metrics/vcf_by_chr/truth.chr3.vcf metrics/vcf_by_chr/truth.chr4.vcf metrics/vcf_by_chr/truth.chr5.vcf metrics/vcf_by_chr/truth.chr6.vcf metrics/vcf_by_chr/truth.chr7.vcf metrics/vcf_by_chr/truth.chr8.vcf metrics/vcf_by_chr/truth.chr9.vcf metrics/vcf_by_chr/truth.chr10.vcf metrics/vcf_by_chr/truth.chr11.vcf metrics/vcf_by_chr/truth.chr12.vcf metrics/vcf_by_chr/truth.chr13.vcf metrics/vcf_by_chr/truth.chr14.vcf metrics/vcf_by_chr/truth.chr15.vcf metrics/vcf_by_chr/truth.chr16.vcf metrics/vcf_by_chr/truth.chr17.vcf metrics/vcf_by_chr/truth.chr18.vcf metrics/vcf_by_chr/truth.chr19.vcf metrics/vcf_by_chr/truth.chr20.vcf metrics/vcf_by_chr/truth.chr21.vcf metrics/vcf_by_chr/truth.chr22.vcf \
  | tee metrics/GNNPhase.hapcut2_statistics.txt
```

指标解释：

- `switch rate/error`：相邻 heterozygous SNP 的 phase 关系与 truth 不一致的比例。
- `mismatch rate` 或 `Hamming error`：短 switch 或连续错误造成的单点 haplotype mismatch 比例。
- `N50`：phased block 按 genomic span 统计的 N50。
- `AN50`：adjusted N50，会考虑未分相 SNP 对 block 连续性的影响。

---

## 11. Haplotype-specific assembly with hifiasm

模型分相后，可以根据 read score 将 reads 划分为 hap0/hap1，然后分别组装。

### 11.1 Extract haplotype-specific reads

如果原始 reads 是 FASTQ：

```bash
seqtk subseq sample.reads.fastq.gz results/GNNPhase.hap0.reads.txt > results/hap0.reads.fastq
seqtk subseq sample.reads.fastq.gz results/GNNPhase.hap1.reads.txt > results/hap1.reads.fastq
pigz -f results/hap0.reads.fastq
pigz -f results/hap1.reads.fastq
```

如果原始 reads 是 FASTA：

```bash
seqtk subseq sample.reads.fasta results/GNNPhase.hap0.reads.txt > results/hap0.reads.fasta
seqtk subseq sample.reads.fasta results/GNNPhase.hap1.reads.txt > results/hap1.reads.fasta
```

### 11.2 hifiasm assembly commands

HiFi reads：

```bash
THREADS=32
mkdir -p asm

hifiasm -o asm/hap0 -t ${THREADS} results/hap0.reads.fastq.gz
hifiasm -o asm/hap1 -t ${THREADS} results/hap1.reads.fastq.gz
```

ONT reads：

```bash
THREADS=32
mkdir -p asm

hifiasm -o asm/hap0 --ont -t ${THREADS} results/hap0.reads.fastq.gz
hifiasm -o asm/hap1 --ont -t ${THREADS} results/hap1.reads.fastq.gz
```

将 GFA contigs 转为 FASTA：

```bash
awk '/^S/{print ">"$2"\n"$3}' asm/hap0.bp.p_ctg.gfa > asm/hap0.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' asm/hap1.bp.p_ctg.gfa > asm/hap1.p_ctg.fa
```

---

## 12. Synthetic data generation

`simulate.py` 中合成数据流程包括：

1. 从 1000 Genomes phased panel 中提取目标样本的 biallelic heterozygous SNPs；
2. 用 `bcftools consensus -H 1` 和 `-H 2` 生成 hap1/hap2 reference；
3. 用 PBSIM3 分别从 hap1/hap2 模拟 reads；
4. 合并两个 haplotype 的 reads；
5. 用 minimap2/samtools 生成 BAM；
6. 用 minimap2 all-vs-all 生成 read-read PAF；
7. 可选：用 samtools 下采样到 30x、15x、10x、5x。

运行前需要修改 `simulate.py` 顶部配置：

```python
BASE_DIR = Path("/path/to/project")
REF_FA = BASE_DIR / "G38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
PANEL_DIR = BASE_DIR / "G38"
OUT_ROOT = BASE_DIR / "gnnphase_ont_synth"
SAMPLES = ["HG00142"]
PBSIM = Path("/path/to/pbsim3/src/pbsim")
QSHMM_MODEL = Path("/path/to/pbsim3/data/QSHMM-ONT-HQ.model")
```

然后运行：

```bash
python simulate.py
```

---

## 13. Coverage estimation and downsampling

估计 BAM 覆盖度并生成下采样命令：

```bash
python estimate_bam_coverage.py \
  --bam sample.bam \
  --mode autosome \
  --threads 16 \
  --targets 5,10,15,30
```

脚本会输出类似：

```bash
samtools view -@ 16 -b --subsample-seed 1 --subsample <fraction> -o sample.5x.bam sample.bam
samtools index -@ 16 sample.5x.bam
```

---

## 14. End-to-end example

下面是一个真实 ONT 数据的最小化流程示例。真实路径需按本地环境修改。

```bash
# 0. Activate environment
mamba activate gnnphase

# 1. Optional: if no real SNP VCF is available, generate VCF using PECAT + Clair3
pecat.pl correct pecat.cfg
pecat.pl assemble pecat.cfg

THREADS=32
CORR_READS=${PROJECT}/1-correct/corrected_reads.fasta
CONTIGS=${PROJECT}/3-assemble/primary.fasta
BAM=sample.corrected_to_contigs.bam
CLAIR3_OUT=clair3_out
CLAIR3_MODEL=${CONDA_PREFIX}/bin/models/ont_r10_dorado_sup_5khz

minimap2 -ax map-ont -t ${THREADS} ${CONTIGS} ${CORR_READS} \
  | samtools sort -@ ${THREADS} -o ${BAM}
samtools index -@ ${THREADS} ${BAM}

run_clair3.sh \
  --bam_fn=${BAM} \
  --ref_fn=${CONTIGS} \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=${CLAIR3_MODEL} \
  --output=${CLAIR3_OUT} \
  --include_all_ctgs

VCF=${CLAIR3_OUT}/merge_output.vcf.gz

# 2. Rename reads and build PAF/BAM
python - <<'PY'
from snp import rename_fasta_with_ids
rename_fasta_with_ids("sample.reads.fasta", "sample_id.fasta", "sample.id2name")
PY

minimap2 -x ava-ont -t ${THREADS} sample_id.fasta sample_id.fasta > sample.paf
minimap2 -ax map-ont -t ${THREADS} ${CONTIGS} sample_id.fasta \
  | samtools sort -@ ${THREADS} -o sample.bam
samtools index -@ ${THREADS} sample.bam

# 3. Build graph
python scripts/build_graph.py \
  --vcf ${VCF} \
  --bam sample.bam \
  --paf sample.paf \
  --id2name sample.id2name \
  --fasta sample_id.fasta \
  --out data/test_graphs/sample_graph.bin \
  --prefix sample

# 4. Test with trained model
python train.py

# 5. Convert read score to phased VCF
python evaluate_saved_graph.py \
  --graph results/test_graphs/test_0_graph.bin \
  --bam sample.bam \
  --vcf ${VCF} \
  --out results/GNNPhase.phased.vcf \
  --id2name sample.id2name \
  --threshold 0.6

# 6. Dump reads and assemble each haplotype
python dump_reads_by_score.py \
  --graph results/test_graphs/test_0_graph.bin \
  --id2name sample.id2name \
  --threshold 0.6 \
  --prefix results/GNNPhase \
  --use_original_name

seqtk subseq sample.reads.fastq.gz results/GNNPhase.hap0.reads.txt > results/hap0.reads.fastq
seqtk subseq sample.reads.fastq.gz results/GNNPhase.hap1.reads.txt > results/hap1.reads.fastq
pigz -f results/hap0.reads.fastq results/hap1.reads.fastq

hifiasm -o asm/hap0 --ont -t ${THREADS} results/hap0.reads.fastq.gz
hifiasm -o asm/hap1 --ont -t ${THREADS} results/hap1.reads.fastq.gz
```

---

## 15. Notes and troubleshooting

1. **Coordinate consistency**：VCF、BAM 和 reference/contig 必须来自同一坐标系。GRCh38 truth VCF 需要 reads-to-GRCh38 BAM；Clair3-on-PECAT-contigs VCF 需要 reads-to-PECAT-contigs BAM。
2. **Read name consistency**：`sample_id.fasta`、PAF、BAM、DGL graph 和 `id2name` 必须对应。若 BAM 使用原始 read name，`evaluate_saved_graph.py` 可通过 `--id2name` 兼容。
3. **Contig naming**：`chr1` 和 `1` 不一致会导致 pileup/VCF fetch 失败。
4. **Large graphs**：大图训练可能内存较高，建议调小 `num_parts_metis_train`、`batch_size_train`，或在构图阶段过滤低质量 PAF 边。
5. **Ambiguous reads**：默认 `0.4 < score < 0.6` 的 reads 不参与 SNP 投票，也不进入 hap0/hap1 read list；可以调整 `--threshold`。
6. **Haplotype flip**：二倍体分相存在整体 hap0/hap1 标签互换问题。评估或输出 VCF 时若发现整体方向相反，可使用 `--flip`。
7. **Current code cleanup**：建议把 `snp.py` 底部样例流程拆到单独脚本，否则 import `snp.py` 时会直接运行样例数据。
