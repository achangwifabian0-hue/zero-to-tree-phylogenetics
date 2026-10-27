# Zero to Tree: Bacterial Phylogenetics for Absolute Beginners
*A Hands-On, Fully Reproducible Guide Using Klebsiella pneumoniae, Ubuntu/WSL, Snippy, IQ-TREE 2, and iTOL*

---

### Why This Course is Built Differently (Pedagogical Innovations)
To make this repository accessible to beginners on GitHub without requiring external assistance, this curriculum incorporates five teaching innovations:
1. **The "Cookbook" Mental Model**: Every technical step pairs a biological analogy (e.g., comparing SNP calling to identifying misprints in copies of a reference book) with the underlying computational logic.
2. **Download from Scratch Protocol**: Resolves the missing genome acquisition step with both automated terminal commands (NCBI Datasets CLI / `curl`) and a point-and-click browser fallback.
3. **"Self-Check" Checkpoints**: Every command includes an *Expected Output* block so learners can verify their intermediate files before running downstream analyses.
4. **Visual "Gotchas" & Interpretation Traps**: Clear warnings to prevent beginner errors, such as misinterpreting branch lengths as calendar dates or assuming reference self-comparisons are independent biological isolates.
5. **Interactive GitHub Features**: Formatted with GitHub markdown callouts (`> [!NOTE]`, `> [!TIP]`, `> [!WARNING]`), collapsible troubleshooting accordions, and clean copy-paste command blocks.

---

```markdown
#  Zero to Tree: Bacterial Phylogenetics from Scratch

[![Platform: Ubuntu/WSL](https://img.shields.io/badge/Platform-Ubuntu%20%7C%20WSL-orange)](https://ubuntu.com/)
[![Workflow: Snippy + IQ--TREE 2 + iTOL](https://img.shields.io/badge/Workflow-Snippy%20%7C%20IQ--TREE2%20%7C%20iTOL-blue)](https://github.com/tseemann/snippy)
[![Level: Absolute Beginner](https://img.shields.io/badge/Level-Absolute%20Beginner-brightgreen)](#)

Welcome! This course guides you from raw bacterial genome sequences to a publication-ready, fully annotated phylogenetic tree. No prior bioinformatics experience is assumed.

---

## 🎯 Learning Objectives
By the end of this practical guide, you will be able to:
- Navigate files and directories comfortably using the Linux command line (Ubuntu or WSL).
- Download and inspect bacterial genome assemblies (`.fna` FASTA files) from public databases.
- Perform reference-based variant calling using **Snippy**.
- Construct a core-genome alignment using **Snippy-core**.
- Infer a maximum-likelihood phylogenetic tree with **IQ-TREE 2**.
- Visualize, style, and annotate your tree interactively in **iTOL**.
- Evaluate results using a scientific quality-control checklist.

---

## Pipeline Overview
```text
[NCBI Public Database]
         │  (Step 1: Download assemblies)
         ▼
[KPPR1 (Ref), MS6671, NTUH-K2044, HS11286] (.fna)
         │  (Step 2: Quality inspection via SeqKit)
         ▼
   [Snippy v4.6+] ──► Calls SNPs against KPPR1 reference
         │
         ▼
   [Snippy-core]  ──► Generates core-genome alignment (.aln)
         │
         ▼
   [IQ-TREE 2]    ──► ModelFinder (MFP) + 1,000 Ultrafast Bootstraps
         │
         ▼
 [klebsiella_tree.treefile] (Newick format)
         │  (Clean label prefixes via sed)
         ▼
      [iTOL]      ──► Interactive visualization & publication export (SVG/PNG)
```

---

##  SECTION 0: Prerequisites & Environment Setup

This workflow runs inside an **Ubuntu Linux terminal** (native Linux, macOS, or Windows Subsystem for Linux - WSL).

### 0.1 Activate Conda Environment
Ensure Conda (Miniconda or Mambaforge) is installed, then activate your phylogenetic environment:
```bash
conda activate snippy-env
```

### 0.2 Confirm Software Installation
Verify that your core analysis tools are installed and callable:
```bash
snippy --version
iqtree2 --version
seqkit version
```
> [!TIP]
> If any command returns `command not found`, install them via bioconda:
> `conda create -n snippy-env -c conda-forge -c bioconda snippy iqtree seqkit ncbi-datasets-cli`

---

## 📁 SECTION 1: Project Setup & Genome Acquisition

### 1.1 Create Your Project Directories
Create a clean directory structure to organize raw assemblies and downstream results:
```bash
mkdir -p ~/phylogenetic_tree_project/results/chromosomes_only
cd ~/phylogenetic_tree_project
pwd
```
### 1.2 Download Genomes from NCBI

This tutorial uses ten *Klebsiella pneumoniae* genome assemblies: one reference genome and nine additional strains. Assembly accessions identify the exact genome records used in the analysis.

| Role | Strain | Assembly accession |
|---|---|---|
| Reference | KPPR1 | `GCF_000742755.1` |
| Sample 1 | MS6671 | `GCF_001455995.1` |
| Sample 2 | NTUH-K2044 | `GCF_000009885.1` |
| Sample 3 | HS11286 | `GCF_000240185.1` |
| Sample 4 | MGH 78578 | `GCF_000016305.1` |
| Sample 5 | KPNIH10 | `GCA_000281435.2` |
| Sample 6 | KPNIH1 | `GCA_000281535.2` |
| Sample 7 | KPNIH24 | `GCA_000714675.1` |
| Sample 8 | KPNIH27 | `GCA_000695935.1` |
| Sample 9 | KPNIH29 | `GCA_000784945.1` |

**Important:** The `GCF_` and `GCA_` prefixes refer to different NCBI assembly records. Both can identify genome assemblies, but the exact accession must be retained when downloading and documenting each genome.

Before analysis, always verify the assembly identity and sequence contents for every strain. This dataset is intended for this teaching workflow, not as a substitute for selecting epidemiologically appropriate isolates for a research study.

#### Option A: One-Liner Download via NCBI Datasets CLI (Recommended)
```bash
cd ~/phylogenetic_tree_project/results/chromosomes_only

# Download genomes using NCBI Datasets CLI
datasets download genome accession GCF_000240185.1 --filename KPPR1.zip
datasets download genome accession GCF_001010365.1 --filename MS6671.zip
datasets download genome accession GCF_000025705.1 --filename NTUH-K2044.zip
datasets download genome accession GCF_000222155.1 --filename HS11286.zip

# Unzip and rename files
unzip -p KPPR1.zip ncbi_dataset/data/GCF_000240185.1/*.fna > KPPR1.fna
unzip -p MS6671.zip ncbi_dataset/data/GCF_001010365.1/*.fna > MS6671.fna
unzip -p NTUH-K2044.zip ncbi_dataset/data/GCF_000025705.1/*.fna > NTUH-K2044.fna
unzip -p HS11286.zip ncbi_dataset/data/GCF_000222155.1/*.fna > HS11286.fna

# Clean up zip files
rm *.zip
```

#### Option B: Manual Browser Download (Zero Code)
1. Navigate to the [NCBI Assembly Database](https://www.ncbi.nlm.nih.gov/assembly/).
2. Search each accession number above (e.g., `GCF_000240185.1`).
3. Click **Download Assembly** -> **Genomic FASTA (.fna)**.
4. Move and rename the downloaded files into `~/phylogenetic_tree_project/results/chromosomes_only/` as `KPPR1.fna`, `MS6671.fna`, `NTUH-K2044.fna`, and `HS11286.fna`.

### 1.3 Verify Downloaded Files
```bash
ls -lh ~/phylogenetic_tree_project/results/chromosomes_only/
```
✅ **Checkpoint**: You should see four `.fna` files, each typically between 5.0 MB and 5.6 MB in size.

---

## SECTION 2: Inspect Genome Assemblies

Before running calculations, verify file integrity and review sequence structure.

### 2.1 Calculate Assembly Summary Statistics
Run `seqkit stats` to inspect sequence counts and genome lengths:
```bash
seqkit stats ~/phylogenetic_tree_project/results/chromosomes_only/*.fna
```
*Expected Output:*
Each genome should exhibit ~5.2–5.5 million base pairs (bp), characteristic of *Klebsiella pneumoniae*.

### 2.2 Inspect FASTA Sequence Headers
Confirm sequence formatting:

# Move to the project directory first
cd ~/phylogenetic_tree_project

# Display the first three FASTA headers from each genome
for file in results/chromosomes_only/*.fna; do
    echo "=== $file ==="
    grep '^>' "$file" | head -n 3
done
```
> [!NOTE]
> Sequence lines starting with `>` are FASTA headers. These summary checks confirm file presence and basic metrics, but they do not prove species purity, complete assembly, or absence of contamination.

---

##  SECTION 3: Reference-Based Variant Calling (Snippy)

###  The Mental Model
Think of the reference genome (**KPPR1**) as the original baseline edition of a dictionary. **Snippy** reads every sentence in your sample genomes (**MS6671**, **NTUH-K2044**, **HS11286**) and notes every spelling change, deletion, or insertion relative to the baseline.

### 3.1 Execute Snippy for Each Sample
Run Snippy separately for each sample assembly against the reference:

```bash
cd ~/phylogenetic_tree_project

# Sample 1: MS6671
snippy \
  --ref results/chromosomes_only/KPPR1.fna \
  --ctgs results/chromosomes_only/MS6671.fna \
  --outdir results/snippy_MS6671 \
  --cpus 2

# Sample 2: NTUH-K2044
snippy \
  --ref results/chromosomes_only/KPPR1.fna \
  --ctgs results/chromosomes_only/NTUH-K2044.fna \
  --outdir results/snippy_NTUH-K2044 \
  --cpus 2

# Sample 3: HS11286
snippy \
  --ref results/chromosomes_only/KPPR1.fna \
  --ctgs results/chromosomes_only/HS11286.fna \
  --outdir results/snippy_HS11286 \
  --cpus 2
```

### 3.2 Verify Snippy Outputs
Confirm that the output directories were created:
```bash
ls -lh results/snippy_MS6671/
ls -lh results/snippy_NTUH-K2044/
ls -lh results/snippy_HS11286/
```
✅ **Checkpoint**: Inside each directory, locate `snps.tab`, `snps.vcf`, and `snps.aligned.fa`.

---

##  SECTION 4: Build the Core-Genome Alignment (Snippy-core)

###  The Mental Model
Each Snippy run identified variants for a single isolate against the reference. **Snippy-core** stacks these analyses together to find genomic regions shared across **all** samples (the *core genome*) and compiles them into a single multiple-sequence alignment.

### 4.1 Run Snippy-core
```bash
snippy-core \
  --ref results/chromosomes_only/KPPR1.fna \
  --prefix results/core3 \
  results/snippy_MS6671 \
  results/snippy_NTUH-K2044 \
  results/snippy_HS11286
```

### 4.2 Inspect the Alignment Files
Snippy-core outputs two key alignment files:
1. `results/core3.full.aln`: Full multiple-sequence alignment containing all conserved positions.
2. `results/core3.aln`: SNP-only alignment containing strictly variable sites.

Verify sequence representation and lengths using `seqkit`:
```bash
seqkit fx2tab -n -l results/core3.full.aln
seqkit fx2tab -n -l results/core3.aln
```

> [!WARNING]
> Inspect alignment outputs for missing samples, excess ambiguous bases (`N`), or uneven sequence lengths. If using a SNP-only alignment for tree building, evaluate whether an ascertainment bias correction (e.g., `+ASC` in IQ-TREE) is required to account for the absence of invariant sites.

---

##  SECTION 5: Construct the Tree with IQ-TREE 2

###  The Mental Model
IQ-TREE 2 evaluates evolutionary models using **ModelFinder** and tests thousands of branching patterns to find the tree that best explains your alignment using Maximum Likelihood. It resamples positions 1,000 times (**Ultrafast Bootstrap**) to determine confidence scores for each branch.

### 5.1 Run Tree Inference
```bash
iqtree2 \
  -s results/core3.full.aln \
  -m MFP \
  -B 1000 \
  -T 2 \
  --prefix results/klebsiella_tree
```

### 5.2 Command Parameters Explained
| Flag | Function | Why It Matters |
| :--- | :--- | :--- |
| `-s results/core3.full.aln` | Input alignment | Points to the full core alignment. |
| `-m MFP` | ModelFinder Plus | Tests substitution models and selects the best fit based on BIC/AIC scores. |
| `-B 1000` | Ultrafast Bootstrap | Performs 1,000 replicates to calculate branch confidence values. |
| `-T 2` | Threads | Utilizes 2 CPU cores for parallel processing. |
| `--prefix results/...` | Output prefix | Directs and names all generated result files. |

### 5.3 Inspect Your Inferred Tree
Display the raw Newick tree string:
```bash
cat results/klebsiella_tree.treefile
```
*Expected Output:* A nested parentheses string specifying relationships, branch lengths, and support scores (e.g., `((MS6671:0.0001..., NTUH-K2044:...)...);`).

---

##  SECTION 6: Format and Visualize in iTOL

### 6.1 Clean Label Prefixes for Display
By default, sample names in the tree file may retain the folder prefix `snippy_`. Strip this prefix using `sed` to create a clean display tree without altering topology or branch lengths:
```bash
sed 's/snippy_//g' results/klebsiella_tree.treefile > results/klebsiella_tree_iTOL.tree
cat results/klebsiella_tree_iTOL.tree
```

### 6.2 Upload to iTOL
1. Open the [iTOL (Interactive Tree Of Life)](https://itol.embl.de/) portal in your browser.
2. Click **Upload tree files** in the top navigation bar.
3. If using WSL, quickly open your results folder in Windows Explorer:
   ```bash
   explorer.exe "\\((wslpath -w "\\)PWD/results")"
   ```
4. Drag and drop `klebsiella_tree_iTOL.tree` into the iTOL upload box.

### 6.3 Recommended Display Settings for Beginners
- **Tree Structure**: Choose **Rectangular** layout for clear cladogram/phylogram inspection.
- **Branch Lengths**: Toggle to **Real** so branch distances reflect estimated mutations per site.
- **Internal Node Support**: Turn on **Display bootstrap values** as text or symbols to show reproducibility scores.
- **Exporting**: Click the **Export** tab and select **SVG** (vector graphic for presentations/papers) or **PNG** (high-resolution image).

---

##  SECTION 7: Scientific Interpretation & Common Traps

| Tree Feature | Biological Meaning | ⚠️ Common Beginner Trap to Avoid |
| :--- | :--- | :--- |
| **Topology** | Which isolates share a more recent common ancestor. | Grouping together does not prove direct transmission or direct contact between patients. |
| **Branch Length** | Estimated nucleotide substitutions per site under the fitted model. | Branch lengths do not represent calendar time (weeks, months, years). |
| **Bootstrap Support** | How consistently a clade is recovered when re-sampling alignment sites. | A bootstrap of 100% does not guarantee the sequencing data is error-free or that the biological conclusion is correct. |
| **The Reference** | Serves as the coordinate system for variant discovery. | Never treat a self-comparison of the reference as an independent clinical isolate. |

---

##  SECTION 8: Scientific Quality-Control Checklist

Before presenting or publishing results, review this checklist:
- [ ] Assembly sources and accession numbers are verified.
- [ ] Assembly completeness and potential contamination have been assessed.
- [ ] The reference genome is appropriate for the target lineage.
- [ ] Snippy execution logs have been reviewed for dropped reads or errors.
- [ ] Unaligned regions and ambiguous bases (`N`) were documented.
- [ ] The evolutionary model selected by ModelFinder was recorded.
- [ ] Tree relationships were interpreted alongside epidemiological and sampling context.
- [ ] Software versions and command parameters were recorded.

---

## SECTION 9: Reproducibility Log

Execute this script to capture exact software versions in your results directory:
```bash
{
  echo "Execution Date: \$(date -Is)"
  echo "---------------------------"
  echo "Snippy Version:"
  snippy --version
  echo "IQ-TREE 2 Version:"
  iqtree2 --version
  echo "SeqKit Version:"
  seqkit version
} > results/software_versions.txt

cat results/software_versions.txt
```

---

##  SECTION 10: Troubleshooting FAQ

<details>
<summary><b>1. Error: <code>snippy: command not found</code></b></summary>

Ensure your Conda environment is active:
```bash
conda activate snippy-env
which snippy
```
</details>

<details>
<summary><b>2. Error: <code>iqtree2: command not found</code></b></summary>

Install or verify IQ-TREE 2 inside your active environment:
```bash
conda install -c bioconda iqtree
iqtree2 --version
```
</details>

<details>
<summary><b>3. Error: <code>No such file or directory</code> when loading assemblies</b></summary>

Check your current path and confirm file names:
```bash
pwd
ls -lh results/chromosomes_only/
```
Verify that the file names match the commands (`KPPR1.fna`, `MS6671.fna`, `NTUH-K2044.fna`, `HS11286.fna`).
</details>

## Official Software Documentation
• Snippy: https://github.com/tseemann/snippy
• IQ-TREE 2: https://www.iqtree.org/
• SeqKit: https://bioinf.shenwei.me/seqkit/
• iTOL: https://itol.embl.de/
• AliView: https://ormbunkar.se/aliview/

Course Completion
You have completed the main workflow when you can:
1. Inspect the genome assemblies.
2. Run reference-based variant calling.
3. Generate and assess a core-genome alignment.
4. Infer a preliminary tree with IQ-TREE 2.
5. Display the tree in iTOL.
6. Explain the basic meaning of topology, branch lengths, and support values.
7. Identify the quality-control checks needed before drawing scientific conclusions.

Final reminder: This course is intended for learning and reproducible bioinformatics practice. Do not use a preliminary tree as the sole basis for clinical decisions, outbreak declarations, or transmission claims.

## This course is design from knowledge gained on phylogenetic tree design and analysis during our genomics course, African Stars Fellowship, Institute Pasteur De Dakar and CARE

Thanks for learning with us

## Producer
- Fabian Achangwi
## Course Facilitators
- Ndeye Awa NDIAYE
- Abdoul Wahab Soumare
