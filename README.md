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

mkdir -p ~/phylogenetic_tree_project/results/genomes
mkdir -p ~/phylogenetic_tree_project/results/genome_downloads/extracted
cd ~/phylogenetic_tree_project
pwd

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


```bash
# Move to the project directory
cd ~/phylogenetic_tree_project

# Create the download directory
mkdir -p results/genome_downloads

# Download KPPR1
datasets download genome accession GCF_000742755.1 \
  --include genome \
  --filename results/genome_downloads/KPPR1.zip

# Download MS6671
datasets download genome accession GCF_001455995.1 \
  --include genome \
  --filename results/genome_downloads/MS6671.zip

# Download NTUH-K2044
datasets download genome accession GCF_000009885.1 \
  --include genome \
  --filename results/genome_downloads/NTUH-K2044.zip

# Download HS11286
datasets download genome accession GCF_000240185.1 \
  --include genome \
  --filename results/genome_downloads/HS11286.zip

# Download MGH 78578
datasets download genome accession GCF_000016305.1 \
  --include genome \
  --filename results/genome_downloads/MGH78578.zip

# Download KPNIH10
datasets download genome accession GCA_000281435.2 \
  --include genome \
  --filename results/genome_downloads/KPNIH10.zip

# Download KPNIH1
datasets download genome accession GCA_000281535.1 \
  --include genome \
  --filename results/genome_downloads/KPNIH1.zip

# Download KPNIH24
datasets download genome accession GCA_000714675.1 \
  --include genome \
  --filename results/genome_downloads/KPNIH24.zip

# Download KPNIH27
datasets download genome accession GCA_000695935.1 \
  --include genome \
  --filename results/genome_downloads/KPNIH27.zip

# Download KPNIH29
datasets download genome accession GCA_000784945.1 \
  --include genome \
  --filename results/genome_downloads/KPNIH29.zip
```
### 1.3 Extract and Rename Genome FASTA Files

NCBI Datasets downloads each assembly as a ZIP archive. Extract the archives and copy the genomic FASTA files into a single directory with consistent strain names.

**Step 1: Create the output directory**

```bash
mkdir -p results/genomes
mkdir -p results/genome_downloads/extracted
```

**Step 2: Extract and rename all ten genome files**

Run these commands from the root of the project directory.

```bash
cd ~/phylogenetic_tree_project

# Process all ten genome archives
for strain in KPPR1 MS6671 NTUH-K2044 HS11286 \
              MGH78578 KPNIH10 KPNIH1 KPNIH24 KPNIH27 KPNIH29
do
    echo "Processing ${strain}..."

    archive="results/genome_downloads/${strain}.zip"
    extract_dir="results/genome_downloads/extracted/${strain}"

    # Check that the downloaded archive exists
    if [ ! -f "$archive" ]; then
        echo "ERROR: Archive not found: $archive"
        continue
    fi

    # Extract the archive
    mkdir -p "$extract_dir"
    unzip -oq "$archive" -d "$extract_dir"

    # Locate the genomic FASTA file
    fasta=$(find "$extract_dir/ncbi_dataset/data" \
        -type f -name '*_genomic.fna' -print -quit)

    # Copy it to the standard strain filename
    if [ -n "$fasta" ]; then
        cp "$fasta" "results/genomes/${strain}.fna"
        echo "Saved: results/genomes/${strain}.fna"
    else
        echo "ERROR: Genomic FASTA not found for ${strain}"
    fi
done
```

**Step 3: Verify the renamed files**

```bash
ls -lh results/genomes/
```

The directory should contain these ten FASTA files:

```text
KPPR1.fna
MS6671.fna
NTUH-K2044.fna
HS11286.fna
MGH78578.fna
KPNIH10.fna
KPNIH1.fna
KPNIH24.fna
KPNIH27.fna
KPNIH29.fna
```

**Step 4: Inspect the sequence headers**

```bash
for file in results/genomes/*.fna; do
    echo "=== $file ==="
    grep '^>' "$file" | head -n 5
done
```

The headers help you identify the sequences included in each assembly.

**Important quality-control note:** An assembly FASTA may contain both chromosome and plasmid sequences. The directory is therefore named `genomes`, not `chromosomes_only`. Do not assume these files contain only chromosomes. Verify their sequence headers and accession records before using them for chromosome-based phylogenetic analysis.

# Clean up zip files
rm *.zip
```

#### Option B: Manual Browser Download (Zero Code)
1. Navigate to the [NCBI Assembly Database](https://www.ncbi.nlm.nih.gov/assembly/).
2. Search each accession number above (e.g., `GCF_000240185.1`).
3. Click **Download Assembly** -> **Genomic FASTA (.fna)**.
4. Move and rename the downloaded files into `~/phylogenetic_tree_project/results/chromosomes_only/` as `KPPR1.fna`, `MS6671.fna`, `NTUH-K2044.fna`, and `HS11286.fna`.

### 1.4 Verify Downloaded Files
cd ~/phylogenetic_tree_project

ls -lh results/genomes/

# Count the downloaded FASTA files
find results/genomes -maxdepth 1 -type f -name '*.fna' | wc -l
✅ **Checkpoint**: You should see four `.fna` files, each typically between 5.0 MB and 5.6 MB in size.

---

## SECTION 2: Inspect Genome Assemblies

Before running Snippy, check that the ten FASTA files exist and inspect their sequence counts and lengths.

### 2.1 Calculate Assembly Summary Statistics

```bash
cd ~/phylogenetic_tree_project

seqkit stats results/genomes/*.fna
```

**Checkpoint:** Confirm that ten files are listed. Check the sequence count and total length for each assembly. Genome assemblies may contain multiple sequences, including plasmids, so total lengths and sequence counts can differ.

### 2.2 Inspect FASTA Sequence Headers

```bash
for file in results/genomes/*.fna; do
    echo "=== $file ==="
    grep '^>' "$file" | head -n 5
done
```

FASTA headers begin with `>`. Inspect them to identify the sequences included in each assembly.

**Important:** These checks confirm file presence and basic sequence structure. They do not prove species purity, assembly completeness, or the absence of contamination. Verify the accession and chromosome/plasmid status before analysis.

---

##  SECTION 3: Reference-Based Variant Calling (Snippy)

###  The Mental Model
Think of the reference genome (**KPPR1**) as the original baseline edition of a dictionary. **Snippy** reads every sentence in your sample genomes (**MS6671**, **NTUH-K2044**, **HS11286**) and notes every spelling change, deletion, or insertion relative to the baseline.

3.1 Run Snippy for Each Sample
Run Snippy separately for each of the nine sample genomes, using KPPR1 as the reference.

cd ~/phylogenetic_tree_project

# Sample 1: MS6671
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/MS6671.fna \
  --outdir results/snippy_MS6671 --cpus 2

# Sample 2: NTUH-K2044
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/NTUH-K2044.fna \
  --outdir results/snippy_NTUH-K2044 --cpus 2

# Sample 3: HS11286
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/HS11286.fna \
  --outdir results/snippy_HS11286 --cpus 2

# Sample 4: MGH 78578
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/MGH78578.fna \
  --outdir results/snippy_MGH78578 --cpus 2

# Sample 5: KPNIH10
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/KPNIH10.fna \
  --outdir results/snippy_KPNIH10 --cpus 2

# Sample 6: KPNIH1
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/KPNIH1.fna \
  --outdir results/snippy_KPNIH1 --cpus 2

# Sample 7: KPNIH24
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/KPNIH24.fna \
  --outdir results/snippy_KPNIH24 --cpus 2

# Sample 8: KPNIH27
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/KPNIH27.fna \
  --outdir results/snippy_KPNIH27 --cpus 2

# Sample 9: KPNIH29
snippy --ref results/genomes/KPPR1.fna \
  --ctgs results/genomes/KPNIH29.fna \
  --outdir results/snippy_KPNIH29 --cpus 2
What this does: Each sample is compared with KPPR1, and each analysis is saved in its own output directory. The nine output directories will later be combined to create the core alignment.

Commit with this message:

Update Snippy commands to use the genomes director

### 3.2 Verify Snippy Outputs

After all nine Snippy commands have completed successfully, check that each sample has its own output directory.

```bash
cd ~/phylogenetic_tree_project

for strain in MS6671 NTUH-K2044 HS11286 MGH78578 \
              KPNIH10 KPNIH1 KPNIH24 KPNIH27 KPNIH29
do
    echo "=== $strain ==="

    if [ -d "results/snippy_${strain}" ]; then
        ls "results/snippy_${strain}" | head
    else
        echo "ERROR: Output directory missing for $strain"
    fi
done
```

**Checkpoint:** Confirm that all nine output directories exist. Inspect the Snippy logs and expected output files for errors before proceeding to `snippy-core`.

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
## SECTION 4: Build the Core-Genome Alignment

### 4.1 Run Snippy-core

After all nine Snippy analyses have completed successfully, combine their results into a shared core-genome alignment.

```bash
cd ~/phylogenetic_tree_project

snippy-core \
  --ref results/genomes/KPPR1.fna \
  --prefix results/core10 \
  results/snippy_MS6671 \
  results/snippy_NTUH-K2044 \
  results/snippy_HS11286 \
  results/snippy_MGH78578 \
  results/snippy_KPNIH10 \
  results/snippy_KPNIH1 \
  results/snippy_KPNIH24 \
  results/snippy_KPNIH27 \
  results/snippy_KPNIH29
```

### 4.2 Inspect the Alignment Files

The expected alignment files are:

- `results/core10.full.aln` — the full core alignment.
- `results/core10.aln` — the SNP-only alignment.

Check the sequence counts and alignment lengths:

```bash
grep -c '^>' results/core10.full.aln
seqkit stats results/core10.full.aln results/core10.aln
```

**Checkpoint:** The full alignment should contain ten sequences: the reference plus nine samples. Investigate any missing sequences, unexpected gaps, or excessive ambiguous bases before proceeding.

**Important:** A successful sequence count alone does not establish that the alignment is biologically appropriate. Review the Snippy logs and alignment quality before building the tree.

---

## SECTION 5: Clean the Alignment and Build the Phylogenetic Tree

### 5.1 Clean the Full Core Alignment

Remove problematic alignment positions according to the Snippy cleaning utility's documented behavior.

```bash
cd ~/phylogenetic_tree_project

snippy-clean_full_aln < results/core10.full.aln \
  > results/core10.clean.full.aln
```

Check that the cleaned alignment exists:

```bash
ls -lh results/core10.clean.full.aln
```

### 5.2 Optional Advanced Analysis: Recombination Filtering

*Klebsiella pneumoniae* can undergo homologous recombination, which may affect inferred relationships. For a more advanced analysis, use Gubbins to identify likely recombinant regions.

Gubbins must be installed, and the input alignment must meet its requirements.

```bash
run_gubbins.py \
  --prefix results/gubbins_core10 \
  --threads 2 \
  results/core10.clean.full.aln
```

Inspect the files generated by your installed Gubbins version before selecting the alignment for tree inference. Do not assume every version produces identical filenames.

### 5.3 Build a Preliminary Maximum-Likelihood Tree

For a basic analysis using the cleaned full core alignment, run:

```bash
iqtree2 \
  -s results/core10.clean.full.aln \
  -m MFP \
  -B 1000 \
  -T 2 \
  --prefix results/klebsiella_tree
```

This uses ModelFinder to select a suitable substitution model and performs 1,000 ultrafast bootstrap replicates.

**Expected outputs**

- `results/klebsiella_tree.treefile` — inferred maximum-likelihood tree.
- `results/klebsiella_tree.iqtree` — analysis report.
- `results/klebsiella_tree.contree` — consensus tree with branch support.

Check the log for successful completion, the selected model, and any warnings.

**Important:** This is a preliminary teaching workflow. For a final scientific analysis, assess alignment quality, missing data, recombination, reference suitability, and the effect of analysis choices. If you instead use a SNP-only alignment, assess ascertainment bias and select the appropriate IQ-TREE model options.


### 6. Build the Maximum-Likelihood Phylogenetic Tree

The alignment produced above contains variable SNP sites only. Use an ascertainment-bias correction when inferring the tree from this SNP-only alignment.

```bash
iqtree2 \
  -s results/core10.recombination_filtered_snps.aln \
  -m MFP+ASC \
  -B 1000 \
  -T 2 \
  --prefix results/klebsiella_tree
```

**Expected outputs**

- `results/klebsiella_tree.treefile` — maximum-likelihood tree.
- `results/klebsiella_tree.iqtree` — analysis report, including the selected model.
- `results/klebsiella_tree.contree` — consensus tree with branch support.

Check the IQ-TREE log to confirm that the analysis completed successfully and that the model and ascertainment correction were applied as intended.

**Important:** This workflow requires `snippy-clean_full_aln`, `run_gubbins.py`, `snp-sites`, and `iqtree2` to be installed and available in the active environment. Do not interpret the tree until the alignment quality and genome labels have been checked.
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

## SECTION 6: Format and Visualize the Tree in iTOL

### 6.1 Prepare Readable Labels

First, inspect the inferred tree to confirm its actual tip labels.

```bash
cd ~/phylogenetic_tree_project

cat results/klebsiella_tree.treefile
```

Create a visualization copy with simpler labels:

```bash
sed 's/snippy_//g; s/Reference/KPPR1/g' \
  results/klebsiella_tree.treefile \
  > results/klebsiella_tree_iTOL.tree
```

This assumes the tree uses `Reference` for the reference tip and `snippy_` prefixes for sample tips. Check the tree first; if the actual labels differ, adjust the replacement accordingly.

### 6.2 Upload the Tree to iTOL

1. Open [iTOL](https://itol.embl.de/).
2. Upload `results/klebsiella_tree_iTOL.tree`.
3. Check that all ten expected tips are present.
4. Inspect the topology, branch lengths, and support values.
5. Customize the display and export an SVG or PNG if needed.

If you use Windows Subsystem for Linux (WSL), you can open the results directory in Windows Explorer:

```bash
explorer.exe "$(wslpath -w "$PWD/results")"
```

**Interpretation reminder:** Genetic similarity alone does not prove direct transmission between patients or establish the direction of transmission.

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

cd ~/phylogenetic_tree_project

{
  echo "Execution Date: $(date -Is)"
  echo "---------------------------"
  echo "Snippy Version:"
  snippy --version
  echo "IQ-TREE 2 Version:"
  iqtree2 --version
  echo "SeqKit Version:"
  seqkit version
  echo "NCBI Datasets Version:"
  datasets version
} > results/software_versions.txt

cat results/software_versions.txt
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
