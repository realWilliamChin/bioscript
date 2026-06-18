# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview
This is a comprehensive bioinformatics analysis script repository containing tools for genomic, transcriptomic, metabolomic, and metagenomic data processing, analysis, and visualization. All scripts are designed for high-throughput sequencing data analysis workflows.

## Core Architecture
The repository is organized by functional domains:

- **CommonTools/** - Shared utility libraries and generic processing tools
  - Submodules: Bam, Blast, Fasta, Fastq, alignment, gff, geneid conversion, multithreading utilities
  - Key shared components: `multithreads_task_runner.py` (parallel task execution), `load_input.py` (universal file I/O), `data_check.py` (data cleaning)

- **Analysis/** - Domain-specific analysis pipelines
  - MultiDESeq (differential gene expression), WGCNA (co-expression network), enrichment analysis, time-series analysis, etc.

- **Plot/** - Visualization modules
  - All common bioinformatics plots: PCA, heatmaps, volcano plots, Venn diagrams, correlation plots, GO/KEGG classification plots, etc.

- **kns_annotation/** - Genome annotation pipeline
  - End-to-end annotation workflows: KEGG (KAAS integration), SwissProt, NR, GO, PPI network prediction, summary statistics

- **Specialized Analysis Modules**:
  - MetaboliteAnalysis/ - Metabolomics data processing and enrichment
  - metagenomics/ - 16S/metagenomics analysis
  - transcriptome/ - Transcriptome analysis pipelines
  - reseq/ - Variant calling (GATK) and resequencing workflows
  - chloroplast/ - Chloroplast genome assembly and analysis
  - DenovoGenome/ - De novo genome assembly tools

- **Workflow/** - Snakemake pipelines for complex workflows
- **Rscript/** - R utility functions and microeco analysis
- **lib/** - Shared reference databases and static files
- **Unused/** - Deprecated scripts (do not use for new development)

## Coding Conventions
### Python Scripts
- Use `loguru` for logging, include type hints in function definitions
- Docstrings and comments are in Chinese
- CLI tools use `argparse` for argument parsing
- Parallel processing uses `multithreads_task_runner.py` for both CPU and I/O bound tasks
- Import order: standard libraries, third-party libraries, local modules

### R Scripts
- Use `optparse` for command line argument parsing
- Start with `suppressPackageStartupMessages(library())` for all imports
- Output tabular results as Excel files using `openxlsx`
- Plots are generated in both PDF and PNG formats by default
- Comments and error messages are in Chinese

### Commit Message Format
Based on recent commit history, use a concise Chinese one-line commit message with the changed module or script path up front:

- Preferred pattern: `[TYPE] Module/script.py 中文描述`
- Use uppercase bracketed types:
  - `[NEW]` - 新增脚本或全新的独立工具
  - `[FEAT]` - 为现有脚本、流程或工具新增功能
  - `[UPDATE]` - 优化现有功能、参数、输出或兼容性
  - `[FIX]` - 修复错误、数据类型问题或结果异常
- Put the primary changed file/path immediately after the type, for example: `[NEW] CommonTools/Bam/bam_stat.py 新增BAM文件批量统计工具...`
- Prefer one script per commit to keep history clear and reviewable.
- Multiple scripts may be committed together only when their functions are strongly related, such as one pipeline plus its helper script, paired Python/R plotting outputs, or coordinated changes required for the same feature/fix.
- Descriptions should be in Chinese and clearly summarize the user-visible change; for larger tools, mention key capabilities separated by Chinese punctuation.
- Avoid older inconsistent forms seen in history such as lowercase `[feat]`, `[update]`, or conventional prefixes like `feat:` / `refactor:` for new commits.

## Common Commands
### Running Scripts
Most scripts are standalone executables:
```bash
# Python scripts
python CommonTools/Bam/bam_stat.py --help
python kns_annotation/scripts/kegg_annotation.py --input sequences.fasta --species plant --output ./annotation

# R scripts
Rscript Plot/pca.r --data_file fpkm_matrix.txt --samples samples.txt --output pca_result
```

### Snakemake Workflows
For pipelines with Snakefile:
```bash
# Run with 8 cores
snakemake -s kns_annotation/Snakefile --cores 8 --configfile config.yaml
```

### Parallel Processing
Use the shared multi-thread runner for custom parallel tasks:
```python
from CommonTools.multithreads_task_runner import run_multithreads_tasks
results = run_multithreads_tasks(your_function, list_of_args_tuples, max_workers=8)
```

## Key Notes
- All file paths in scripts use Unix-style forward slashes `/` even on Windows
- Scripts expect standard bioinformatics file formats (FASTQ, BAM, GFF, VCF, etc.)
- Many scripts have built-in file format detection and automatic handling of compressed files
- Chinese characters are widely used in comments, docstrings, and output file headers - ensure proper encoding (UTF-8) when processing outputs

## Development Workflow
1. New scripts should be placed in the appropriate module directory based on functionality
2. Reuse existing utilities from CommonTools/ instead of rewriting similar functionality
3. For complex multi-step pipelines, use Snakemake and place in WorkFlow/ or module-specific Snakefiles
4. Update commit messages following the standard format
5. Deprecated scripts should be moved to Unused/ rather than deleted
