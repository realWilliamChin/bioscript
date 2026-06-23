# 生信软件 PATH 集合 (单行 PATH 类工具，已去重)
# Aspera
PATH=$PATH:/home/train/.aspera/connect/bin/

# 测序质控 / reads 处理
PATH=$PATH:/home/data/opt/biosoft/FastQC/
PATH=$PATH:/home/data/opt/biosoft/FastUniq/source/
PATH=/opt/biosoft/fastp/:$PATH
PATH=/opt/biosoft/sickle/:$PATH
PATH=/opt/biosoft/TrimGalore-0.6.0:$PATH
PATH=/home/data/opt/biosoft/cutadapt:$PATH

# 纠错 / 组装
PATH=$PATH:/home/data/opt/biosoft/BLESS/
alias bless='mkdir -p kmc/bin/; cp /home/data/opt/biosoft/BLESS/kmc/bin/kmc kmc/bin/; bless'
PATH=$PATH:/home/data/opt/biosoft/ALLPATHS-LG/bin/
PATH=$PATH:/home/data/opt/biosoft/lordec-src_0.9/
PATH=$PATH:/home/data/opt/biosoft/idba-1.1.3/bin/
PATH=$PATH:/home/data/opt/biosoft/SOAPdenovo2-r241/
PATH=$PATH:/home/data/opt/biosoft/DBG2OLC
PATH=$PATH:/home/data/opt/biosoft/wtdbg-2.4_x64_linux/
PATH=$PATH:/home/data/opt/biosoft/Platanus_allee_v2.0.2_Linux_x86_64
PATH=$PATH:/home/data/opt/biosoft/quickmerge-0.3
PATH=$PATH:/home/data/opt/biosoft/GapFiller_v1-11_linux-x86_64/
PATH=$PATH:/home/data/opt/biosoft/GapCloser-v1.12-r6/
PATH=$PATH:/home/data/opt/biosoft/SSPACE-STANDARD-3.0_linux-x86_64/
PATH=$PATH:/home/data/opt/biosoft/SPAdes-3.13.0-Linux/bin/
PATH=$PATH:/opt/biosoft/canu-1.8/Linux-amd64/bin/
PATH=/home/data/opt/biosoft/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$PATH

# 比对 / 序列工具
PATH=$PATH:/home/data/opt/biosoft/minimap2-2.17_x64-linux/
PATH=$PATH:/home/data/opt/biosoft/mummer-4.0.0beta2/bin/
PATH=$PATH:/home/data/opt/biosoft/bowtie2-2.3.5-linux-x86_64/
PATH=$PATH:/home/data/opt/biosoft/bowtie-1.2.2-linux-x86_64/
PATH=$PATH:/home/data/opt/biosoft/bwa-0.7.17/
PATH=$PATH:/home/data/opt/biosoft/samtools-1.9/bin/
PATH=$PATH:/home/data/opt/biosoft/bcftools-1.9/bin
PATH=$PATH:/home/data/opt/biosoft/tophat-2.1.1.Linux_x86_64/
PATH=$PATH:/home/data/opt/biosoft/hisat2-2.1.0/
PATH=$PATH:/home/data/opt/biosoft/blat/
PATH=$PATH:/home/data/opt/biosoft/fasta-36.3.8g/bin/
PATH=$PATH:/home/data/opt/biosoft/gmap-2017-11-15/bin/
PATH=$PATH:/home/data/opt/biosoft/blast-2.2.26/bin/
PATH=/home/data/opt/biosoft/ncbi-blast-2.9.0+/bin/:$PATH
PATH=$PATH:/home/data/opt/biosoft/diamond
PATH=/home/data/opt/biosoft/seqtk:$PATH
PATH=/home/data/opt/biosoft/seqkit:$PATH
PATH=/home/data/opt/biosoft/bioawk:$PATH
PATH=$PATH:/home/data/opt/biosoft/subread-2.0.1-Linux-x86_64/bin/

# HMMER / 结构域
PATH=/home/data/opt/biosoft/hmmer-3.2.1/bin/:$PATH
PATH=/home/data/opt/biosoft/RpsbProc-x64-linux:$PATH

# 重复序列
PATH=$PATH:/home/data/opt/biosoft/RepeatMasker
PATH=$PATH:/home/data/opt/biosoft/RepeatModeler-open-1.0.11

# ncRNA / tRNA / 引物
PATH=$PATH:/home/data/opt/biosoft/infernal-1.1.5/bin
PATH=$PATH:/home/data/opt/biosoft/Misa_Primer3/
PATH=$PATH:/home/data/opt/biosoft/primer3-2.4.0/src/
PATH=$PATH:/home/data/opt/biosoft/tRNAscan-SE-2.0/bin/
PATH=/home/data/opt/biosoft/rnammer-1.2:$PATH

# 转录组 / 定量
PATH=$PATH:/home/data/opt/biosoft/Trinity-v2.8.5/
PATH=$PATH:/home/data/opt/biosoft/RSEM-1.3.1/
PATH=$PATH:/home/data/opt/biosoft/kallisto
PATH=$PATH:/home/data/opt/biosoft/TransDecoder-v5.5.0/
PATH=$PATH:/home/data/opt/biosoft/cufflinks-2.2.1.Linux_x86_64/
PATH=$PATH:/home/data/opt/biosoft/stringtie-1.3.5.Linux_x86_64/
PATH=$PATH:/home/data/opt/biosoft/salmon-latest_linux_x86_64/bin/

# 基因预测 / 注释
PATH=$PATH:/home/data/opt/biosoft/PASApipeline-v2.3.3/bin/
PATH=$PATH:/home/data/opt/biosoft/gmes_linux_64/
PATH=$PATH:/home/data/opt/biosoft/exonerate-2.2.0-x86_64/bin/
PATH=$PATH:/home/data/opt/biosoft/snap
PATH=$PATH:/home/data/opt/biosoft/geta-2.4.3/bin/
PATH=$PATH:/home/data/opt/biosoft/SpliceGrapher-0.2.7/scripts/
PATH=$PATH:/home/data/opt/biosoft/isolasso/bin/
PATH=$PATH:/home/data/opt/biosoft/gms2_linux_64//
PATH=$PATH:/home/data/opt/biosoft/sequin/
PATH=$PATH:/home/data/opt/biosoft/tbl2asn/
PATH=$PATH:/opt/biosoft/maker/bin/
PATH=/home/data/opt/biosoft/MetaGeneMark_linux_64/mgm:$PATH
PATH=$PATH:/home/data/opt/biosoft/gmes_linux_64/

# 比较基因组 / 聚类
PATH=$PATH:/home/data/opt/biosoft/mcl-14-137/bin/
PATH=$PATH:/home/data/opt/biosoft/orthomclSoftware-v2.0.9/bin/
PATH=$PATH:/home/data/opt/biosoft/mafft/bin/
PATH=$PATH:/home/data/opt/biosoft/MCScanX/
PATH=$PATH:/home/data/opt/biosoft/mauve_2.4.0/
PATH=/home/data/opt/biosoft/cdhit-4.6.2:$PATH

# 可视化
PATH=$PATH:/home/data/opt/biosoft/IGV_Linux_2.5.3/
PATH=$PATH:/home/data/opt/biosoft/circos-0.69-6/bin/

# Web 服务 (biosoft 下的 lighttpd)
PATH=$PATH:/home/data/opt/biosoft/lighttpd-1.4.54/sbin/

# 细胞器 / 叶绿体
PATH=$PATH:/home/data/opt/biosoft/GetOrganelle-master/
PATH=$PATH:/home/data/opt/biosoft/GetOrganelle-master/Utilities/
PATH=/opt/biosoft/PGA:$PATH

# 群体遗传 / 变异
PATH=$PATH:/home/data/opt/biosoft/vcftools-0.1.16/src/cpp
PATH=$PATH:/home/data/opt/biosoft/plink_1.9
PATH=$PATH:/home/data/opt/biosoft/admixture_linux-1.3.0

# 宏基因组 / 物种分类
PATH=/home/data/opt/biosoft/taxonkit:$PATH
PATH=/home/data/opt/biosoft/csvtk:$PATH
PATH=/home/data/opt/biosoft/FLASH-1.2.11-Linux-x86_64:$PATH
PATH=/home/data/opt/biosoft/centrifuge:$PATH

# 流程 / 工具链
PATH=$PATH:$HOME/local/app/edirect
PATH=$PATH:/home/train/smrtlink/smrtcmds/bin/
PATH=/home/data/opt/biosoft/bpipe-0.9.9.9/bin:$PATH
PATH=/home/train/miniconda3/envs/test/bin:$PATH
PATH=/home/train/miniconda3/envs/test2/bin:$PATH

# htslib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/biosoft/htslib-1.9/lib/
C_INCLUDE_PATH=$C_INCLUDE_PATH:/opt/biosoft/htslib-1.9/include/
PATH=$PATH:/home/data/opt/biosoft/htslib-1.9/bin/

# GeneWise
PATH=$PATH:/home/data/opt/biosoft/wise2.4.1/src/bin/
export WISECONFIGDIR=/home/data/opt/biosoft/wise2.4.1/wisecfg/

# GATK / lncRNA
PATH=$PATH:/home/data/opt/biosoft/gatk-4.4.0.0/
export CPC_HOME=/home/data/opt/biosoft/CPC2-beta

# SRA toolkit
export SRA_TOOLKIT_PATH=/home/data/opt/biosoft/sratoolkit.3.1.1-centos_linux64/bin
PATH=$PATH:$SRA_TOOLKIT_PATH
