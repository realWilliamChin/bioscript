# ============================================================
# general.sh —— 历史遗留 / 待确认配置留存档
# ============================================================
# 原 general.sh 中所有“生效”的环境变量已按 sys/bio/other 分类拆分到
# 同目录的子文件夹中，并由 env.sh 统一加载 (zshrc 现 source env.sh)。
#
# 本文件仅保留原 general.sh 中【被注释掉的死代码】【失效/被覆盖的配置】
# 【标注 ??? 存疑】的内容，供你后续逐项确认是否彻底删除或恢复。
# 默认本文件不被 env.sh 加载，source 它不会改变环境 (内容均为注释)。
# ============================================================

# ---- 被注释的旧版本 / 弃用配置 ----
# PATH=$PATH:/home/data/opt/biosoft/ncbi-blast-2.6.0+/bin/
# PATH=$PATH:/home/data/opt/biosoft/infernal-1.1.2/bin
# PATH=$PATH:/home/data/opt/biosoft/tRNAscan-SE-1.3.1/bin/
# export PERL5LIB=$PERL5LIB:/home/data/opt/biosoft/tRNAscan-SE-1.3.1/bin/  # ????
# GATK4.1.2 对 MNPs 支持不太好，修改为 4.4.0.0 版本 - 2023_09_21
# PATH=$PATH:/home/data/opt/biosoft/gatk-4.1.2.0/

# ---- 标注 ??? 存疑 / 路径可能不存在 ----
# PATH=$PATH:/opt/biosoft/Gblocks_0.91b/ ??? 空
# PATH=$PATH:/opt/biosoft/exonerate-2.2.0-x86_64/bin/  ??? 没有
# PATH=/home/data/opt/biosoft/MCScanX/downstream_analyses:$PATH ??? 没有
# PATH=/home/data/opt/biosoft/MCScanX/downstream_analyses/:$PATH  ??? 没有
# PATH=$PATH:/opt/biosoft/eggnog-mapper-1.0.3  ??? 没找到
# PATH=/home/data/opt/biosoft/GetOrganelle-master/:$PATH  # 可能移动有问题

# ---- 被新版本覆盖而失效的 Augustus 3.3.2 (CONFIG 已被 3.4.0 覆盖) ----
# export AUGUSTUS_CONFIG_PATH=/home/data/opt/biosoft/augustus-3.3.2/config/
#   注: 3.3.2 的 bin/scripts 仍在 PATH 中生效，已迁移到 bio/augustus.sh

# ---- 被删除/注释的 JAVA / make / sra 旧写法 ----
#export JAVA_HOME=/opt/biosoft/openjdk17/jdk-17.0.2  # 删了
#export PATH=$JAVA_HOME/bin:$PATH
# make4
# export PATH=/opt/biosoft/make-4.4/build:$PATH
# sra toolkit
# PATH=$PATH:/opt/biosoft/sratoolkit/bin/

# ---- 注释掉的 R-4.2.2 直接加载写法 (现由 sys/rlang.sh 的 r422 alias 提供) ----
# PATH=/home/data/opt/biosoft/R40/bin:$PATH
# export R_HOME=/home/data/opt/biosoft/R-422/lib64/R
# export R_LIBS=/home/data/opt/biosoft/R-422/lib64/R/library
# export LD_LIBRARY_PATH=$R_HOME/lib64/R/lib:$LD_LIBRARY_PATH
# PATH=$R_HOME/bin:$PATH

# ---- 其他注释 ----
# source /home/data/opt/biosoft/root/bin/thisroot.sh
#export PROTHINT_PATH=/home/data/opt/biosoft/ProtHint-2.6.0/bin
#export GENEMARK_PATH=/opt/biosoft/gmes_linux_64_4/
#export PATH=/home/train/software/ProtHint/bin/:$PATH
# source $(dirname "${BASH_SOURCE[0]}")/ngrok.sh
# 原 freetype 的错误路径写法 (liblib 笔误):
# export PKG_CONFIG_PATH=/home/data/opt/syssoft/freetype-2.14.1/build/liblib/pkgconfig:$PKG_CONFIG_PATH
