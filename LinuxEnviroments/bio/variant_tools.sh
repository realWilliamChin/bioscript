# 结构变异 / motif 等 (Manta, ROOT, CNVnator, MEME, SURVIVOR)

# Manta (结构变异检测)
export MANTA_HOME=/home/data/opt/biosoft/manta-1.6.0
export PATH=$MANTA_HOME/bin:$PATH

# ROOT (CNVnator 依赖的数据分析框架)
export ROOTSYS=/home/data/opt/biosoft/root
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

# CNVnator (拷贝数变异)
export CNVNATOR_HOME=/home/data/opt/biosoft/CNVnator_v0.4.1/src
export PATH=$CNVNATOR_HOME:$PATH

# MEME (motif 分析)
export MEME_BIN=/home/data/opt/biosoft/meme-5.5.9/build/bin
export MEME_HOME=/home/data/opt/biosoft/meme-5.5.9/build
export PATH=$MEME_HOME/bin:$MEME_HOME/scripts:$PATH

# SURVIVOR (SV 合并)
export PATH=/home/data/opt/biosoft/SURVIVOR/Debug/:$PATH
