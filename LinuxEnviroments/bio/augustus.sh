# Augustus
# 注: AUGUSTUS_CONFIG_PATH 由 3.4.0 生效 (原 3.3.2 的同名设置被覆盖，已移至 general.sh 留存)。
#     PATH 中的可执行/脚本目录沿用原 general.sh 中实际生效的 3.3.2 路径。
export AUGUSTUS_CONFIG_PATH=/home/data/opt/biosoft/augustus-3.4.0/config
export AUGUSTUS_BIN_PATH=/home/data/opt/biosoft/augustus-3.4.0/bin
export AUGUSTUS_SCRIPTS_PATH=/home/data/opt/biosoft/augustus-3.4.0/scripts

PATH=$PATH:/home/data/opt/biosoft/augustus-3.3.2/bin/
PATH=$PATH:/home/data/opt/biosoft/augustus-3.3.2/scripts/
