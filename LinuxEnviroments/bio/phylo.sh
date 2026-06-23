# 系统发育 / 进化分析 (prottest, RAxML, r8s, beast, tracer, CAFE, beagle-lib)
export PROTTEST_HOME=/home/data/opt/biosoft/prottest-3.4.2
PATH=$PATH:/home/data/opt/biosoft/RAxML-8.2.12/
PATH=$PATH:/home/data/opt/biosoft/r8s1.81/src
PATH=$PATH:/home/data/opt/biosoft/beast/bin/
PATH=$PATH:/home/data/opt/biosoft/Tracer_v1.7.1/bin/
PATH=$PATH:/home/data/opt/biosoft/CAFE-4.2.1/release/

# beagle-lib (beast 依赖)
export PKG_CONFIG_PATH=/home/data/opt/biosoft/beagle-lib-3.1.2/lib/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/home/data/opt/biosoft/beagle-lib-3.1.2/lib/:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/home/data/opt/biosoft/beagle-lib-3.1.2/include:$C_INCLUDE_PATH
