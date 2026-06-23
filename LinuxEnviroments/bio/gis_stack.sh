# GIS 栈 (SQLite, PROJ, GDAL, GEOS) —— R 空间分析包 (sf/terra/rgdal) 编译运行依赖

# 根目录
export SQLITE_HOME=/home/data/opt/syssoft/sqlite-autoconf-3450100/build
export PROJ_HOME=/home/data/opt/biosoft/proj-9.3.1/build
export GDAL_HOME=/home/data/opt/biosoft/gdal-3.8.4/build

# 命令行工具 (gdal-config, proj, sqlite3)
export PATH=$GDAL_HOME/bin:$PROJ_HOME/bin:$SQLITE_HOME/bin:$PATH

# 运行时动态库
export LD_LIBRARY_PATH=$GDAL_HOME/lib:$PROJ_HOME/lib:$SQLITE_HOME/lib:$LD_LIBRARY_PATH

# 头文件 (R 包编译查找)
export C_INCLUDE_PATH=$GDAL_HOME/include:$PROJ_HOME/include:$SQLITE_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$GDAL_HOME/include:$PROJ_HOME/include:$SQLITE_HOME/include:$CPLUS_INCLUDE_PATH

# pkg-config
export PKG_CONFIG_PATH=$GDAL_HOME/lib/pkgconfig:$PROJ_HOME/lib/pkgconfig:$SQLITE_HOME/lib/pkgconfig:$PKG_CONFIG_PATH

# GIS 资源文件
export GDAL_DATA=$GDAL_HOME/share/gdal
export PROJ_LIB=$PROJ_HOME/share/proj

# lib64 / geos 补充库路径 (部分库装在 lib64)
export LD_LIBRARY_PATH=/home/data/opt/biosoft/gdal-3.8.4/build/lib64:/home/data/opt/biosoft/gdal-3.8.4/build/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/data/opt/biosoft/proj-9.3.1/build/lib64:/home/data/opt/biosoft/proj-9.3.1/build/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/data/opt/biosoft/geos-3.13.0/build/lib64:/home/data/opt/biosoft/geos-3.13.0/_build/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/data/opt/syssoft/libwebp-1.3.2/build/lib/:$LD_LIBRARY_PATH
