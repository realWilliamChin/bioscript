# 系统库 (R 包编译依赖): freetype, libwebp, curl, openssl, hdf5, udunits2
# 注: sqlite 在 bio/gis_stack.sh 中与 GDAL/PROJ 一起配置

# freetype (安装 R 包时提示版本太低)
export PKG_CONFIG_PATH=/home/data/opt/syssoft/freetype-2.14.1/build/lib/pkgconfig:$PKG_CONFIG_PATH
export CPATH=/home/data/opt/syssoft/freetype-2.14.1/build/include/freetype2:/home/data/opt/syssoft/freetype-2.14.1/build/include:$CPATH
export LIBRARY_PATH=/home/data/opt/syssoft/freetype-2.14.1/build/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/data/opt/syssoft/freetype-2.14.1/build/lib:$LD_LIBRARY_PATH

# libwebp (安装 tidyverse 时提示版本太低)
export PKG_CONFIG_PATH=/home/data/opt/syssoft/libwebp-1.3.2/build/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/home/data/opt/syssoft/libwebp-1.3.2/build/lib:$LD_LIBRARY_PATH
export LDFLAGS="-L/home/data/opt/syssoft/libwebp-1.3.2/build/lib"

# curl
export PATH=/home/data/opt/syssoft/curl-8.6.0/build/bin:$PATH
export LD_LIBRARY_PATH=/home/data/opt/syssoft/curl-8.6.0/build/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/home/data/opt/syssoft/curl-8.6.0/build/lib/pkgconfig:$PKG_CONFIG_PATH
export CPATH=/home/data/opt/syssoft/curl-8.6.0/build/include:$CPATH

# OpenSSL 1.1.1w
export OPENSSL_HOME=/home/data/opt/syssoft/openssl-1.1.1w/build
export PATH=$OPENSSL_HOME/bin:$PATH
export LD_LIBRARY_PATH=$OPENSSL_HOME/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$OPENSSL_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$OPENSSL_HOME/include:$CPLUS_INCLUDE_PATH
export PKG_CONFIG_PATH=$OPENSSL_HOME/lib/pkgconfig:$PKG_CONFIG_PATH

# HDF5 (安装 loupeR 包时提示版本太低)
export HDF5_HOME=/home/data/opt/syssoft/hdf5-1.14.6/build
export PATH=$HDF5_HOME/bin:$PATH
export LD_LIBRARY_PATH=$HDF5_HOME/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$HDF5_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$HDF5_HOME/include:$CPLUS_INCLUDE_PATH
export PKG_CONFIG_PATH=$HDF5_HOME/lib/pkgconfig:$PKG_CONFIG_PATH
export HDF5_DIR=$HDF5_HOME

# UDUNITS2
export UDUNITS2_HOME=/home/data/opt/syssoft/udunits2/build
export PATH=$UDUNITS2_HOME/bin:$PATH
export LD_LIBRARY_PATH=$UDUNITS2_HOME/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$UDUNITS2_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$UDUNITS2_HOME/include:$CPLUS_INCLUDE_PATH
export UDUNITS2_INCLUDE=$UDUNITS2_HOME/include
export UDUNITS2_LIBS=$UDUNITS2_HOME/lib
