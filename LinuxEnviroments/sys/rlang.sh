# R 多版本切换 alias (默认不激活，按需执行)
alias r422='export R_HOME=/home/data/opt/biosoft/R-422/lib64/R; \
export R_LIBS=/home/data/opt/biosoft/R-422/lib64/R/library; \
export LD_LIBRARY_PATH=$R_HOME/lib64/R/lib:$LD_LIBRARY_PATH; \
export PATH=$R_HOME/bin:$PATH'

alias r443='export R_HOME=/home/data/opt/biosoft/R-4.4.3/build/lib64/R; \
export R_LIBS=/home/data/opt/biosoft/R-4.4.3/build/lib64/R/library; \
export LD_LIBRARY_PATH=$R_HOME/lib:$LD_LIBRARY_PATH; \
export PATH=$R_HOME/bin:$PATH'
