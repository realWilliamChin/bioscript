#!/bin/bash
# ============================================================
# 环境变量汇总入口
# 由原 general.sh 拆分而来，按 sys / bio / other 分类组织。
# 在 zshrc 中 source 本文件即可加载全部环境变量。
# 新增软件: 在对应分类目录下新建 .sh 文件即可被自动加载。
# ============================================================

# 本文件所在目录 (兼容 bash 与 zsh)
if [ -n "$ZSH_VERSION" ]; then
    ENV_ROOT="$(cd "$(dirname "${(%):-%x}")" && pwd)"
else
    ENV_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# 按 sys -> bio -> other 顺序加载各分类下的所有 .sh
for _category in sys bio other; do
    _dir="$ENV_ROOT/$_category"
    [ -d "$_dir" ] || continue
    for _f in "$_dir"/*.sh; do
        [ -f "$_f" ] && source "$_f"
    done
done
unset _category _dir _f

# 统一导出 PATH 等可能在子脚本中以非 export 形式赋值的变量
export PATH
export LD_LIBRARY_PATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH PKG_CONFIG_PATH CPATH LIBRARY_PATH PERL5LIB
