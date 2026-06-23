#!/bin/zsh
# 功能: tmux会话命令历史持久化脚本（zsh兼容版）
# 描述: 在tmux会话中自动记录命令历史，分别存储到全局tmux历史文件和当前工作目录的历史文件中
# 特性:
#   1. 排除ls/l/pwd/clear/exit/history等常用无需记录的命令
#   2. 每条记录包含执行用户、时间、工作目录和完整命令
#   3. 支持动态检测tmux环境，进入tmux后自动生效
#   4. 使用add-zsh-hook注册钩子，避免与其他插件冲突

# 确保zsh钩子模块加载
autoload -Uz add-zsh-hook

# 全局配置
HISTORY_EXCLUDE_COMMANDS=("ls" "l" "ll" "la" "pwd" "cd" ".." "." "clear" "exit" "history" "htop" "top" "df" "du" "free" "whoami" "date" "cal" "man" "info" "ps" "uptime" "w" "id" "groups" "uname" "hostname" "tty" "which" "whereis" "whatis" "alias" "jobs" "fg" "bg")
export HISTTIMEFORMAT="%F %T "

# 临时变量
LAST_CMD=""
LAST_CMD_TIME=""
LAST_CMD_CWD=""

# 目录切换钩子
_tmux_chpwd() {
  # 只有在tmux环境中才生效
  if [[ -n $TMUX ]]; then
    export CWD_HISTFILE="$PWD/.${USER}_history"
  fi
}

# 命令执行前钩子
_tmux_preexec() {
  # 只有在tmux环境中才生效
  if [[ -z $TMUX ]]; then
    LAST_CMD=""
    return
  fi
  
  [[ -z "$1" ]] && return
  local cmd="${1%% *}"
  for exclude in "${HISTORY_EXCLUDE_COMMANDS[@]}"; do
    if [[ "$cmd" == "$exclude" ]]; then
      LAST_CMD=""
      return
    fi
  done
  
  # 保存命令执行前信息
  LAST_CMD="$1"
  LAST_CMD_TIME=$(date "+%F %T")
  LAST_CMD_CWD="$PWD"
}

# 命令执行后钩子
_tmux_precmd() {
  # 只有在tmux环境中才生效
  if [[ -z $TMUX ]]; then
    LAST_CMD=""
    return
  fi
  
  local exit_code=$?
  # 没有需要记录的命令直接返回
  [[ -z "$LAST_CMD" ]] && return
  
  # 动态获取当前tmux会话信息
  local TMUX_SESSION_INFO=$(tmux display-message -p '#S-#I-#P' 2>/dev/null)
  local GLOBAL_HISTFILE="$HOME/.tmux_history/bash_history_$TMUX_SESSION_INFO.txt"
  local CWD_HISTFILE="${CWD_HISTFILE:-$PWD/.${USER}_history}"
  
  # 确保历史文件目录存在
  mkdir -p "$(dirname "$GLOBAL_HISTFILE")" 2>/dev/null
  
  # 判断命令执行状态
  local status_tag status_color
  if [[ $exit_code -eq 0 ]]; then
    status_tag="[Success]"
    status_color="\e[32m"
  else
    status_tag="[Error:${exit_code}]"
    status_color="\e[31m"
  fi

  # 写入历史文件
  echo -e "\e[35m${USER}\e[0m ${status_color}${status_tag}\e[0m \e[32m${LAST_CMD_TIME}\e[0m \e[34m${LAST_CMD_CWD}\e[0m: ${LAST_CMD}" >> "$GLOBAL_HISTFILE"
  echo -e "\e[35m${USER}\e[0m ${status_color}${status_tag}\e[0m \e[32m${LAST_CMD_TIME}\e[0m \e[34m${LAST_CMD_CWD}\e[0m: ${LAST_CMD}" >> "$CWD_HISTFILE"
  
  # 重置临时变量
  LAST_CMD=""
}

# 注册zsh钩子
add-zsh-hook chpwd _tmux_chpwd
add-zsh-hook preexec _tmux_preexec
add-zsh-hook precmd _tmux_precmd

# 初始化CWD_HISTFILE
_tmux_chpwd