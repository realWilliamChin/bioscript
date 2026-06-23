#================ qinqiang alias================
alias l='ls -lh'
alias les='less -Sm'
alias table='column -t'
alias csv="column -t -s ','"
function cpcd() {
  cp "$1" "$2" && cd "$2"
}
function cpbak() {
  cp "$1" "${1}.bak"
}
alias xlsx=in2csv

# Shell working directory reporting on zsh for Tabby
precmd() { echo -n "\x1b]1337;CurrentDir=$(pwd)\x07"; }
