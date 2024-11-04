
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/davneetkaur/miniconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/davneetkaur/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/davneetkaur/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/davneetkaur/miniconda3/bin:$PATH"
    fiB
fi
unset __conda_setup
# <<< conda initialize <<<
export PATH="/usr/local/bin:/opt/homebrew/bin:$PATH"
export PATH=$PATH:/Applications/FastQC.app/Contents/MacOS
