# CTA-containers/module/anaconda-3.6
# HOME_ANACONDA=/opt/anaconda3
HOME_ANACONDA=$HOME/anaconda3

export PATH=$HOME_ANACONDA/bin:$PATH
unset PYTHONPATH
export LD_LIBRARY_PATH=$HOME_ANACONDA/lib:$LD_LIBRARY_PATH
source $HOME_ANACONDA/etc/profile.d/conda.sh
