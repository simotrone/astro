CTOOLS_HOME=${1:-/opt}
# ctools activation:
export GAMMALIB=$CTOOLS_HOME
export CTOOLS=$CTOOLS_HOME
source $GAMMALIB/bin/gammalib-init.sh
source $CTOOLS/bin/ctools-init.sh

