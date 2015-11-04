export ANALYSISFW_PATH=/afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/NA62Analysis
export ANALYSISFW_USERDIR=/afs/cern.ch/work/r/rmarchev/private/New_NA62Fw/analysis_new
export NA62MCSOURCE=/afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/NA62MC

if [[ ! $LD_LIBRARY_PATH =~ $ANALYSISFW_PATH/lib ]]; then
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANALYSISFW_PATH/lib
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANALYSISFW_USERDIR/lib
	export PATH+=:$ANALYSISFW_PATH
fi

source $ANALYSISFW_PATH/scripts/env.sh
