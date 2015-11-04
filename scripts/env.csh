setenv ANALYSISFW_PATH /afs/cern.ch/user/r/rmarchev/work/NA62Analysis
setenv ANALYSISFW_USERDIR /afs/cern.ch/work/r/rmarchev/private/New_NA62Fw/analysis_new
setenv NA62MCSOURCE /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/NA62MC

if ( "$LD_LIBRARY_PATH" !~ "*$ANALYSISFW_PATH/lib*" ) then
	setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ANALYSISFW_PATH}/lib
	setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ANALYSISFW_USERDIR}/lib
	setenv PATH ${PATH}:${ANALYSISFW_PATH}
endif

source $ANALYSISFW_PATH/scripts/env.csh
