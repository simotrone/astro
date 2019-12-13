SEED=${1:-1}
OPTIONS=${@:2}
export PYTHONPATH=.

CRAB=( 83.6331  22.0145 )
model="only_background.xml"
caldb="prod3b-v2"
irf="South_z20_0.5h"
irf_file="irf_prod3b_v2_South_z20_0.5h.fits"
energy=( 0.025  150.0 )
time=( 0  1200 ) # total observation time
TIME_STEP=100
obs_time_limit=1200 # the single observation time limit
# pointing wrt source
# { 'name': 'R', 'ra': -0.5, 'dec':  0.0 }
# { 'name': 'U', 'ra':  0.0, 'dec': -0.5 }
# { 'name': 'L', 'ra': +0.5, 'dec':  0.0 }
# { 'name': 'B', 'ra':  0.0, 'dec': +0.5 }
names=(p0_R  p1_U  p2_L  p3_B)
ras=(  +0.5   0.0  -0.5   0.0)
decs=(  0.0  -0.5   0.0  +0.5)

HEADER='write_it'
OUTPUT_FILE="results/results_$SEED.tsv"
mkdir -p results
# check on output file
[ -e $OUTPUT_FILE ] && echo "cannot overwrite existing file '$OUTPUT_FILE'" && exit 1

running_time_start=${time[0]}
# loop over parameters
for i in ${!names[@]}; do
	# the simulation base
	pnt_ra=`python  -c "print(${CRAB[0]}+${ras[$i]})"`
	pnt_dec=`python -c "print(${CRAB[1]}+${decs[$i]})"`
	working_dir="${names[$i]}/$SEED"
	events_file="$working_dir/events.fits"
	# simulations if events file is missing
	if [ ! -e $events_file ]; then
		mkdir -p $working_dir
		final_simulation_time=$(( $running_time_start + $obs_time_limit ))
		echo "generating $events_file [$running_time_start - $final_simulation_time sec]..."
		ctobssim seed=$SEED ra=$pnt_ra dec=$pnt_dec rad=5.0 \
		   tmin=$running_time_start tmax=$final_simulation_time emin=${energy[0]} emax=${energy[1]} \
		   caldb=$caldb irf=$irf inmodel=$model outevents=$events_file logfile="$working_dir/ctobssim.log"
	fi

	while [ $running_time_start -lt $final_simulation_time ]; do
		running_time_stop=$(( $running_time_start + $TIME_STEP ))
		select_working_dir="$working_dir/sel_${running_time_stop}"
		select_events_file="$select_working_dir/events.fits"
		# selecting simulated events
		if [ ! -e $select_events_file ]; then
		  mkdir -p $select_working_dir
		  echo "selecting $select_events_file [$running_time_start - $running_time_stop sec]"
		  ctselect inobs=$events_file outobs=$select_events_file ra=INDEF dec=INDEF rad=INDEF \
		     tmin=$running_time_start tmax=$running_time_stop emin=INDEF emax=INDEF \
		     logfile="$select_working_dir/ctselect.log"
		fi

		[ ! -r $select_events_file ] && echo "cannot analyze '$select_events_file' because not readable" && exit 1
		results=`python extract_photometric_data.py --source-ra ${CRAB[0]} --source-dec ${CRAB[1]} --pointing-ra $pnt_ra --pointing-dec $pnt_dec -irf $irf_file -emin ${energy[0]} -emax ${energy[1]} --livetime $TIME_STEP $select_events_file`
		[ $? ] || exit 1
		if [ "x$HEADER" == "xwrite_it" ]; then
		  echo "ra	dec	seed	tmax	on	off	alpha	excess	li_ma	aeff	flux" >> $OUTPUT_FILE
		  HEADER='wrote'
		fi
		echo "$pnt_ra	$pnt_dec	$SEED	$running_time_stop	$results" >> $OUTPUT_FILE
		running_time_start=$running_time_stop
		[[ $running_time_stop -ge ${time[1]} ]] && echo "exit selecting loop" && break # the end of the selecting loop
	done
	[[ $running_time_stop -ge ${time[1]} ]] && echo "exit from main loop" && break # the end of the main loop
	# choose the time_step smartly :-)
	running_time_start=$final_simulation_time
done

for dir in ${names[@]}; do
  test $SEED -gt 5 && [ -d $dir/$SEED ] && rm -r "$dir/$SEED"
done
