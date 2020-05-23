PYTHONPATH=../../
python ../../grb_simulations/export_run.py ../../grb_simulations/run0406_ID000126.fits ../../grb_simulations/run0406_ID000126.xml --tmax 1800 --dir grb_sim --save grb_sim/timeslices.tsv
python events_generation_from_slices.py grb_sim/timeslices.tsv --dir grb_sim --model ../../grb_simulations/source_model.xml --save
python plot_spectrum.py grb_sim/1/sim_obs_list.xml --model ../../grb_simulations/source_model.xml --name run0406_ID000126 --save --dir grb_sim/1/
