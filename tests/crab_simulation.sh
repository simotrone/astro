PYTHONPATH=../../
python events_generation.py --model ../../crab_simulations/crab.xml --save --dir crab_sim --seed 3 --tmax 1800 --name Crab --ra 83.6331 --dec 22.0145
python plot_spectrum.py crab_sim/3/test_events.fits --model ../../crab_simulations/crab.xml --name Crab --save --dir crab_sim/3
