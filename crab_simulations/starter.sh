SEED=${1:-1}
OPTIONS=${@:2}
MODEL="crab.xml"

python run_simulations.py $MODEL --dir dec_0.5 --tmax 1800 --dec-shift 0.5 --model $MODEL --seed $SEED $OPTIONS &
python run_simulations.py $MODEL --dir dec_1.0 --tmax 1800 --dec-shift 1.0 --model $MODEL --seed $SEED $OPTIONS &
python run_simulations.py $MODEL --dir ra_0.5 --tmax 1800 --ra-shift 0.5 --model $MODEL --seed $SEED $OPTIONS &
#python run_simulations.py $TS --dir ra_0.6 --tmax 1800 --ra-shift 0.6 --model source_model.xml --seed $SEED &
python run_simulations.py $MODEL --dir ra_1.0_dec_1.0 --tmax 1800 --ra-shift 1.0 --dec-shift 1.0 --model $MODEL --seed $SEED $OPTIONS &
python run_simulations.py $MODEL --dir ra_2.0_dec_2.0 --tmax 1800 --ra-shift 2.0 --dec-shift 2.0 --model $MODEL --seed $SEED $OPTIONS 
