python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800 --dec-shift 0.5 --dir dec_0.5 --model source_model.xml --save &
python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800 --dec-shift 1.0 --dir dec_1.0 --model source_model.xml --save &
# ra +0.5 doesn't provide an on/off analysis
#python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800 --ra-shift 0.5  --dir ra_0.5  --model source_model.xml --save &
#python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800 --ra-shift 0.6  --dir ra_0.6  --model source_model.xml --save &
python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800 --ra-shift 1.0 --dec-shift 1.0 --dir ra_1.0_dec_1.0 --model source_model.xml --save &
python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800 --ra-shift 2.0 --dec-shift 2.0 --dir ra_2.0_dec_2.0 --model source_model.xml --save &
