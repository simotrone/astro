import gammalib
import os
import sys

input_model = sys.argv[1] 
models = gammalib.GModels(input_model)
print(models)
print("TS:", models[0].ts())

