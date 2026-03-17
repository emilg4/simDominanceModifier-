# Imports necessary modules
import json
import itertools

# Load parameter ranges from json file
# Insert directory of json file:

with open('*.json') as f:
    ranges = json.load(f)

# Prepares keys and values for combinations
keys = list(ranges.keys())
values = [ranges[k] for k in keys]

# Generates all combinations and write PARAM files
for i, combo in enumerate(itertools.product(*values)):
    params = dict(zip(keys, combo))
    params['NAME'] = f"sim1_{i}" 
    params['rep'] = i
    all_keys = list(params.keys())

# Insert directory of parameter files:

    param_filename = f'/gpfs/home/yyb21bru/scratch/simDominanceModifier-/parameters_v2/PARAM_{i}.csv'
    with open(param_filename, 'w') as csvf:
        csvf.write(','.join(all_keys) + '\n')
        csvf.write(','.join(str(params[k]) for k in all_keys) + '\n')
'EOF'
