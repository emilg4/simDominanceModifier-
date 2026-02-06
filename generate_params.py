# Imports necessary modules
import json
import itertools

# Load parameter ranges from json file
with open('/gpfs/home/yyb21bru/Project/param_ranges_test.json') as f:
    ranges = json.load(f)

# Prepares keys and values for combinations
keys = list(ranges.keys())
values = [ranges[k] for k in keys]

# Generates all combinations and write PARAM files
for i, combo in enumerate(itertools.product(*values)):
    params = dict(zip(keys, combo))
    param_filename = f'/gpfs/home/yyb21bru/scratch/simDominanceModifier-/PARAM_{i}.csv'
    with open(param_filename, 'w') as csvf:
        csvf.write(','.join(keys) + '\n')
        csvf.write(','.join(str(params[k]) for k in keys) + '\n')
'EOF'
