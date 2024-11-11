Run with 
```snakemake -d <wdir> --conda-prefix <conda> --sdm conda```
using an absolute path to a working directory for output files and an absolute path to a location in which to build the conda environments.

Pipeline can run on any panman available on web or saved as \<tree>.panman in resources/panmans. Specify \<tree> in config file. 

Here's how the `panmap` wildcards handling works step-by-step:

### 1. **Define and Expand `panmap_params` Combinations in Config**
   The `panmap_params` section in the config file provides a mix of lists and single integer values for each `panmap` parameter:
   ```yaml
   panmap_params:
     len_k_mer: [10, 15]
     len_s_mer: 5
     check_frequency: 5
     redo_read_threshold: 0
     remove_iteration: 5
     cpus: 4
   ```
   - `len_k_mer` is a list (`[10, 15]`), which means you want to test two different values for this parameter.
   - All other parameters are single values.

### 2. **Generate All Parameter Combinations with `generate_panmap_param_wildcards`**
   The function `generate_panmap_param_wildcards` reads the config parameters and generates all combinations of parameters to map each to a unique wildcard identifier:

   ```python
   from itertools import product

   def generate_panmap_param_wildcards(config):
       param_keys = list(config['panmap_params'].keys())
       param_values = [config['panmap_params'][key] if isinstance(config['panmap_params'][key], list) else [config['panmap_params'][key]] for key in param_keys]
       param_combinations = list(product(*param_values))
       wildcard_map = {f"panmap{i+1}": dict(zip(param_keys, combination)) for i, combination in enumerate(param_combinations)}
       return wildcard_map

   panmap_params_wildcards = generate_panmap_param_wildcards(config)
   ```
   Here’s how it works:
   - `param_keys`: List of parameter names, e.g., `['len_k_mer', 'len_s_mer', ...]`.
   - `param_values`: Lists of parameter values, ensuring that single values are wrapped in lists for compatibility.
   - `param_combinations`: Generates all possible combinations across the parameter lists. Here, two combinations are generated: one for `len_k_mer = 10` and another for `len_k_mer = 15` (since other parameters are single values).
   - `wildcard_map`: Each combination is assigned a unique name, like `panmap1` or `panmap2`.

   For the given config:
   ```python
   panmap_params_wildcards = {
       "panmap1": {'len_k_mer': 10, 'len_s_mer': 5, 'check_frequency': 5, 'redo_read_threshold': 0, 'remove_iteration': 5, 'cpus': 4},
       "panmap2": {'len_k_mer': 15, 'len_s_mer': 5, 'check_frequency': 5, 'redo_read_threshold': 0, 'remove_iteration': 5, 'cpus': 4}
   }
   ```

### 3. **Use Wildcards in the `all` Rule for File Expansion**
   In the `all` rule, `panmap_params_wildcards.keys()` (i.e., `["panmap1", "panmap2"]`) is used to generate multiple outputs:
   ```python
   expand("panmap_outputs/{tree_name}_{num_sequences}_{simulator}_{n_reads}_{panmap_params}_abundance_{i}.tsv",
          tree_name=config['tree'], 
          num_sequences=abundance_wildcards.keys(),
          simulator=config['simulator'],
          n_reads=config['n_reads'],
          panmap_params=panmap_params_wildcards.keys(),
          i=reps)
   ```
   This creates separate output paths for each `panmap_params` combination (e.g., `panmap1` and `panmap2`).

### 4. **Access Parameter Values in the `run_panmap` Rule Using Wildcards**
   In the `run_panmap` rule, `panmap_params_wildcards[wildcards.panmap_params]` is used to retrieve the specific parameter combination:
   ```python
   params:
       panmap_config = lambda wildcards: panmap_params_wildcards[wildcards.panmap_params]
   ```
   The `params` dictionary provides access to the specific parameter values for each wildcard combination, allowing those values to be used in the shell command.

### 5. **Using Parameters in the Shell Command**
   Each `panmap_params` combination is passed to `panmap` as unique arguments, enabling tailored runs:
   ```bash
   -k {params.panmap_config[len_k_mer]} \
   -s {params.panmap_config[len_s_mer]} \
   --check-frequency {params.panmap_config[check_frequency]} \
   --redo-read-threshold {params.panmap_config[redo_read_threshold]} \
   --remove-iteration {params.panmap_config[remove_iteration]} \
   --cpus {params.panmap_config[cpus]}
   ```
   This section ensures the command uses the correct parameters for each `panmap` run based on the specified wildcard combination, e.g., `panmap1` or `panmap2`.

### Summary
This setup enables you to:
   - Define multiple parameter combinations in the config file.
   - Automatically generate unique wildcards for each combination.
   - Pass specific parameter values to `panmap` in each run based on the selected combination. 

Each combination is assigned a unique name (`panmap1`, `panmap2`, etc.) and incorporated into filenames to distinguish outputs. This makes it easy to trace outputs to specific parameter combinations.