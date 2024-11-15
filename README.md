# Panmama-Snakemake Workflow Documentation

## Quickstart

1. **Install Required Software and Packages**

   To run the Panmama-Snakemake workflow, the following software and packages are required:

   - **Snakemake**: Install into conda environment following instructions at : https://snakemake.readthedocs.io/
   - **Panmap**: Install following the instructions at: https://github.com/amkram/panmap
   - **Swampy**: Install from the relevant repository: https://github.com/your-swampy-repo

2. **Enable slurm (optional)**

   To run the Panmama-Snakemake workflow with SLURM, install SLURM plugin for Snakemake. 
   ```
   conda activate snakemake
   pip install snakemake-executor-plugin-slurm matplotlib
   ```
   SLURM job parameters are set in `workflow/profiles/slurm/config.yaml`.

3. **Set Up Config File**
   
   Set up the configuration file (`workflow/config.yaml`) and save in the desired directory for pipeline outputs.

4. **Run Pipeline**

   With slurm:

   ```bash
   snakemake -d <out_dir> --conda-prefix <conda_dir> --workflow-profile /workflow/profiles/slurm
   ```

   Without slurm:

   ```bash
   snakemake -d <out_dir> --conda-prefix <conda_dir> --sdm conda --cores 8
   ```

   - `<out_dir>` is the absoluate path to the output directory containing the `config.yaml` file with parameters for the run. Every pipeline run will have a different `<out_dir>`.
   - `<conda_dir>` is an absolute path to a location to build conda environments. All pipeline runs should use the same `<conda_dir>` to avoid rebuilding conda environments unnecessarily. 

---

## Config File Details

1. **Paths to External Tools**:
   - **`panmanUtils_path`**: Location of the `panmanUtils` executable (used for some panmap utilities).
   - **`panmap_path`**: Location of the `panmap` executable, which handles the main panmap processing.
   - **`swampy_path`**: Location of the `simulate_metagenome.py` script for the "swampy" simulator.

2. **Replicates**:
   - **`n_replicates`**: Defines the number of random sequence selections or replicates to generate for each condition

3. **Tree Selection**:
   - **`tree`**: Defines the panman for sequence selection. 
   - Use `"ecoli1000"`, `"sars20000"`, `"hiv20000"`, `"klebs1000"`, `"rsv4000"`, or `"tb400"` to access panmans published online.
   - Or, save a custom panman in `resources/panmans` and use the file prefix in config file. File prefixes cannot contain underscores.
   - One or multiple trees allowed.

4. **Abundances**:
   - **`abundances`**: Specifies the relative abundance of haplotypes, given as a list of fractions that sum to 1.
   - One or multiple lists allowed. 

5. **Simulator**:
   - **`simulator`**: Specifies which simulator to use for sequencing reads (`"iss"` (for InSilicoSeq) or `"swampy"`). 
   - One or multiple simulators allowed.  

6. **Reads Configuration**:
   - **`n_reads`**: The number of reads, or average depth (InSilicoSeq only), to simulate.
   - When using depth, set **`use_depth`** to `True`.
   - One or multiple values allowed.

7. **InSilicoSeq Parameters** (`iss_params`):
   - **`model`**: Specifies the model used by the `iss` simulator (`"miseq"`, `"hiseq"`, or `"novaseq"`).
   - **`use_depth`**: A boolean flag that, if set to `True`, will compute `n_reads` based on depth and genome size.
   - **`cpus`**: The number of CPU cores to use for InSilicoSeq simulations.

8. **Swampy Parameters** (`swampy_params`):
   - **`primer_set`**: Defines which set of primers to use for the swampy simulator (`"a1"`, `"a4"`, `"a5"`, or `"n2"`).

9. **Panmap Parameters** (`panmap_params`):
   - **`len_k_s_mers`**: Two integers separated by a space specifying the lengths for k-mers and s-mers (e.g., `10 5`). These represent the lengths of the sequencing fragments.
   - **`check_frequency`**: Integer value for frequency checking in panmap simulations.
   - **`redo_read_threshold`**: Threshold for redoing read iterations.
   - **`remove_iteration`**: Number of iterations after which to remove certain data.
   - **`cpus`**: The number of CPUs to be used by panmap.
   - One or multiple values are allowed for all panmap paramters. 

Example `config.yaml`:
```yaml
panmanUtils_path: "~/installs/panmap/build/bin/panmanUtils"
panmap_path: "~/installs/panmap/build/bin/panmap"
swampy_path: "~/installs/SWAMPy/src/simulate_metagenome.py"
n_replicates: 20
tree: ["rsv4000", "sars20000"]
abundances:
  - [1]
  - [.5, .25, .25]
simulator: "iss"
n_reads: 100
iss_params:
  model: "miseq"
  use_depth: True # Will calculate actual n_reads required to achieve average depth of 100x
  cpus: 8
swampy_params: 
  primer_set: "a5"
panmap_params:
  len_k_s_mers: [10 5, 17 9]
  check_frequency: 5
  redo_read_threshold: 0
  remove_iteration: 5
  cpus: 4

```

---
## `panmap` Wildcard Handling

This document explains how the Snakemake workflow dynamically handles wildcards for `panmap` parameter combinations and abundance lists, incorporating recent updates to the script. Add new parameters by updating the config file and the `run_panmap` rule in the snakefile. 

### 1. **Define and Expand `panmap_params` in Config**

The `panmap_params` section in the config file allows flexible parameter configurations. These parameters can be specified as single values, lists, or pairs (e.g., `len_k_s_mers` for `k` and `s` values). Here’s an example:

```yaml
panmap_params:
  len_k_s_mers: [10 5, 15 7]  # Pairs of k-mer and s-mer lengths
  check_frequency: [5, 10]
  redo_read_threshold: [0]
  remove_iteration: [5, 10]
  cpus: 4
```

- In this example:
  - `len_k_s_mers` defines two pairs of values: `(10, 5)` and `(15, 7)`.
  - `check_frequency` specifies two options: `5` and `10`.
  - All combinations of these parameters will be tested.

---

### 2. **Generate All Parameter Combinations with `generate_panmap_param_wildcards`**

The `generate_panmap_param_wildcards` function creates a unique wildcard for each combination of `panmap_params`:

```python
from itertools import product

def generate_panmap_param_wildcards(config):
    param_keys = list(config['panmap_params'].keys())
    param_values = [
        config['panmap_params'][key] if isinstance(config['panmap_params'][key], list) else [config['panmap_params'][key]]
        for key in param_keys
    ]
    param_combinations = list(product(*param_values))
    wildcard_map = {f"panmap{i+1}": dict(zip(param_keys, combination)) for i, combination in enumerate(param_combinations)}
    return wildcard_map
```

For the example config:
```yaml
panmap_params_wildcards = {
    "panmap1": {'len_k_s_mers': "10 5", 'check_frequency': 5, 'redo_read_threshold': 0, 'remove_iteration': 5, 'cpus': 4},
    "panmap2": {'len_k_s_mers': "10 5", 'check_frequency': 10, 'redo_read_threshold': 0, 'remove_iteration': 5, 'cpus': 4},
    ...
}
```

Each combination is assigned a unique name, like `panmap1`, `panmap2`, etc.

---

### 3. **Map Abundance Lists to Wildcards**

Haplotype abundance lists are mapped to wildcards using the `map_abundance_wildcards` function:

```python
def map_abundance_wildcards(abundances):
    length_count = {}
    abundance_mapping = {}

    for abundance in abundances:
        length = len(abundance)
        suffix = chr(97 + length_count.get(length, 0))  # 'a', 'b', etc.
        unique_name = f"{length}hap-{suffix}"
        length_count[length] = length_count.get(length, 0) + 1
        abundance_mapping[unique_name] = abundance

    return abundance_mapping
```

For abundances like `[[1], [0.5, 0.5], [0.3, 0.7]]`, the wildcards generated are:
```python
abundance_wildcards = {
    "1hap-a": [1],
    "2hap-a": [0.5, 0.5],
    "2hap-b": [0.3, 0.7]
}
```

---

### 4. **Rule Expansion Using Wildcards**

The `all` rule expands outputs based on `panmap_params_wildcards`, `abundance_wildcards`, and replicates (`reps`):

```python
rule all:
    input:
        expand("panmap_outputs/{i}/{tree_name}_{num_sequences}_{simulator}_{n_reads}_{panmap_params}_abundance_{i}.tsv",
               tree_name=config['tree'],
               num_sequences=abundance_wildcards.keys(),
               simulator=config['simulator'],
               n_reads=config['n_reads'],
               panmap_params=panmap_params_wildcards.keys(),
               i=reps),
        ...
```

This creates unique outputs for each combination of parameters and abundances.

---

### 5. **Log Wildcard Mappings**

To maintain a record of parameter mappings, the `log_wildcard_mapping` rule logs `abundance_wildcards` and `panmap_params_wildcards`:

```python
rule log_wildcard_mapping:
    output:
        "results/parameter.log"
    params:
        abundance_wildcards=abundance_wildcards,
        panmap_params_wildcards=panmap_params_wildcards
    shell:
        """
        echo "Haplotype Abundances:\n{params.abundance_wildcards}\n" >> {output}
        echo "Panmap Parameters:\n{params.panmap_params_wildcards}\n" >> {output}
        """
```

---

### 6. **Access Parameters in `run_panmap`**

The `run_panmap` rule accesses parameter combinations through wildcards:
```python
params:
    panmap_config = lambda wildcards: panmap_params_wildcards[wildcards.panmap_params]
```

This ensures the correct parameters are used in the shell command:
```bash
-k {params.panmap_config[len_k_s_mers].split()[0]} \
-s {params.panmap_config[len_k_s_mers].split()[1]} \
--check-frequency {params.panmap_config[check_frequency]} \
...
```

These updates streamline parameter exploration and ensure reliable file management for complex Snakemake workflows.