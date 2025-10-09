# VPF-Class 2

**VPF-Class 2** VPF-Class2 is a hybrid tool combining reference alignment and machine learning for viral taxonomic classification. It takes raw viral sequences as input, predicts protein-coding genes, identifies marker protein families through HMMER searches, and uses a neural network to assign taxonomic labels down to the genus level. The tool is designed to be modular, reproducible, and easy to use both for benchmarking against ICTV releases and for real-world applications where users submit novel viral contigs for classification.


## 1 - Setup and environment installation
Clone the repository and move into:

```bash
git clone https://github.com/luisvidalj/VPFClass2.git
git clone https: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//@github.com/luisvidalj/VPFClass2.git
cd VPFClass2
```

Dependencies installation is guided via Conda, which you will need to have installed. Select either option depending on whether your device has a GPU.

```bash
# For CPU-only systems
conda env create -f environment_cpu.yml
conda activate vpfclass

# For GPU-enabled systems
conda env create -f environment_gpu.yml
conda activate vpfclass-gpu
``` 
### Marker profiles selection

Before downloading the data, choose which marker profiles you want to use:

- **Virus markers** → smaller dataset (~17 GB), only includes profiles associated with viruses (**recommended**, see [paper reference] for details).  
- **Complete markers** → larger dataset (~30 GB), includes all GeNomad profiles (viruses, plasmids, chromosomes).


```bash
# Viral markers setup
bash scripts/setup_tool_data.sh virus

# Complete markers setup
bash scripts/setup_tool_data.sh all
```

> **Note:** If one dataset has already been installed, the other can also be added later — however, note that the total storage required will be approximately 60 GB.



## Usage

VPF-Class 2 provides two main commands: `check` and `predict`.

### Check setup

Verify that the required models and marker profiles are correctly installed.  
With the default setup, you should be able to successfully run:

```bash
vpfclass2 check --msl 40 --markers virus
```
Where

- `--markers {all,virus}`:  
  - `virus`: only use profiles associated with viruses.  
  - `all`: use all profiles from GeNomad (viruses, plasmids, chromosomes).  

- `--msl {40,39,38,35,33}`: specify the ICTV MSL release for which the model was trained (this allows reproducing classifications as if using older releases). 

### Run predictions

```bash
vpfclass2 predict --fasta test/mini_test.fna --outdir toy_example/readme --markers virus --msl 40
``` 

Where the arguments:
- `--fasta`: input FASTA (.fasta, .fa or .fna) file containing viral contigs or genomes
- `--outdir`: output directory where results will be saved.
- `markers`: see above.
- `msl`: model version (ICTV release) to use. See above.

Parallelization is automatically configured to use 8 CPUs by default.  
This can be adjusted with the `--num-cpus` flag (e.g., `--num-cpus 28`).

### Output structure

After running `vpfclass2 predict`, results will be saved in the specified `--outdir`.  
For example:

```bash
toy_example/readme/
├── 0_preds
│   └── prediction.csv
├── 1_prodigal
│   ├── logs
│   │   ├── shard_0001.stderr.log
│   │   └── shard_0001.stdout.log
│   ├── output.faa
│   ├── shard_outputs
│   └── shards
├── 2_hmmer
│   ├── hmmsearch.stderr.log
│   ├── split_1.faa
│   ├── split_1.tbl
│   ├── split_2.faa
│   └── split_2.tbl
└── run_meta.json
```