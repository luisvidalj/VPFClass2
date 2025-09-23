# VPF-Class 2

**VPF-Class 2** is a machine learning–based tool for viral taxonomic classification. It takes raw viral sequences as input, predicts protein-coding genes, identifies marker protein families through HMMER searches, and uses a sparse neural network to assign taxonomic labels down to the genus level. The tool is designed to be modular, reproducible, and easy to use both for benchmarking against ICTV releases and for real-world applications where users submit novel viral contigs for classification.


## Setup and environment installation
Clone the repository and move into:

```bash
git clone https://github.com/luisvidalj/VPFClass2.git
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


Download the protein profiles and trained models (note: the archive is large, ~35 GB).
```bash
wget https://bioinfo.uib.es/~recerca/vpfclass2/tool_data_v0.4.0.tar.zst
```
Decompress the archive (requires `zstd` installed).
If you do not have `zstd` installed, you can install it as follows:
- **Linux (Debian/Ubuntu)**:
```bash
sudo apt update && sudo apt install zstd
```

- **masOS (Homebrew)**:
```bash
brew install zstd
```

- **Using conda (any platform)**
```bash
conda install -c conda-forge zstd
```

When `zstd` is available:
```bash
zstd -d --long=31 -T4 < tool_data_v0.4.0.tar.zst | tar -xvf -
pip install -e .
```

## Usage

VPF-Class 2 provides two main commands: `check` and `predict`.

### Check setup

Verify that the required models and marker profiles are correctly available:

```bash
vpfclass2 check --msl 40 --markers all
```
Where

- `--markers {all,virus}`:  
  - `virus`: only use profiles associated with viruses.  
  - `all`: use all profiles from GeNomad (viruses, plasmids, chromosomes).  

- `--msl {40,39,38,35,33}`: specify the ICTV MSL release for which the model was trained (this allows reproducing classifications as if using older releases). 

### Run predictions

```bash
vpfclass2 predict --fasta test/mini_test.fna --outdir test/output_test --markers virus --msl 39
``` 

Where the arguments:
- `--fasta`: input FASTA (.fasta, .fa or .fna) file containing viral contigs or genomes
- `--outdir`: output directory where results will be saved.
- `markers`: see above.
- `msl`: model version (ICTV release) to use. See above.

### Output structure

After running `vpfclass2 predict`, results will be saved in the specified `--outdir`.  
For example:

```bash
$ tree test/output_test/
test/output_test/
├── fasta
│   └── fasta_headers.csv
├── features
│   ├── accessions.txt
│   ├── features_counts_sparse.npz
│   ├── features_stats.json
│   └── vpf_to_index.used.json
├── hmmer
│   ├── hmm_hits.csv
│   ├── hmmsearch.stderr.log
│   ├── split_1.faa
│   ├── split_1.tbl
│   ├── split_2.faa
│   └── split_2.tbl
├── logs
├── preds
│   └── prediction.csv
├── prodigal
│   ├── output.faa
│   └── prodigal_proteins.csv
└── run_meta.json
