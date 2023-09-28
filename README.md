# PHRIC

Prediction of complementary regions between RNA-RNA contacts derived from RIC-seq data. Developed by Sergei Margasyuk (smargasyuk@gmail.com) and Dmitri Pervouchine (pervouchine@gmail.com).

## Description

This package contains a pipeline for prediction of complementary regions between RNA-RNA contacts derived from RIC-seq data. (Cai et al., [2020](https://doi.org/10.1038/s41586-020-2249-1)). 

## Usage

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis.

    git clone https://github.com/smargasyuk/PHRIC.git
    cd PHRIC

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup (described in [Settings](config/README.md)).

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    # install mamba package manager if you don't have it
    conda install -n base -c conda-forge mamba
    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### Test run

To make a test run, type

```
snakemake --use-conda --configfile config/test_config.yaml -c8 test
```

The rule will execute the pipeline for toy dataset (RIC-seq data for `chr21`) and the output files in `results/test_hg19/S16.tsv` will be compared to reference version. The comparison is a Snakemake rule, so the correct completion of the test is indicated by a Snakemake success message (`.... of ... steps (100%) completed' in green text), and a test failure is indicated by a Snakemake error message in red text.


#### Run on full RIC-seq data

PHRIC pipeline starts from the compressed intermediate output of [RNAcontacts](https://github.com/smargasyuk/RNAcontacts.git) pipeline. Run the RNAcontacts pipeline on some RIC-seq data (we have used [GSE190214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190214)) and copy or softlink its output, `resources` folder, to `resources/RNAcontacts-output` of this pipeline. Then put the same sample sheet used in RNAcontacts (described in [Settings](config/README.md)) to `config/samples.tsv` and run PHRIC pipeline.

### Step 5: Investigate results

The output of the pipeline is the file `results/{genome}/S16.tsv` that contains the list of pairs of complementary regions between RIC-seq contacts. It contains the following fields:
1. first 12 field: record in BED12 describing pair of complementary regions;
2. `'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'`: coordinates of the left and right complementary regions;
3. `cell_lines_1`: list of experiments supporting the inner contact of the pair;
4. `n_cell_lines_1`: number of experiments supporting the inner contact;
5. `n_junctions_1`: number of different pairs of contacting points supporting the inner contact;
6. `n_reads_1`: number of reads supporting the inner contact;
7. fields similar to 3-6 with `_2` suffix correspond to the outer contact instead;
8. `strand_`: strand of the transcript containing the pair of contacts;
9. `gene_name`: name of the gene containing the pair of contacts;
10. `pair_id`: unique identifier of the pair of contacts; 
11. `structure`: dot-bracket string describing the structure between the complementary regions; 
12. `handle1_coord`: coordinates of complementary region inside the left folding region;
13. `energy`: free energy of pairing, kcal/mol;
14. `fold_interval_1_seq`: full sequence of the left folding region;
15. `handle1_seq`: sequence of left complementary region;
16. fields similar to 12,14,15 with `_2` suffix correspond to the right folding region instead;