from pathlib import Path
import pandas as pd


def foldint_subseq(s, ind):
    ind_int = [int(i) for i in ind.split(',')]
    return s[ind_int[0]:ind_int[1]+1]


S16_bed12_columns = 'chrom	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb	blockCount	blockSizes	blockStarts'.split('\t')
S16_other_columns = ['chrom1', 'start1', 'end1','chrom2', 'start2', 'end2',
                 'cell_lines_1', 'n_cell_lines_1', 'n_junctions_1', 'n_reads_1',
                 'cell_lines_2', 'n_cell_lines_2', 'n_junctions_2', 'n_reads_2',
                 'strand_', 'gene_name', 'pair_id', 'structure', 'handle1_coord',
                 'handle2_coord', 'energy', 'fold_interval_1_seq', 'fold_interval_2_seq']
S16_columns = S16_bed12_columns + S16_other_columns


samples = (
    pd.read_csv(config["samples"], sep="\t")
    .applymap(lambda x: x.strip() if isinstance(x, str) else x)
    .set_index("sample_name", drop=False)
)


def get_junctions(wildcards):
    relevant_samples = samples.loc[(samples.genome == wildcards["genome"]) & (samples.treatment == "experiment")]
    return [f"resources/RNAcontacts-output/{wildcards['genome']}/{r.project}/junctions/{r.sample_name}/{jtype}.tsv.gz"
        for r in relevant_samples.itertuples()
        for jtype in ["Neo", "Chimeric"]]


def get_chroms(wildcards):
    with open(f'resources/genomes/{wildcards.genome}/chromSizes') as f:
        ls = [l.strip() for l in f.readlines() if l.strip()]
        return [l.split('\t')[0] for l in ls]


def get_all_genomes(wildcards):
    return samples.genome.unique().tolist()


def get_all_outputs(wildcards):
    return expand('results/{genome}/S16.bed', genome=get_all_genomes(wildcards))