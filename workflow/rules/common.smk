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
    return [f'results/{wildcards["genome"]}/S0/{r.project}/{r.sample_name}/{jtype}.tsv.gz' for r in relevant_samples.itertuples() for jtype in ["Neo", "Chimeric"]]


def get_chroms(wildcards):
    with open(f'resources/genomes/{wildcards.genome}/chromSizes') as f:
        ls = [l.strip() for l in f.readlines() if l.strip()]
        return [l.split('\t')[0] for l in ls]


def get_all_genomes(wildcards):
    return samples.genome.unique().tolist()


def get_all_outputs(wildcards):
    return expand('results/{genome}/S16.bed', genome=get_all_genomes(wildcards))


def postprocess_preph_table(t0):
    t1 = t0.drop(columns=31)
    t1 = t1.iloc[:, 12:]
    t1.columns = S16_other_columns
    t1 = pd.concat([t1,
                    t1.handle1_coord.str.split(',', expand=True).astype(int),
                    t1.handle2_coord.str.split(',', expand=True).astype(int).rename(columns={0:2, 1:3})], axis=1)
    t1['handle1_seq'] = t1.apply(lambda r: r['fold_interval_1_seq'][r[0]:r[1]+1], axis=1)
    t1['handle2_seq'] = t1.apply(lambda r: r['fold_interval_2_seq'][r[2]:r[3]+1], axis=1)
    t2 = t1.copy()
    t2['start1'] = t1.apply(lambda r: r.start1 + r[0] if r.strand_ == '+' else r.end1 - r[1], axis=1)
    t2['end1'] = t1.apply(lambda r: r.start1 + r[1] if r.strand_ == '+' else r.end1 - r[0], axis=1) + 1
    t2['start2'] = t1.apply(lambda r: r.start2 + r[2] if r.strand_ == '+' else r.end2 - r[3], axis=1)
    t2['end2'] = t1.apply(lambda r: r.start2 + r[3] if r.strand_ == '+' else r.end2 - r[2], axis=1) + 1
    t2 = t2.drop(columns=[0,1,2,3])

    t2_bed = pd.DataFrame.from_dict({
        'chrom': t2.chrom1,
        'chromStart': t2.start1,
        'chromEnd': t2.end2,
        'name': t2.pair_id,
        'score':  (- t2.energy * 10).astype(int),
        'strand': t2.strand_,
        'thickStart': t2.start1,
        'thickEnd': t2.end2,
        'itemRgb': '0,0,0',
        'blockCount': 2,
        'blockSizes': (t2.end1 - t2.start1).astype(str) + ',' + (t2.end2 - t2.start2).astype(str),
        'blockStarts':'0,' + (t2.start2 - t2.start1).astype(str)
    })

    return pd.concat([t2_bed, t2], axis=1)
