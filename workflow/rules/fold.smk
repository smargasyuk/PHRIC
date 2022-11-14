rule filter_pairs_get_intercontacts:
    input:  'results/{genome}/S9.tsv'
    output: 'results/{genome}/S10.tsv'
    params: 
        min_support = config['min_contact_support']
    conda: "../envs/common.yaml"
    shell: """
cat {input} |\
awk -v 'OFS=\t' '($10>{params.min_support})&&($22>{params.min_support})' |\
awk -v 'OFS=\t' '{{print $1,$15,$2,$13,$6,$17,$7,$8,$9,$10,$19,$20,$21,$22,$11,$12,$25}}' > {output}
"""


rule filter_intercontacts:
    input: 
        pairs = 'results/{genome}/S10.tsv',
        rmsk = 'resources/genomes/{genome}/RMSK.bed'
    output: 'results/{genome}/S11.tsv'
    params:
        min_radius = config['distances']['nested_pass2']['rmin'],
        max_radius = config['distances']['nested_pass2']['rmax'],
    conda: "../envs/common.yaml"
    shell: """
paste <(cut -f1-6 {input.pairs}) {input.pairs} |\
awk -v 'OFS=\t' '($3-$2<={params.max_radius}) && ($3-$2>={params.min_radius})' |\
sort-bed - |\
bedops -n 1 - <(sort-bed {input.rmsk}) |\
cut -f4- |\
awk -v 'OFS=\t' '($3-$2<={params.max_radius}) && ($3-$2>={params.min_radius})' |\
sort-bed - |\
bedops -n 1 - <(sort-bed {input.rmsk}) |\
cut -f 4- > {output}
"""

rule get_intercontacts_fasta:
    input: 
        pairs = 'results/{genome}/S11.tsv',
        genome_fasta = 'resources/genomes/{genome}/genome.fa'
    output:
        seq_L = 'results/{genome}/S12L.fa',
        seq_R = 'results/{genome}/S12R.fa',
    conda: "../envs/common.yaml"
    shell: """
cat {input.pairs} |\
awk -v 'OFS=\t' '{{print $1,$2,$3,$17,".",$15}}' |\
bedtools getfasta -s -bed stdin -fi {input.genome_fasta} > {output.seq_L}

cat {input.pairs} |\
awk -v 'OFS=\t' '{{print $4,$5,$6,$17,".",$15}}' |\
bedtools getfasta -s -bed stdin -fi {input.genome_fasta} > {output.seq_R}
"""

rule PrePH:
    input:
        seq_L = 'results/{genome}/S12L.fa',
        seq_R = 'results/{genome}/S12R.fa',
    output: 'results/{genome}/S14_PrePH.tsv'
    conda: "../envs/PrePH.yaml"
    threads: 16
    shell: """
paste -d '\n' <(grep -v '>' {input.seq_L}) <(grep -v '>' {input.seq_R}) |\
grep -v "^$" |  xargs -n2 |\
parallel 2>/dev/null --keep-order -j{threads} --colsep ' ' "python workflow/scripts/PrePH/src/fold_SM.py -f {{1}} -s {{2}} \
-k 3 -a 3 -e -1 -u False -d 2" > {output}   
"""

rule merge_w_preph_results:
    input: 
        pairs = 'results/{genome}/S11.tsv',
        preph = 'results/{genome}/S14_PrePH.tsv',
        seq_L = 'results/{genome}/S12L.fa',
        seq_R = 'results/{genome}/S12R.fa',
    output: 'results/{genome}/S15.bed'
    conda: "../envs/common.yaml"
    shell: """
paste {input.pairs} {input.preph} <(grep -v '>' {input.seq_L}) <(grep -v '>' {input.seq_R}) |\
awk -v 'OFS=\t' '{{gsub(/[\(\)]/, "", $22)}}1' |\
awk -v 'OFS=\t' 'NF==24' |\
awk -v 'OFS=\t' '{{print $1,$2,$6,$17,$10+$14,$15,$2,$6,"0,0,0",2,$3-$2","$6-$5,0","$5-$2,$0}}'|\
sort-bed - > {output}
"""

rule table_postprocess:
    input: 'results/{genome}/S15.bed'
    output: 
        bed = 'results/{genome}/S16.bed',
        tsv = 'results/{genome}/S16.tsv'
    params:
        max_energy = config['max_energy']
    run:
        t0 = pd.read_table(str(input), header=None)
        t1 = t0.drop(columns=31)
        t1.columns = S16_columns
        t1 = t1.loc[t1.energy < params.max_energy].copy()
        t1['handle1_seq'] = t1.apply(lambda r: foldint_subseq(r['fold_interval_1_seq'], r['handle1_coord']), axis=1)
        t1['handle2_seq'] = t1.apply(lambda r: foldint_subseq(r['fold_interval_2_seq'], r['handle2_coord']), axis=1)
        t1.to_csv(str(output.tsv), index=None, sep='\t')
        t1.to_csv(str(output.bed), index=None, sep='\t', header=None)