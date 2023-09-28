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
    params:
        preph_parameters = config['preph_parameters']
    threads: 16
    shell: """
paste  <(grep -v '>' {input.seq_L}) <(grep -v '>' {input.seq_R}) |\
grep -v "^$" |\
python workflow/scripts/PrePH/src/fold2.py {params.preph_parameters} -j{threads} \
> {output}
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
    run:
        t0 = pd.read_table(str(input), header=None)
        t3 = postprocess_preph_table(t0)
        t3.to_csv(str(output.tsv), index=None, sep='\t')
        t3.to_csv(str(output.bed), index=None, sep='\t', header=None)