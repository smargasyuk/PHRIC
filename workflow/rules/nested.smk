PREFIX = 'results/'
experiments = ['HepG2', 'H1', 'K562', 'IMR90', 'GM12878', 'HeLa', 'hNPC', 'neuron']
merge_min_radius = 20
merge_max_radius = 100
min_support = 8
merge_min_second_radius = 20
merge_max_second_radius = 50



rule add_file_info:
    input:
        'results/{genome}/{sample}/junctions/RIC-seq_{sample}_rep{rep}/{t}.tsv.gz',
    output:
        PREFIX + '{genome}/S0/{sample}_r{rep}_{t}.tsv.gz'
    shell: """
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $0,"{wildcards.sample}","{wildcards.rep}","{wildcards.t}"}}' |\
pigz - > {output}    
"""

rule merge_junction_files:
    input: expand(PREFIX + '{{genome}}/S0/{sample}_r{rep}_{t}.tsv.gz', sample=experiments, rep=[1,2], t=["Neo", "Chimeric"])
    output: PREFIX + '{genome}/S1.tsv.gz'
    shell: """
unpigz -c {input} | pigz - > {output}    
"""

rule filter_intrachrom:
    input: PREFIX + '{genome}/S1.tsv.gz'
    output: PREFIX + '{genome}/S2.tsv.gz'
    params:
        min_spread = 100
    shell: """
unpigz -c {input} |\
awk -v 'OFS=\t' '$1==$5' |\
awk -v 'OFS=\t' '{{if ($2>$6){{s1=$2;s2=$3;s3=$4;$2=$6;$3=$7;$4=$8;$6=s1;$7=s2;$8=s3}};print}}' |\
awk -v 'OFS=\t' '$6-$2>{params.min_spread}' |\
pigz > {output}
"""

rule count_support:
    input: PREFIX + '{genome}/S2.tsv.gz'
    output: PREFIX + '{genome}/S3.tsv.gz'
    threads: 8
    resources:
        mem_mb=16384
    shell: """
unpigz -c {input} |\
cut -f1,2,3,5,6,7,11 |\
sort --buffer-size={resources.mem_mb}M --parallel={threads} |\
uniq -c |\
sed -E 's/^ *//; s/ /\t/' |\
awk -v 'OFS=\t' '{{print $2,$3,$4,$5,$6,$7,$8,$1}}' |\
pigz > {output}
"""

rule cluster_lr_points:
    input: PREFIX + '{genome}/S3.tsv.gz'
    output: PREFIX + '{genome}/S4.tsv.gz'
    threads: 8
    resources:
        mem_mb=16384
    params:
        clustering_distance = 10
    shell: """
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $0,$0}}' |\
sort-bed --tmpdir {resources.tmpdir} --max-mem {resources.mem_mb}M - |\
bedtools cluster -d {params.clustering_distance} -i stdin |\
cut -f 4- |\
sort-bed --tmpdir {resources.tmpdir} --max-mem {resources.mem_mb}M - |\
bedtools cluster -d {params.clustering_distance} -i stdin |\
cut -f 6- |\
awk -v 'OFS=\t' '{{print $0,"id_"NR}}' |\
sort --parallel={threads}  --buffer-size={resources.mem_mb}M -k9,9 -k10,10 |\
datamash groupby 9,10 first 1 min 2 max 3 first 4 min 5 max 6 unique 7 countunique 7 countunique 11 sum 8 |\
cut -f 3- |\
pigz > {output}
"""

rule build_contact_bed:
    input: PREFIX + '{genome}/S4.tsv.gz'
    output: PREFIX + '{genome}/S4.bed.starch'
    resources:
        mem_mb=16384
    shell: """
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $1,$2,$6,".", ".", "+",$0}}' |\
sort-bed --tmpdir {resources.tmpdir} --max-mem {resources.mem_mb}M - |\
starch - > {output}
"""


rule filter_intragene_annotate:
    input: 
        contacts = PREFIX + '{genome}/S4.bed.starch',
        annotation = PREFIX + '{genome}/annotation/annotation.bed.starch'
    output: PREFIX + '{genome}/S5/{chrom}.tsv.gz'
    shell: """
mkdir -p $(dirname {output})
bedmap --chrom {wildcards.chrom} --fraction-ref 1 --echo --echo-map --delim "\t" --multidelim "\t" --skip-unmapped {input.contacts} {input.annotation} |\
awk -v "OFS=\t" '{{$17=$22;$18=$20; print $0}}' | cut -f7-18 |\
awk -v "OFS=\t" '{{print $0,"id_"NR}}' |\
pigz > {output}
"""

rule get_left_handle:
    input: PREFIX + '{genome}/S5/{chrom}.tsv.gz'
    output: PREFIX + '{genome}/S6/L/{chrom}.bed'
    shell: """
mkdir -p $(dirname {output})
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $1,$2,$3,$13}}' |\
sort-bed - > {output}  
""" 

rule get_right_handle:
    input: PREFIX + '{genome}/S5/{chrom}.tsv.gz'
    output: PREFIX + '{genome}/S6/R/{chrom}.bed'
    shell: """
mkdir -p $(dirname {output})
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $4,$5,$6,$13}}' |\
sort-bed - > {output}  
""" 

rule expand_left_handle:
    input: PREFIX + '{genome}/S6/L/{chrom}.bed'
    output: PREFIX + '{genome}/S6/LI/{chrom}.bed'
    params:
        min_radius = merge_min_radius,
        max_radius = merge_max_radius
    shell: """
mkdir -p $(dirname {output})
cat {input} |\
awk -v 'OFS=\t' '{{$2=$3+{params.min_radius};$3=$3+{params.max_radius};print $0}}' | sort-bed -  > {output}  
""" 

rule expand_right_handle:
    input: PREFIX + '{genome}/S6/R/{chrom}.bed'
    output: PREFIX + '{genome}/S6/RI/{chrom}.bed'
    params:
        min_radius = merge_min_radius,
        max_radius = merge_max_radius
    shell: """
mkdir -p $(dirname {output})
cat {input} |\
awk -v 'OFS=\t' '{{$3=$2-{params.min_radius};$2=$2-{params.max_radius};print $0}}' | sort-bed -  > {output}  
""" 

rule intersect_handles:
    input:
        handles = PREFIX + '{genome}/S6/{t}/{chrom}.bed',
        windows = PREFIX + '{genome}/S6/{t}I/{chrom}.bed'
    output: PREFIX + '{genome}/S7/{t}/{chrom}.tsv.gz'
    threads: 2
    resources:
        mem_mb=4096
    shell: """
mkdir -p $(dirname {output})
intersectBed -wa -wb -sorted \
-a {input.handles} \
-b {input.windows} |\
awk -v 'OFS=\t' '{{print "{wildcards.chrom}_"$4"_"$8,$4,$8}}' |\
sort --parallel={threads}  --buffer-size={resources.mem_mb}M  -T $TEMPDIR |\
pigz - > {output}
"""

rule merge_is_handles:
    input:
        right = PREFIX + '{genome}/S7/R/{chrom}.tsv.gz',
        left = PREFIX + '{genome}/S7/L/{chrom}.tsv.gz'
    output: PREFIX + '{genome}/S7/B/{chrom}.tsv'
    shell: """
mkdir -p $(dirname {output})
comm -12 <(unpigz -c {input.left}) <(unpigz -c {input.right}) | awk -v 'OFS=\t' '$2!=$3' > {output}
"""

rule annotate_pairs:
    input:
        pairs = PREFIX + '{genome}/S7/B/{chrom}.tsv',
        contacts = PREFIX + '{genome}/S5/{chrom}.tsv.gz'
    output: 
        l = PREFIX + '{genome}/S8/S1/{chrom}.tsv',
        r = PREFIX + '{genome}/S8/S2/{chrom}.tsv'
    shell: """
mkdir -p $(dirname {output})
join -1 2 -2 13 -t$'\t' <(sort -k2,2 {input.pairs}) <(unpigz -c {input.contacts} | sort -k13,13) > {output.l}
join -1 3 -2 13 -t$'\t' <(sort -k3,3 {output.l})    <(unpigz -c {input.contacts} | sort -k13,13) | cut -f 4- > {output.r}
"""

rule merge_pairs:
    input: lambda wildcards: expand(PREFIX + '{{genome}}/S8/S2/{chrom}.tsv', chrom=get_chroms(wildcards))
    output: PREFIX + '{genome}/S9.tsv'
    shell: """
awk -v 'OFS=\t' '{{print $0,"id_"NR}}' {input} > {output}    
"""

rule filter_pairs_get_intercontacts:
    input:  PREFIX + '{genome}/S9.tsv'
    output: PREFIX + '{genome}/S10.tsv'
    params: 
        min_support = min_support
    shell: """
cat {input} |\
awk -v 'OFS=\t' '($10>{params.min_support})&&($22>{params.min_support})' |\
awk -v 'OFS=\t' '{{print $1,$15,$2,$13,$6,$17,$7,$8,$9,$10,$19,$20,$21,$22,$11,$12,$25}}' > {output}
"""


rule filter_intercontacts:
    input: 
        pairs = PREFIX + '{genome}/S10.tsv',
        rmsk = 'resources/genomes/{genome}/RMSK.bed'
    output: PREFIX + '{genome}/S11.tsv'
    params:
        min_radius = merge_min_second_radius,
        max_radius = merge_max_second_radius
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
        pairs = PREFIX + '{genome}/S11.tsv',
        genome_fasta = 'resources/genomes/{genome}/genome.fa'
    output:
        seq_L = PREFIX + '{genome}/S12L.fa',
        seq_R = PREFIX + '{genome}/S12R.fa',
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
        seq_L = PREFIX + '{genome}/S12L.fa',
        seq_R = PREFIX + '{genome}/S12R.fa',
    output: PREFIX + '{genome}/S14_PrePH.tsv'
    conda: "../envs/PREPH.yaml"
    threads: 16
    shell: """
paste -d '\n' <(grep -v '>' {input.seq_L}) <(grep -v '>' {input.seq_R}) |\
grep -v "^$" |  xargs -n2 |\
parallel 2>/dev/null --keep-order -j{threads} --colsep ' ' "python workflow/scripts/PrePH/src/fold_SM.py -f {{1}} -s {{2}} \
-k 3 -a 3 -e -1 -u False -d 2" > {output}   
"""

rule merge_w_preph_results:
    input: 
        pairs = PREFIX + '{genome}/S11.tsv',
        preph = PREFIX + '{genome}/S14_PrePH.tsv',
        seq_L = PREFIX + '{genome}/S12L.fa',
        seq_R = PREFIX + '{genome}/S12R.fa',
    output: PREFIX + '{genome}/S15.bed'
    shell: """
paste {input.pairs} {input.preph} <(grep -v '>' {input.seq_L}) <(grep -v '>' {input.seq_R}) |\
awk -v 'OFS=\t' '{{gsub(/[\(\)]/, "", $22)}}1' |\
awk -v 'OFS=\t' 'NF==24' |\
awk -v 'OFS=\t' '{{print $1,$2,$6,$17,$10+$14,$15,$2,$6,"0,0,0",2,$3-$2","$6-$5,0","$5-$2,$0}}'|\
sort-bed - > {output}
"""


# rule all:
#     input: 'results/hg19/S11.tsv' 


# rule all:
#     input: expand('results/hg19/S12{t}.fa', t=['L', 'R'])

# rule all:
#     input: 'results/hg19/S14_PrePH.tsv' 

rule all:
    input: 'results/hg19/S15.bed' 