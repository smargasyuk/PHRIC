rule annotation_gtf_to_bed:
    input: 'resources/genomes/{genome}/annotation.gtf.gz'
    output: 
        gp1 = 'results/{genome}/annotation/annotation.gp',
        gp2 = 'results/{genome}/annotation/annotation.gene_name.gp',
        bed1 = 'results/{genome}/annotation/annotation.unsorted.bed',
        starch = 'results/{genome}/annotation/annotation.bed.starch'
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output.starch})
gtfToGenePred {input} {output.gp1} -geneNameAsName2 -ignoreGroupsWithoutExons -genePredExt
awk '{{print $12"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"}}' {output.gp1} > {output.gp2}
genePredToBed {output.gp2} {output.bed1}
sort-bed {output.bed1} | starch -  > {output.starch}
"""


rule add_file_info:
    input: "resources/RNAcontacts-output/{genome}/{project}/junctions/{sample}/{t}.tsv"
    output:
        'results/{genome}/S0/{sample}_{t}.tsv.gz'
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $0,"{wildcards.sample}","{wildcards.t}"}}' |\
pigz - > {output}    
"""

rule merge_junction_files:
    input: get_junctions
    output: 'results/{genome}/S1.tsv.gz'
    conda: "../envs/common.yaml"
    shell: """
unpigz -c {input} | pigz - > {output}    
"""

rule filter_intrachrom:
    input: 'results/{genome}/S1.tsv.gz'
    output: 'results/{genome}/S2.tsv.gz'
    params:
        min_spread = 2 * config['distances']['nested_pass1']['rmax']
    conda: "../envs/common.yaml"
    shell: """
unpigz -c {input} |\
awk -v 'OFS=\t' '$1==$5' |\
awk -v 'OFS=\t' '{{if ($2>$6){{s1=$2;s2=$3;s3=$4;$2=$6;$3=$7;$4=$8;$6=s1;$7=s2;$8=s3}};print}}' |\
awk -v 'OFS=\t' '$6-$2>{params.min_spread}' |\
pigz > {output}
"""

rule count_support:
    input: 'results/{genome}/S2.tsv.gz'
    output: 'results/{genome}/S3.tsv.gz'
    threads: 8
    resources:
        mem_mb=16384
    conda: "../envs/common.yaml"
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
    input: 'results/{genome}/S3.tsv.gz'
    output: 'results/{genome}/S4.tsv.gz'
    threads: 8
    resources:
        mem_mb=16384
    params:
        clustering_distance = config['distances']['cluster']
    conda: "../envs/common.yaml"
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
    input: 'results/{genome}/S4.tsv.gz'
    output: 'results/{genome}/S4.bed.starch'
    resources:
        mem_mb=16384
    conda: "../envs/common.yaml"
    shell: """
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $1,$2,$6,".", ".", "+",$0}}' |\
sort-bed --tmpdir {resources.tmpdir} --max-mem {resources.mem_mb}M - |\
starch - > {output}
"""


rule filter_intragene_annotate:
    input: 
        contacts = 'results/{genome}/S4.bed.starch',
        annotation = 'results/{genome}/annotation/annotation.bed.starch'
    output: 'results/{genome}/S5/{chrom}.tsv.gz'
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
bedmap --chrom {wildcards.chrom} --fraction-ref 1 --echo --echo-map --delim "\t" --multidelim "\t" --skip-unmapped {input.contacts} {input.annotation} |\
awk -v "OFS=\t" '{{$17=$22;$18=$20; print $0}}' | cut -f7-18 |\
awk -v "OFS=\t" '{{print $0,"id_"NR}}' |\
pigz > {output}
"""

rule get_left_handle:
    input: 'results/{genome}/S5/{chrom}.tsv.gz'
    output: 'results/{genome}/S6/L/{chrom}.bed'
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $1,$2,$3,$13}}' |\
sort-bed - > {output}  
""" 

rule get_right_handle:
    input: 'results/{genome}/S5/{chrom}.tsv.gz'
    output: 'results/{genome}/S6/R/{chrom}.bed'
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
unpigz -c {input} |\
awk -v 'OFS=\t' '{{print $4,$5,$6,$13}}' |\
sort-bed - > {output}  
""" 

rule expand_left_handle:
    input: 'results/{genome}/S6/L/{chrom}.bed'
    output: 'results/{genome}/S6/LI/{chrom}.bed'
    params:
        min_radius = config['distances']['nested_pass1']['rmin'],
        max_radius = config['distances']['nested_pass1']['rmax'],
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
cat {input} |\
awk -v 'OFS=\t' '{{$2=$3+{params.min_radius};$3=$3+{params.max_radius};print $0}}' | sort-bed -  > {output}  
""" 

rule expand_right_handle:
    input: 'results/{genome}/S6/R/{chrom}.bed'
    output: 'results/{genome}/S6/RI/{chrom}.bed'
    params:
        min_radius = config['distances']['nested_pass1']['rmin'],
        max_radius = config['distances']['nested_pass1']['rmax'],
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
cat {input} |\
awk -v 'OFS=\t' '{{$3=$2-{params.min_radius};$2=$2-{params.max_radius};print $0}}' | sort-bed -  > {output}  
""" 

rule intersect_handles:
    input:
        handles = 'results/{genome}/S6/{t}/{chrom}.bed',
        windows = 'results/{genome}/S6/{t}I/{chrom}.bed'
    output: 'results/{genome}/S7/{t}/{chrom}.tsv.gz'
    threads: 2
    resources:
        mem_mb=4096
    conda: "../envs/common.yaml"
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
        right = 'results/{genome}/S7/R/{chrom}.tsv.gz',
        left = 'results/{genome}/S7/L/{chrom}.tsv.gz'
    output: 'results/{genome}/S7/B/{chrom}.tsv'
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
comm -12 <(unpigz -c {input.left}) <(unpigz -c {input.right}) | awk -v 'OFS=\t' '$2!=$3' > {output}
"""

rule annotate_pairs:
    input:
        pairs = 'results/{genome}/S7/B/{chrom}.tsv',
        contacts = 'results/{genome}/S5/{chrom}.tsv.gz'
    output: 
        l = 'results/{genome}/S8/S1/{chrom}.tsv',
        r = 'results/{genome}/S8/S2/{chrom}.tsv'
    conda: "../envs/common.yaml"
    shell: """
mkdir -p $(dirname {output})
join -1 2 -2 13 -t$'\t' <(sort -k2,2 {input.pairs}) <(unpigz -c {input.contacts} | sort -k13,13) > {output.l}
join -1 3 -2 13 -t$'\t' <(sort -k3,3 {output.l})    <(unpigz -c {input.contacts} | sort -k13,13) | cut -f 4- > {output.r}
"""

rule merge_pairs:
    input: lambda wildcards: expand('results/{{genome}}/S8/S2/{chrom}.tsv', chrom=get_chroms(wildcards))
    output: 'results/{genome}/S9.tsv'
    conda: "../envs/common.yaml"
    shell: """
awk -v 'OFS=\t' '{{print $0,"id_"NR}}' {input} > {output}    
"""
