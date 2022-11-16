# General settings

To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the following columns have to be defined:

+ `sample_name`: unique sample name;
+ `project`: project name for a group of RIC-seq experiments and control RNA-seq;
+ `genome`: genome name. Used to select the set of annotations for the sample:
  + `resources/genomes/{genome}/genome.fa`: genome sequence fasta;
  + `resources/genomes/{genome}/chromSizes`: two-column table with chromosome names and sizes;
  + `resources/genomes/{genome}/annotation.gtf.gz`: transcript annotation for the genome;
  + `resources/genomes/{genome}/RMSK.bed`: repeats annotation in BED format;
+ `treatment`: 'experiment' or 'control';
+ `fq1`, `fq2`: paired-end fastq read files.
