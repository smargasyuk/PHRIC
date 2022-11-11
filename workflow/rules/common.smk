from pathlib import Path

def get_chroms(wildcards):
    with open(f'resources/genomes/{wildcards.genome}/chromSizes') as f:
        ls = [l.strip() for l in f.readlines() if l.strip()]
        return [l.split('\t')[0] for l in ls]