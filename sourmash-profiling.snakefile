import os, sys
import pandas as pd
import glob

configfile: "inputs/ralstonia.dna.conf"

out_dir = config.get('output_dir', 'output.sourmash-profiling')
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")

sample_info = pd.read_csv(config['sample_info'])
SAMPLES = sample_info["name"].to_list()
# set name as index for easy access
sample_info.set_index('name', inplace=True)

search_databases = config['search_databases'] # must be dictionary
ksize = config.get("ksize", [51])
if not isinstance(ksize, list):
    ksize=[ksize]
scaled = config.get("scaled", [1000])
if not isinstance(scaled, list):
    scaled=[scaled]


onstart:
    print("------------------------------")
    print("sourmash taxonomic profiling workflow")
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


rule all:
    input:
        expand(os.path.join(out_dir, 'gather', '{sample}.k{ks}-sc{scaled}.gather.kreport.txt'), sample=SAMPLES, ks=ksize, scaled=scaled),
#        expand(os.path.join(out_dir, 'gather', '{sample}.k{ks}.gather.with-lineages.csv'),sample=SAMPLES, ks=ksize),


rule sourmash_sketch_dna:
    input: ancient(lambda w: sample_info.at[w.sample, "reads"]) 
    output:
        os.path.join(out_dir, "reads", "{sample}.dna.sig.zip")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sketch", "{sample}.sketch_dna.log")
    benchmark: os.path.join(benchmarks_dir, "sketch", "{sample}.sketch_dna.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch dna {input} -p k=21,k=31,k=51,dna,scaled=5,abund \
                                    --name {wildcards.sample} -o {output} 2> {log}
        """

rule sourmash_gather:
    input:
        query=os.path.join(out_dir, "reads", "{sample}.dna.sig.zip"),
        databases = lambda w: search_databases[f"k{w.ksize}"],
    output:
        gather_csv=os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}.gather.txt'),
    params:
        threshold_bp = config.get('sourmash_database_threshold_bp', '50000'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        time=10000,
        partition="bmh",
    log: os.path.join(logs_dir, "gather", "{sample}.k{ksize}-sc{scaled}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{sample}.k{ksize}-sc{scaled}.gather.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {input.databases}"
        echo "DB(s): {input.databases}" > {log}

        sourmash gather {input.query} {input.databases} --dna --ksize {wildcards.ksize} \
                 --threshold-bp {params.threshold_bp} --scaled {wildcards.scaled} \
                 -o {output.gather_csv} > {output.gather_txt} 2>> {log}
        
        touch {output.gather_txt}
        touch {output.gather_csv}
        """

rule tax_metagenome:
    input:
        gather = os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}.gather.krona.tsv'),
        os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}.gather.summarized.csv'),
        os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}.gather.kreport.txt'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "bml",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'gather'),
        out_base= lambda w: f'{w.sample}.k{w.ksize}-sc{scaled}.gather',
    conda: "conf/env/sourmashLIN.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome -g {input.gather} -t {input.lineages} -o {params.out_base} \
                                --output-dir {params.outd} --output-format krona csv_summary kreport \
                                --LIN-taxonomy --LIN-position 19
        """

#rule tax_annotate:
#    input:
#        gather = os.path.join(out_dir, 'gather', '{sample}.k{ksize}.gather.csv'),
#        lineages = config['database_lineage_files'],
#    output:
#        os.path.join(out_dir, 'gather', '{sample}.k{ksize}.gather.with-lineages.csv'),
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *3000,
#        time=240,
#        partition = "low2",
#    params:
#        outd= lambda w: os.path.join(out_dir, f'gather'),
#    conda: "conf/env/sourmashLIN.yml"
#    shell:
#        """
#        mkdir -p {params.outd}
#        sourmash tax annotate -g {input.gather} -t {input.lineages} -o {params.outd}
#        """


