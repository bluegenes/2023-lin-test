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

threshold_bp = config.get('sourmash_gather_threshold_bp', [50000])
if not isinstance(threshold_bp, list):
    threshold_bp=[threshold_bp]
print(threshold_bp)



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
        expand(os.path.join(out_dir, 'gather', '{sample}.k{ks}-sc{sc}-thr{thr}.gather.kreport.txt'), sample=SAMPLES, ks=ksize, sc=scaled, thr=threshold_bp),
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
        #databases = lambda w: search_databases[f"k{w.ksize}"],
    output:
        gather_csv=os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.txt'),
    params:
        # use threshold as wildcard instead
        #threshold_bp = config.get('sourmash_database_threshold_bp', '5000'),
        # since I actually made diff scaled databases for this test set, let's just load 
        # the exact resolution we need instead of downsampling the lower scaled dbs
        database = lambda w: [x for x in search_databases[f"k{w.ksize}"] if f".sc{w.scaled}." in x],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        time=10000,
        partition="bmm",
    log: os.path.join(logs_dir, "gather", "{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {params.database}"
        echo "DB(s): {params.database}" > {log}

        sourmash gather {input.query} {params.database} --dna --ksize {wildcards.ksize} \
                 --threshold-bp {wildcards.thresh} --scaled {wildcards.scaled} \
                 -o {output.gather_csv} > {output.gather_txt} 2>> {log}
        
        touch {output.gather_txt}
        touch {output.gather_csv}
        """

rule tax_metagenome:
    input:
        gather = os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.krona.tsv'),
        os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.summarized.csv'),
        os.path.join(out_dir, 'gather', '{sample}.k{ksize}-sc{scaled}-thr{thresh}.gather.kreport.txt'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "bml",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'gather'),
        out_base= lambda w: f'{w.sample}.k{w.ksize}-sc{scaled}-thr{w.thresh}.gather',
    log: os.path.join(logs_dir, "tax", "{sample}.k{ksize}-sc{scaled}-thr{thresh}.tax-metagenome.log")
    benchmark: os.path.join(benchmarks_dir, "tax", "{sample}.k{ksize}-sc{scaled}-thr{thresh}.tax-metagenome.benchmark")
    conda: "conf/env/sourmashLIN.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome -g {input.gather} -t {input.lineages} -o {params.out_base} \
                                --output-dir {params.outd} --output-format krona csv_summary kreport \
                                --LIN-taxonomy --LIN-position 19 2> {log}
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


