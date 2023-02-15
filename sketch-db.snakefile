import pandas as pd


out_dir = 'output.database'
logs_dir = os.path.join(out_dir, "logs")
taxinfo = "inputs/ralstonia-lin.taxonomy.with-picklist-options.csv"
tx = pd.read_csv(taxinfo, sep = ',')

tx['ASM'] = tx['FileName'].str.split('_', expand=True)[2]
tx['name'] = tx['Species'] + " strain=" + tx["Strain"] + ', ' + tx['ASM']
tx['genome_filename'] = "inputs/genomes/" + tx['Accession'] + "_genomic.fna.gz"
tx['name'] = tx['Accession']
tx['protein_filename'] = ''

ACC = tx['Accession'].tolist()
print(ACC)



rule all:
    input: expand("{out_dir}/ralstonia.dna.sc{scaled}.zip", out_dir=out_dir, scaled=[1,5,10,100])

rule write_fromfile:
    output: "inputs/genomes/ralstonia.fromfile.csv"
    run:
        ffDF = tx[["name", "genome_filename", "protein_filename"]]
        ffDF.to_csv(str(output), index=False)


rule sourmash_sketch_dna:
    input: 
        ff="inputs/genomes/ralstonia.fromfile.csv",
        genomes=expand("inputs/genomes/{acc}_genomic.fna.gz", acc=ACC)
    output:
        os.path.join(out_dir, "ralstonia.dna.sc{scaled}.zip")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sketch", "sketch_dna.sc{scaled}.log")
    benchmark: os.path.join(logs_dir, "sketch_dna.sc{scaled}.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch fromfile {input.ff} -p k=21,k=31,k=51,dna,scaled={wildcards.scaled},abund \
                                 -o {output} 2> {log}
        """



