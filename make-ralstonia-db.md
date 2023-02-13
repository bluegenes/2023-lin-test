# Make this small Ralstonia database

## Extract existing signatures from GenBank/GTDB

26 of 27 are available. In `GenBank`, they have `GCA` accessions ; in `GTDB`, `GCF. First make some additional columns to allow matching either GCF or GCA (done in pandas, output file: `ralstonia-lin.taxonomy.with-picklist-options.csv`).

Extract from GenBank
```
sourmash sig extract --picklist ralstonia-lin.taxonomy.with-picklist-options.csv:identprefix:identprefix /group/ctbrowngrp/sourmash-db/genbank-2022.03/genbank-2022.03-bacteria-k21.zip --ksize 21 -o ralstonia-k21.from-genbank.zip
```

output:
```
== This is sourmash version 4.6.1. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

picking column 'identprefix' of type 'identprefix' from 'ralstonia-lin.taxonomy.with-picklist-options.csv'
loaded 27 distinct values into picklist.
loaded 26 total that matched ksize & molecule type
extracted 26 signatures from 1 file(s)
for given picklist, found 26 matches to 27 distinct values
WARNING: 1 missing picklist values.
```

> Note: could alternatively extract from gtdb:
> GTDB:
> ```sourmash sig extract --picklist ralstonia-lin.taxonomy.with-picklist-options.csv:Accession:ident /group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.genomic.k51.zip --ksize 51 -o ralstonia-k51.from-gtdb.zip
> ```

Repeated the extractions for k=21, k=31 just to have all ksizes available.

## Find the missing genome

```
sourmash sig check ralstonia-k51.from-genbank.zip --picklist ralstonia-lin.taxonomy.with-picklist-options.csv:identprefix:identprefix -o missing.genbank.txt
```

In `missing.genbank.txt`, I see: 

```
LIN,Species,Strain,FileName,Accession,ident,identprefix
14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;1;0;5;0;0,Ralstonia solanacearum,CIP238_UW490,GCF_023076135.1_ASM2307613v1_genomic.fna,GCF_023076135.1,GCA_023076135.1,GCA_023076135
```

## Download and sketch missing genome 

Date Added: 2022/04/20, which explains why it's not in our March 2022 GenBank db or GTDB-rs207.
Info: https://www.ncbi.nlm.nih.gov/search/all/?term=GCA_023076135


RefSeq accession (`GCF`) is suppressed, so download the `GenBank` instead. For the rest of the signatures, use sigs extracted from `GenBank` rather than `GTDB` so all accession will start with `GCA`.

# sketch missing genome using same name format as other signatures
```
sourmash sketch dna GCA_023076135.1_ASM2307613v1_genomic.fna.gz -p k=21,k=31,k=51,scaled=1000,abund --name "GCA_023076135.1 Ralstonia solanacearum strain=CIP238_UW490, ASM2307613v1" -o GCA_023076135.sig.zip
```

# cat all signatures into a single Ralstonia database
```
sourmash sig cat ralstonia-k*.from-genbank.zip GCA_023076135.sig.zip -o ralstonia.GCA.zip
```

output:
```
== This is sourmash version 4.6.1. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

loaded 81 signatures total, from 4 files
loaded 81 signatures total.
output 81 signatures
```

sanity check: 27 genomes * 3 ksizes = 81 sigs.

sourmash check:
```
sourmash sig fileinfo ralstonia.GCA.zip

== This is sourmash version 4.6.1. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

** loading from 'ralstonia.GCA.zip'
path filetype: ZipFileLinearIndex
location: /home/ntpierce/2023-lin-test/ralstonia.GCA.zip
is database? yes
has manifest? yes
num signatures: 81
** examining manifest...
total hashes: 445041
summary of sketches:
   27 sketches with DNA, k=21, scaled=1000, abund     148324 total hashes
   27 sketches with DNA, k=31, scaled=1000, abund     148111 total hashes
   27 sketches with DNA, k=51, scaled=1000, abund     148606 total hashes
```
