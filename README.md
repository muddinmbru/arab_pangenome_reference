# Arab Pangenome Reference


Comprehensive Resource for the Arab World Delve into the intricate genetic tapestry that defines the Arab people. A Comprehensive Resource for the Arab World built with 106 Assemblies (53 Individuals)

106 Assemblies: Double the insight with assemblies from 53 individuals, ensuring a broad spectrum of genetic representation.
High-Definition Data: Utilizing 30x HiFi coverage and 50x ONT Ultra-long data, we offer unparalleled depth and accuracy in our genome sequences

# Data filtering and QC
Kraken was used to filter out non-human sequences from the assemblies.

# Assembly building
To build the assemblies use Hifiasm (v0.19.5-r603) and run the following command:

If you have HiC data available:
hifiasm -o ${sample} --dual-scaf -t128 --ul ${sample}.filtered.ont.fasta --h1 ${sample}_r1.sampled.fastq --h2 ${sample}_r2.sampled.80.fastq ${sample}.filtered.pb.fastq  

If only Ultra-long and Hifi data is available:
hifiasm -o ${sample} -t128 --ul ${sample}.filtered.ont.fasta ${sample}.filtered.pb.fastq 

Please see for more details:
https://github.com/chhylp123/hifiasm


# Assembly QC
## Quast
To assess some general quality parameters of assembly, such as n50, coverage, run the following command:
quast.py -o APR043.chm13v2.0 -r chm13v2.0.fa -t 16 APR043.bp.hap1.p_ctg.fa APR043.bp.hap2.p_ctg.fa --large -e 


## Flagger

To assess small-scale assembly errors use flagger.

Concatenate the assemblies into one file per sample and map the HiFi reads to the assembly using either [minimap2](https://github.com/lh3/minimap2) or [Winnowmap](https://github.com/marbl/Winnowmap). 
minimap2 --cs -L -Y -t 32 -ax map-hifi $sample.concat.fa.gz $sample.hifi.fq.gz 
Make sure to sort the resulting sam files. [Samtools](https://github.com/samtools/samtools) can be used for that.

Map the individual assemblies to a reference such as [CHM13v2](https://github.com/marbl/CHM13)
minimap2 -t 64 -L --eqx --cs -ax asm5 chm13v2.0.fa $sample.polished.fa.gz  

The output of the previous files can be used in [flagger](https://github.com/mobinasri/flagger)


# Pangenome

Once the QC is complete, the pangenome can be built using the [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 

# Panscan

The complex regions in the pangenome were assessed using the [Panscan](https://github.com/muddinmbru/panscan) tool. 


The data can be found here:

https://www.mbru.ac.ae/the-arab-pangenome-reference/


# Analyses Notebooks

There are three main notebooks:

**Plotting pipeline** was used to plot the bandage plots 
**CSV analysis** was used to analyse complex sites in the pangenome
**gene analysis** was used to anotate genes in the complex regions

The other notebooks are supplemantary and were used to produce plots and stats for the assembly and the pangenome.

Since the analyses is being performed on graph pangenome, the run time for each segment of is around 1 minute to produce the sub-gfa files. 

In case of any questions, please open an issue or contact muhammad.kumail@mbru.ac.ae

