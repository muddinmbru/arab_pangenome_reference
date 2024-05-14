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
To evaluate the quality and structural integrity of 53 primary assemblies and 106 haplotype assemblies, we used QUAST (v5.2.0)​62​, which provides metrics such as completeness, N50, and number of contigs. QUAST was run with the following extensive parameter set, as detailed: 

quast.py -o APR043.chm13v2.0 -r chm13v2.0.fa -t 16 APR043.bp.hap1.p_ctg.fa APR043.bp.hap2.p_ctg.fa --large -e 

 
## Yak
The yak suite (v0.1-r66)​24​, which includes yak count, yak qv and yak trioeval, was employed for genome quality validation.  

yak count -t16 -b37 -o APR043.pb.yak APR043.pb.fastq.gz  

was first used to build kmer database and yak qv was used to calculate the QV score with the following parameters:  

-t32 -p -K 3.2g -l 100k APR043.pb.yak APR043.bp.hap1.p_ctg.fa > APR043.hap1.pb.yak.qv.txt 



## Unialigner
All assembled contigs were aligned to the CHM13v2 reference genome using Minimap2 with -L --eqx option. After successful mapping, the chromosome to which each contig was predominantly aligned was identified for further analysis. The centromeric region corresponding to this chromosome was then extracted from the reference genome using SAMtools faidx. To assess the alignment quality of the assembly contig to the extracted centromeric region, we employed UniAligner​31​, specifically using the tandem_aligner command: 

tandem_aligner --first chm13v2.chr1.fa --second $sample{}.1.polished.part_h1tg000053l.fa -o ${sample}.1_h1tg000053l_chr1" 

The output of this alignment is a CIGAR string and each contig's alignment to the centromere of its respective chromosome was further analyzed. A sliding window approach was applied to each CIGAR string, using a window size of 100 base pairs, to calculate the alignment percentage within each window. Finally, the alignment percentages from all windows were plotted to visualize the distribution of alignment quality across different sections of the contigs. 

 ## PStools

To calculate the Hamming rate and switch errors in assemblies constructed with Hi-C data, we used pstools (v0.1)​30​. The phasing errors were assessed employing the following command: 
pstools phasing_error -t 96 $sample.1.polished.fa $sample.2.polished.fa $sample_r1..20fastq $sample_r2.20.fastq 

## Flagger

To assess small-scale assembly errors use flagger.

Concatenate the assemblies into one file per sample and map the HiFi reads to the assembly using either [minimap2](https://github.com/lh3/minimap2) or [Winnowmap](https://github.com/marbl/Winnowmap). 
minimap2 --cs -L -Y -t 32 -ax map-hifi $sample.concat.fa.gz $sample.hifi.fq.gz 
Make sure to sort the resulting sam files. [Samtools](https://github.com/samtools/samtools) can be used for that.

Map the individual assemblies to a reference such as [CHM13v2](https://github.com/marbl/CHM13)
minimap2 -t 64 -L --eqx --cs -ax asm5 chm13v2.0.fa $sample.polished.fa.gz  

The output of the previous files can be used in [flagger](https://github.com/mobinasri/flagger)

## Gene duplication 

We used [Liftoff](https://github.com/agshumate/Liftoff) v1.6.3 (ref.​32​), a tool that accurately maps gene annotations between genome assemblies, to identify gene duplications in our dataset. Liftoff aligns gene sequences from a reference (GENCODE v38) to a target genome and finds the mapping that maximizes sequence identity while preserving the gene structure. An identity threshold of 90% (-sc 0.9) was used, and the following command was used: 

liftoff -p 64 -sc 0.90 -copies -g GENCODE.V38.gff3 -u unmapped.txt -o gene_dup.gff3 -polish {$ASSEMBLY}  {$REFERENCE} 

To remove partial matches from the analysis, we used the -exclude_partial option in Liftoff. We then quantitatively assessed the frequency of copy number variations (CNVs) for each gene in the target genome by comparing the number of gene copies to that in the reference (GENCODE GRCh38.p14 (GENCODE v38)). APR-specific gene duplications were compared against HPRC and CPC gene duplication matrices, as reported in studies utilizing similar methodologies and tool versions. We have constructed a downstream analytical tool PanScan (https://github.com/muddinmbru/panscan) for analyses of complex structural variant analysis. Gene duplication analysis and novel sequence estimation from pangenome graph vcf file.  


# Pangenome

Once the QC is complete, the pangenome can be built using the [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 



# Panscan

The complex regions in the pangenome were assessed using the [Panscan](https://github.com/muddinmbru/panscan) tool. 
You can use the following command once the pipeline is available in your path:

cactus-pangenome ./apr-js-chm13.2908 ./seq_hifiasm.seq --latest --disableCaching --outName apr_review_v1_2902_chm13 --outDir apr_review_v1_2902_chm13 --reference CHM13 GRCh38 --filter 9 --giraffe clip filter --vcf  --chrom-vg clip filter --gbz clip filter full --viz --gfa clip full --vcf --logFile apr-chm13.log --mgCores 160 --mapCores 160 --consCores 160 --indexCores 160 --binariesMode local --workDir work
 


# The data can be found here:

https://www.mbru.ac.ae/the-arab-pangenome-reference/


# Analyses Notebooks

There are three main notebooks:

**Plotting pipeline** was used to plot the bandage plots 
**CSV analysis** was used to analyse complex sites in the pangenome
**gene analysis** was used to anotate genes in the complex regions

The other notebooks are supplemantary and were used to produce plots and stats for the assembly and the pangenome.

Since the analyses is being performed on graph pangenome, the run time for each segment of is around 1 minute to produce the sub-gfa files. 

In case of any questions, please open an issue or contact muhammad.kumail@mbru.ac.ae

