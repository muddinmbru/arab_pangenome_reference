# Arab Pangenome Reference


A Comprehensive Resource for the Arab World built with 106 Assemblies (53 Individuals)

106 Assemblies: Double the insight with assemblies from 53 individuals, ensuring a broad spectrum of genetic representation.
High-Definition Data: Utilizing 30x HiFi coverage and 50x ONT Ultra-long data, we offer unparalleled depth and accuracy in our genome sequences
18 samples with high-coverage Hi-C data. 

# Data filtering and QC
Kraken was used to filter out non-human sequences from the assemblies.

# Assembly building
To build the assemblies use Hifiasm (v0.19.5-r603) and run the following command:

If you have HiC data available:
```
hifiasm -o ${sample} --dual-scaf -t128 --ul ${sample}.filtered.ont.fasta --h1 ${sample}_r1.sampled.fastq --h2 ${sample}_r2.sampled.80.fastq ${sample}.filtered.pb.fastq  
```
If only Ultra-long and Hifi data is available:
```
hifiasm -o ${sample} -t128 --ul ${sample}.filtered.ont.fasta ${sample}.filtered.pb.fastq 
```
Please see for more details:
https://github.com/chhylp123/hifiasm


# Assembly QC

## Quast
To evaluate the quality and structural integrity of 53 primary assemblies and 106 haplotype assemblies, we used QUAST (v5.2.0)​62​, which provides metrics such as completeness, N50, and number of contigs. QUAST was run with the following extensive parameter set, as detailed: 
```
quast.py -o APR043.chm13v2.0 -r chm13v2.0.fa -t 16 APR043.bp.hap1.p_ctg.fa APR043.bp.hap2.p_ctg.fa --large -e 
```

 
## Yak
The yak suite (v0.1-r66)​24​, which includes yak count, yak qv and yak trioeval, was employed for genome quality validation.  

```
yak count -t16 -b37 -o APR043.pb.yak APR043.pb.fastq.gz  
```

was first used to build kmer database and yak qv was used to calculate the QV score with the following parameters:  

```
-t32 -p -K 3.2g -l 100k APR043.pb.yak APR043.bp.hap1.p_ctg.fa > APR043.hap1.pb.yak.qv.txt 
```



## Unialigner
All assembled contigs were aligned to the CHM13v2 reference genome using Minimap2 with -L --eqx option. After successful mapping, the chromosome to which each contig was predominantly aligned was identified for further analysis. The centromeric region corresponding to this chromosome was then extracted from the reference genome using SAMtools faidx. To assess the alignment quality of the assembly contig to the extracted centromeric region, we employed UniAligner​31​, specifically using the tandem_aligner command: 

```
tandem_aligner --first chm13v2.chr1.fa --second $sample{}.1.polished.part_h1tg000053l.fa -o ${sample}.1_h1tg000053l_chr1" 
```

The output of this alignment is a CIGAR string and each contig's alignment to the centromere of its respective chromosome was further analyzed. A sliding window approach was applied to each CIGAR string, using a window size of 100 base pairs, to calculate the alignment percentage within each window. Finally, the alignment percentages from all windows were plotted to visualize the distribution of alignment quality across different sections of the contigs. 

 ## PStools

To calculate the Hamming rate and switch errors in assemblies constructed with Hi-C data, we used pstools (v0.1)​30​. The phasing errors were assessed employing the following command: 

```
pstools phasing_error -t 96 $sample.1.polished.fa $sample.2.polished.fa $sample_r1..20fastq $sample_r2.20.fastq 
```

## Flagger

To assess small-scale assembly errors use flagger.

Concatenate the assemblies into one file per sample and map the HiFi reads to the assembly using either [minimap2](https://github.com/lh3/minimap2) or [Winnowmap](https://github.com/marbl/Winnowmap). 
```
minimap2 --cs -L -Y -t 32 -ax map-hifi $sample.concat.fa.gz $sample.hifi.fq.gz 
```
Make sure to sort the resulting sam files. [Samtools](https://github.com/samtools/samtools) can be used for that.

Map the individual assemblies to a reference such as [CHM13v2](https://github.com/marbl/CHM13)
```
minimap2 -t 64 -L --eqx --cs -ax asm5 chm13v2.0.fa $sample.polished.fa.gz
```

The output of the previous files can be used in [flagger](https://github.com/mobinasri/flagger)

## Gene duplication 

We used [Liftoff](https://github.com/agshumate/Liftoff) v1.6.3 (ref.​32​), a tool that accurately maps gene annotations between genome assemblies, to identify gene duplications in our dataset. Liftoff aligns gene sequences from a reference (GENCODE v38) to a target genome and finds the mapping that maximizes sequence identity while preserving the gene structure. An identity threshold of 90% (-sc 0.9) was used, and the following command was used: 

```
liftoff -p 64 -sc 0.90 -copies -g GENCODE.V38.gff3 -u unmapped.txt -o gene_dup.gff3 -polish {$ASSEMBLY}  {$REFERENCE} 
```

To remove partial matches from the analysis, we used the -exclude_partial option in Liftoff. We then quantitatively assessed the frequency of copy number variations (CNVs) for each gene in the target genome by comparing the number of gene copies to that in the reference (GENCODE GRCh38.p14 (GENCODE v38)). APR-specific gene duplications were compared against HPRC and CPC gene duplication matrices, as reported in studies utilizing similar methodologies and tool versions. We have constructed a downstream analytical tool PanScan (https://github.com/muddinmbru/panscan) for analyses of complex structural variant analysis. Gene duplication analysis and novel sequence estimation from pangenome graph vcf file.  


# Pangenome

Once the QC is complete, the pangenome can be built using the [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 
You can use the following command once the pipeline is available in your path:

```
cactus-pangenome ./apr-js-chm13.2908 ./seq_hifiasm.seq --latest --disableCaching --outName apr_review_v1_2902_chm13 --outDir apr_review_v1_2902_chm13 --reference CHM13 GRCh38 --filter 9 --giraffe clip filter --vcf  --chrom-vg clip filter --gbz clip filter full --viz --gfa clip full --vcf --logFile apr-chm13.log --mgCores 160 --mapCores 160 --consCores 160 --indexCores 160 --binariesMode local --workDir work
```


# Panscan

The complex regions in the pangenome were assessed using the [Panscan](https://github.com/muddinmbru/panscan) tool. 


## Mitochondrial pangenome construction 

To construct a mitochondrial Arab pangenome (mtAPR) that captures the diversity of Arab mitochondrial DNA, we used high-quality long reads from 53 individuals. We first mapped the reads to ChrM of the CHM13v2 reference genome using minimap2 (v2.26)​69​ with 90% similarity and retained only reads that were longer than 15 kb; this threshold was set to substantially reduce the chances of inadvertent nuclear DNA contamination. This resulted in a total of 20,520 reads (19,251 ONT reads and 1,385 HiFi reads). We used PGGB (v0.5.4)​70​ to construct a mitochondrial pangenome graph from the mitochondrial contigs of HiFi reads of all individuals, and each read of >15 kb was concatenated in one fasta file along with the CHM13v2 mitochondrial chromosome. [PGGB](https://github.com/pangenome/pggb) was then run on these samples with the following command: 
```
pggb -i pggb.in.90.no_dup.chm13.fasta.gz -o output_chm13_local -n 53 -t 90 -p 90 -s 100 -V 'chm13v2:#:50' 
```

We visualized the mtAPR graph using Bandage NG version (v2022.09)​67​ and displayed the nodes, edges and variants in different colors and shapes. 

 

## Mitochondrial pangenome variant identification 

To process the VCF files of the mtAPR pangenome, we employed a multistage approach to ensure data accuracy and integrity. The initial step involved segregation of multiallelic variant sites into biallelic records, which was accomplished using the BCFtools normalization function (BCFtools norm -m). Next, we implemented the RTG tool vcfdecompose (version 3.12.1) with the parameters --break-mnps and --break-indels to further resolve complex variants into their constituent single-nucleotide polymorphisms (SNPs), indels and SVs. To merge genotype information, identical variants identified across the dataset were consolidated into singular records. This step was carried out using a custom Perl script. Parallel to the APR data, the HPRC mitochondrial pangenome VCF file was subjected to the same rigorous processing methodology to maintain consistency across datasets. The small variants (<10 bp) were further filtered to obtain those that were concordant with DeepVariant calls of the mitochondrial reads. To identify novel variants unique to the APR mitochondrial pangenome, we systematically removed variants that were shared with the HPRC mitochondrial pangenome. We then extracted insertions of 10 base pairs or more that were subsequently clustered using the cd-hit-est program, which applied a stringent 90% similarity threshold (-c 0.9) to discern and characterize novel insertions. 

# Population Structure Analysis

## ADMIXTURE

To run ADMIXTURE analysis on our cohort by first variant calling it on a set of variants found in the Human Origins array dataset downloaded from the Allen Ancient DNA Resource V54.1.p1 using the command 
```
bcftools mpileup -f ${REFGEN} -I -E -T ${VCF} -r chr{i} -b {bam.list} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT} -f GQ
```
The joint called APR VCF was merged with a set of 1040 samples from the Human Origins (HO) dataset. The merged file was filtered for genotype quality (GQ) greater than 20, minor allele frequency (MAF) greater than 0.05, and linkage disequilibrium pruning using the following command: 
```
plink --bfile data --maf 0.05 --indep-pairwise 50 50 0.5 --make-bed --out file 
```
and then ADMIXTURE was run on the filtered file with the following command
```
admixture PCA.file.bed {5..9) —cv -j 20
```

## PCA
For the PCA we used the R library SNPRelate (v1.28.0). We used the same merged and filtered dataset from the ADMIXTURE methodology above. The Rscript to plot the PCA from the vcf file is provided below
```R
library("SNPRelate")
library("rgl")

vcf1.fn<-"{$merged.vcf}"
snpgdsVCF2GDS(vcf1.fn, "{$merged.gds}",  method="biallelic.only")

genofile<- openfn.gds("{$merged.gds}")
ccm_pca<-snpgdsPCA(genofile, num.thread=64)

#Read in your pedigree file showing sample ethnicities for coloring the plot
PED <- read.table('{$samples.ped}', header = TRUE, skip = 0, sep = '\t')
PED <- PED[which(PED$Individual.ID %in% ccm_pca$sample.id), ]
PED <- PED[match(ccm_pca$sample.id, PED$Individual.ID),]
#Making sure all your samples in the dataset match with the pedigree file
all(PED$Individual.ID == ccm_pca$sample.id) == TRUE

#Match population ids with colors you want
PED$Population <- factor(PED$Population, levels=c(
  "ethnicity1-a","ethnicity1-b","ethnicity1-c",
  "ethnicity2-a","ethnicity2-b","ethnicity2-c","ethnicity2-d",
  "ethnicity3-a","ethnicity3-b",
  ))

col <- colorRampPalette(c(
  "yellow","yellow","yellow",
  "forestgreen","forestgreen","forestgreen","forestgreen",
  "red","red"))(length(unique(PED$Population)))[factor(PED$Population)]

project.pca <- ccm_pca$eigenvect
png('pca.png', 1200,1800)
par(mar=c(11,11,4,1)+.1)
plot(project.pca[,1], project.pca[,2],
  type = 'n',
  main = '',
  adj = 0.5,
  xlab = '',
  ylab = '',
  cex.xlab= 10,
  cex.ylab= 10,
  font = 4,
  font.lab = 4,
  cex.lab=4,cex.main=1, cex.axis=2, cex.sub=2)
points(project.pca[,1], project.pca[,2], col = col, pch = 20, cex = 4)
legend('topleft',
  bty = 'n',
  cex = 3,
  title = '',
  c('Ethnicity1', 'Ethnicity2', 'Ethnicity3'),
  fill = c('yellow' , "forestgreen" , "red"))
```

##fineSTRUCTURE
To confirm the unrelated status of all 50 APR samples, we constructed a heatmap using haplotype sharing variance obtained from fineSTRUCTURE​​ runs. The variants called using DeepVariant were filtered for a GQ of ≥20. We then reduced the dataset to one variant per 10 kb window and jointly phased it with Eagle using the following command 
```
for i in {1..22}; do
       eagle --vcfTarget=chr$i.vcf.gz --vcfRef=1000genomes-all-chr-chr$i.vcf.gz --outPrefix=phased_chr$i --geneticMapFile=Eagle_v2.4.1/tables/chrom_genetic_maps/genetic_map_chr$i.txt --numThreads=160 
done
```
The genetic map file was obtained from the website https://alkesgroup.broadinstitute.org/Eagle/ . The resulting phased output is then used to run the fineSTRUCTURE pipeline with the following code
```
for i in {1..22}; do plink --vcf all.merge.chr$i.phased.vcf.gz --recode12 --double-id --out FS_chr$i;done
for i in {1..22}; do plink2chromopainter.pl -p=FS_chr$i.ped -m=FS_chr$i.map -o=FS_chr$i.phasefile;done
for i in {1..22}; do convertrecfile.pl -M hap FS_chr$i.phasefile ../GENETIC MAP/chr$i.genetic.map FS_chr$i.recombfile;done
fs FS.cp -v -import FS.settings -go
```

# The data can be found here:

[APR](https://www.mbru.ac.ae/the-arab-pangenome-reference/)



