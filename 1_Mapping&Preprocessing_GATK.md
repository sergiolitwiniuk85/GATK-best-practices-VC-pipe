
# Mapping & Pre-processing
## GATK best practices

This file provides an introduction to the tools and the good practices workflow used for mapping and preprocessing in gatkv4.0.

### Notes:

* Bwa-mem mapping tool for genomics data.
* Star v2 for RNA-Seq data.
* uBAM recommended over FastQ

(How to) Generate an unmapped BAM from FASTQ or aligned BAM : [Convert FastQ to uBAM](https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM)

(How to) Install all software packages required to follow the GATK Best Practices:[ Install software packages guide ](https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices)

Software packages:

    BWA
    SAMtools
    Picard
    Genome Analysis Toolkit (GATK)
    IGV
    RStudio IDE and R libraries ggplot2 and gsalib



ReadGroup important notes: [Meaning of the read group](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)

Tag @RG must be present for GATK analysis: [@RG error](https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information)

<br />

**Used tools:**

**[bwa-mem:](https://github.com/lh3/bwa)**
BWA is a software package for mapping DNA sequences against a large reference genome, such as the human genome. BWA-MEM is recommended for query sequences longer than ~70bp for a variety of error rates (or sequence divergence).

**[SamToFastq (Picard):](https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq)**
Converts a SAM or BAM file to FASTQ. This tool extracts read sequences and base quality scores from the input SAM/BAM file and outputs them in FASTQ format. This can be used by way of a pipe to run BWA MEM on unmapped BAM (uBAM) files efficiently.

Usage example:

```
java -jar picard.jar SamToFastq \
     I=input.bam \
     FASTQ=output.fastq
```
<br />


**[MergeBamAlignment (Picard):](https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment)**
Merge alignment data from a SAM or BAM with data in an unmapped BAM file. This tool produces a new SAM or BAM file that includes all aligned and unaligned reads and also carries forward additional read attributes from the unmapped BAM (attributes that are otherwise lost in the process of alignment). The purpose of this tool is to use information from the unmapped BAM to fix up aligner output. The resulting file will be valid for use by other Picard tools. For simple BAM file merges, use MergeSamFiles. Note that MergeBamAlignment expects to find a sequence dictionary in the same directory as REFERENCE_SEQUENCE and expects it to have the same base name as the reference FASTA except with the extension ".dict". If the output sort order is not coordinate, then reads that are clipped due to adapters or overlapping will not contain the NM, MD, or UQ tags.

Usage example:
```
java -jar picard.jar MergeBamAlignment \
       ALIGNED=aligned.bam \ 
       UNMAPPED=unmapped.bam \ 
       O=merge_alignments.bam \
       R=reference_sequence.fasta
```
<br />

**[MarkDuplicates (Picard):](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)**
Identifies duplicate reads.

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for additional notes on PCR duplication artifacts. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates.

Usage example:
```
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt
```
<br />

**[BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/21905050792603-BaseRecalibrator)**
First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

Usage example:
 ```
 gatk BaseRecalibrator \
   -I my_reads.bam \
   -R reference.fasta \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table
   ```

<br>

**[ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/21905038144155-ApplyBQSR)**
This tool performs the second pass in a two-stage process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the base qualities of the input reads based on the recalibration table produced by the BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

Usage example:
 ```
 gatk ApplyBQSR \
   -R reference.fasta \
   -I input.bam \
   --bqsr-recal-file recalibration.table \
   -O output.bam
 ```

<br>


### Optional:
**If using FastQ as input:**

**[CleanSam (Picard):](https://broadinstitute.github.io/picard/command-line-overview.html#CleanSam)**
Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads

<br />

**[FixMateInformation (Picard):]()**
Verify mate-pair information between mates and fix if needed.This tool ensures that all mate-pair information is in sync between each read and its mate pair. If no OUTPUT file is supplied then the output is written to a temporary file and then copied over the INPUT file. Reads marked with the secondary alignment flag are written to the output file unchanged.

Usage example:
```
java -jar picard.jar FixMateInformation \
       I=input.bam \ 
       O=fixed_mate.bam
```
<br>

**[AddOrReplaceReadGroups (Picard)]()**
Replace read groups in a BAM file.This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.

For more information about read groups, see the GATK Dictionary entry.

This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH) see [GA4GH](http://ga4gh.org/#/documentation).

Usage example:
```
java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
      O=output.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
```
<br>


**GnomAD v4:**

The major steps taken in v4 sample QC are listed here (with details following):

* Samples with obviously bad quality were removed using hard filters:Applied the novel metric Contamination from Homozygous Alternate Reference Reads [CHARR](https://www.biorxiv.org/content/10.1101/2023.06.28.545801v1)

QC metrics were computed using the Hail sample QC module on all autosomal, biallelic SNVs. In summary, the hard filtering step removed samples that were outliers for the following metrics:

Contamination (mean AB biallelic SNPs) > 1.5%<br>
Number of singletons > 5000<br>
Heterozygous: homozygous ratio > 10<br>
Bases with DP over 1X < 5x107<br>
Bases with DP over 20X < 4x107<br>
Proportion chimera > 0.05<br>
Mean chromosome 20 DP < 10X<br>
Callrate calculated over high quality sites < 0.8<br>
* Inferred exome capture platform for each sample

Performed a principal components analysis (PCA) using missingness rate calculated per sample for each exome calling interval (intersection of Broad and UK Biobank intervals). Using nine principal components (PCs), we inferred that there are 20 unique exome sequencing platforms in gnomAD v4.
* Inferred sex karyotypes for each sample using normalized sex chromosome coverage.

Used [Hail’s impute_sex_chromosome_ploidy](https://hail.is/docs/0.2/vds/hail.vds.impute_sex_chromosome_ploidy.html#hail-vds-impute-sex-chromosome-ploidy) impute_sex_chromosome_ploidy module to calculate mean coverage across chromosomes X (chrX) and Y (chrY).

* Inferred relatedness between samples, including between the v4 exomes and genomes, and identified pairs of first and second degree relatives.

CUDA-based implementation of KING ([cuKING](https://github.com/populationgenomics/cuKING/tree/1d29d72a184da33cf2500f8b537ae56ce66a7e96); implemented by Leonhard Gruenschloss)

* Used genetic similarity to cluster and identify genetic ancestry groups
* Filtered samples based on QC metrics

More details on GnomAD v4.0 [here](https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)

<br>

## Mapping workflow
*Indel Realignment step is deprecated from v4.

```mermaid


flowchart TD;

    subgraph Input:Unmapped-Recommended
    ref----> bam_output
    UM(Unmapped_BAM\nuBAM)----pi(Picard:SamToFastq)--FastQ---bam_output(bwa-mem)--SAM----pik(Picard:MergeBamAlignment\n-Mapped\n-Cleaned\n-Sorted_SAM)
    
    UM --Unmapped\nSAM---- pik
    end



    subgraph Input:FastQ-Alternative
    ref(Ref_Genome GRCh37/38/T2T)
    fastqc(FastQC)-----bwa
        ref-.-> bwa(bwa-mem\nflags:-M -R)--->picard(Optional Picard\n-CleanSam & FixMateInformaMon\nto add readgroup \n as metadata\n-AddOrReplaceReadGroups)
        
           
        picard---sam_compress(SAM\n-Mapped\n-Cleaned\n-Sorted)
    end


    sam_compress -.- BAM
    pik ---> BAM



    subgraph MarkDuplicates
    BAM ---> Mark(MarkDuplicates:Picard)
    end





%%indelRealignment step deprecated


    subgraph Base Quality Score Recalibration
        Mark --BAM--> r(BQSR - GenomeAnalysisTK.jar:\nBaseRecalibrator\nApplyBQSR)
        r--->rc(Recalibrated BAM\n-ready for variant calling analysis)
    end

%%Comment: GATK recommends using uBAM instead of fastQ as it can support metadata (RG,)






```