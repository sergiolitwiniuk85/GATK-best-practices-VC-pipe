
# Mapping & Pre-processing
## GATK best practices
### Notes:

* Bwa-mem mapping tool for genomics data.
* Star v2 for RNA-Seq data.


(How to) Generate an unmapped BAM from FASTQ or aligned BAM : [Convert FastQ to uBAM](https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM)

(How to) Install all software packages required to follow the GATK Best Practices:[ Install software packages guide ](https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices)

ReadGroup important notes: [Meaning of the read group](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)

Tag @RG must be present for GATK analysis: [@RG error](https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information)


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




    subgraph Indel-based:Realignment
        Mark---rtc(GenomeAnalysisTK.jar:\nRealignerTargetCreator\n-Find intervals\nRoche 454 not compatible)

        rtc ---ir(IndelRealigner)-----r_out(Realigned\nBAM)
    end



    subgraph Base Quality Score Recalibration
        r_out --> r(BQSR method identifies bias and applies correction:\nGenomeAnalysisTK.jar:\nBaseRecalibrator\nPrintReads\nAnalyzeCovariates Plots)
        r--->rc(Recalibrated\nBAM\n-ready for variant calling analysis)
    end

%%Comment: GATK recommends using uBAM instead of fastQ as it can support metadata (RG,)






```