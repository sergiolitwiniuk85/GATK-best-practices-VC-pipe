## Using GATK best practices for VC

### 1. Mapping and Preprocessing
Following GATK best practices, the mapping and preprocessing steps are crucial for accurate variant discovery. The process begins with aligning the sequencing reads to the reference genome using BWA-MEM. After alignment, the SAM files are converted to BAM files and sorted using Samtools. Next, duplicates are marked with GATK MarkDuplicates to ensure that duplicate reads do not interfere with downstream analyses. Finally, Base Quality Score Recalibration (BQSR) is performed using GATK's BaseRecalibrator and ApplyBQSR tools to adjust the quality scores of the reads based on known variants, enhancing the accuracy of the subsequent variant calling.

https://github.com/sergiolitwiniuk85/GATK-best-practices-VC-pipe/blob/main/1_Mapping%26Preprocessing_GATK.md

### 2. Variant Calling
Adhering to GATK best practices, the variant calling process involves several key steps to identify SNPs and indels in the preprocessed BAM files. The GATK HaplotypeCaller is used to perform local realignment of reads and call variants, producing a GVCF file. This file is then processed with GATK GenotypeGVCFs to consolidate GVCF files and call genotypes across samples. The final step involves applying hard filters to the called variants using GATK VariantFiltration, ensuring that only high-confidence variants are retained for further analysis. This comprehensive approach ensures robust and reliable variant detection, critical for genomic studies.
