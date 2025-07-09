# GATK Best Practices Variant Calling Pipeline

## Overview
This pipeline implements GATK's best practices for accurate variant discovery, covering read mapping, preprocessing, and variant calling. The workflow ensures reliable identification of SNPs and indels for genomic studies.

![Variant Calling Workflow](https://via.placeholder.com/800x200?text=GATK+Variant+Calling+Workflow)

## Pipeline Steps

### 1. Mapping and Preprocessing
![Mapping Diagram](https://via.placeholder.com/300x150?text=Mapping+%26+Preprocessing)

Following GATK best practices:
1. **Alignment**: Map sequencing reads to reference genome using BWA-MEM
2. **Format Conversion**: Convert SAM → sorted BAM with Samtools
3. **Duplicate Marking**: Identify PCR duplicates with GATK MarkDuplicates
4. **BQSR**: Adjust base quality scores using BaseRecalibrator and ApplyBQSR

[Detailed Documentation](https://github.com/sergiolitwiniuk85/GATK-best-practices-VC-pipe/blob/main/1_Mapping%26Preprocessing_GATK.md)

---

### 2. Variant Calling
![Variant Calling](https://via.placeholder.com/300x150?text=Variant+Calling)

Accurate variant detection:
1. **HaplotypeCaller**: Perform local realignment and call variants (outputs GVCF)
2. **GenotypeGVCFs**: Consolidate GVCF files and call genotypes
3. **Variant Filtration**: Apply hard filters to retain high-confidence variants

[Detailed Documentation](https://github.com/sergiolitwiniuk85/GATK-best-practices-VC-pipe/blob/main/2_VariantCalling_GATK.md)

---

## Key Features
- Implements GATK's gold-standard workflow
- Optimized for accuracy and reproducibility
- Comprehensive quality control at each step
- Compatible with WGS and targeted sequencing data

## Dependencies
- [GATK](https://gatk.broadinstitute.org/) (v4.2+)
- [BWA](https://bio-bwa.sourceforge.net/) (v0.7+)
- [Samtools](https://www.htslib.org/) (v1.12+)
- [Nextflow](https://www.nextflow.io/) (v22.04+)

## References
1. [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices)
2. Van der Auwera, G.A. et al. (2013). From FastQ data to high-confidence variant calls. *Nature Protocols* 8, 2483-2512.
