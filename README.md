# cgap-scripts

#### AddReadGroups.py
Adds read groups to a BAM file using pysam and samtools to multi-thread compression.
Can handle multiple read groups.

#### duplines_vcf.py
Removes duplicate variants from sorted input VCF file while keeping the most informative.

#### portal_reformat_vcf.py
Re-formats a VCF file by adding custom tags and fields used for portal ingestion:

  - expands VEP annotations by appending a bool (0, 1) to each transcript to flag the most severe
  - adds a new tag with maximum DS for SpliceAI if any, `spliceaiMaxds`
  - adds a summary for the most severe transcript / gene using the tag `GENES`
  - adds the type of variant in `variantClass`

`GENES` tag definition

    ##INFO=<ID=GENES,Number=.,Type=String,Description=". Subembedded:'genes':Format:'most_severe_gene|most_severe_transcript|most_severe_feature_ncbi|most_severe_hgvsc|most_severe_hgvsp|most_severe_amino_acids|most_severe_sift_score|most_severe_polyphen_score|most_severe_maxentscan_diff|most_severe_consequence'">

#### samplegeno.py
Adds genotypes information for each variant to INFO field using the tag `SAMPLEGENO`.

Tag definition

    ##INFO=<ID=SAMPLEGENO,Number=.,Type=String,Description="Sample genotype information. Subembedded:'samplegeno':Format:'NUMGT|GT|AD|SAMPLEID'">

Example

    SAMPLEGENO=0/1|C/T|9/4|NA12877,0/0|C/C|34/0|NA12878,0/1|C/T|12/17|NA12879

#### add_altcounts_by_gene.py
Adds alternate alleles count on the most severe gene per sample.
Requires `SAMPLEGENO`.
The count is appended to `SAMPLEGENO` tag.

Tag expansion

    ##INFO=<ID=SAMPLEGENO,Number=.,Type=String,Description="Sample genotype information. Subembedded:'samplegeno':Format:'NUMGT|GT|AD|SAMPLEID|AC'">

#### depth_filter.py
Filters variants called in VCF format based on DP (depth of coverage) from the genotype column.

#### split_bam.py
Splits a bam file in multiple chunks. Each chunk is a bam file that contains the maximum number of reads speficied. The original header is used in writing each bam. Chunks maintain the order of the reads and are ordered with an index in the filename.

#### compare_bams.py
Compares reads in two bam files line by line and print mismatching lines.
