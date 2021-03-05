# cgap-scripts

#### AddReadGroups.py
Add read groups to a BAM file using pysam and samtools to multi-thread compression.
Can handle multiple read groups.

#### duplines_vcf.py
Removes duplicate variants from sorted input VCF file while keeping the most informative one.

#### portal_reformat_vcf.py
Re-formats a VCF file by adding custom tags and fields used for portal ingestion:

  - expands VEP annotations by appending a bool (0, 1) to each transcript to flag the most severe one
  - adds a new tag with maximum DS for SpliceAI if any, `spliceaiMaxds`
  - adds a summary for the most severe transcript / gene using the tag `GENES`
  - adds the type of variant in `variantClass`

`GENES` tag definition

    ##INFO=<ID=GENES,Number=.,Type=String,Description=". Subembedded:'genes':Format:'most_severe_gene|most_severe_transcript|most_severe_feature_ncbi|most_severe_hgvsc|most_severe_hgvsp|most_severe_amino_acids|most_severe_sift_score|most_severe_polyphen_score|most_severe_maxentscan_diff|most_severe_consequence'">

#### samplegeno.py
Adds genotypes information for each variant to info field using the tag `SAMPLEGENO`.

Tag definition

    ##INFO=<ID=SAMPLEGENO,Number=.,Type=String,Description="Sample genotype information. Subembedded:'samplegeno':Format:'NUMGT|GT|AD|SAMPLEID'">

Example

    SAMPLEGENO=0/1|C/T|9/4|NA12877,0/0|C/C|34/0|NA12878,0/1|C/T|12/17|NA12879

#### add_altcounts_by_gene.py
Add alternate allele counts on the most severe gene per sample.
Requires `SAMPLEGENO`.
The count is appended to `SAMPLEGENO` tag.

Tag expansion

    ##INFO=<ID=SAMPLEGENO,Number=.,Type=String,Description="Sample genotype information. Subembedded:'samplegeno':Format:'NUMGT|GT|AD|SAMPLEID|AC'">
