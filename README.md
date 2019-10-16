# sv-tools
Scripts and software for filtering, visualizing and analysing genomic structural variants

## get\_readpair\_type\_summaries\_by\_sv\_no\_parallel.R
A script to quantify the types of read pair orientations that are present near the endpoints of structural variants (SVs).  For each SV in an input file, the script will count the number and types of readpair orientations that occur within some specified distance of the SV endpoints in individuals with at least one SV allele (0/1 or 1/1 genotypes).  The mean of these values is than calculated for each readpair type and the results are reported.  The script is useful in identifying complex SV false positives that occur when signals for multiple SV types (DEL, DUP, INS) are present together.

Readpair types that are considered are:
- *both\_negative*: Both reads in the pair are aligned to the negative DNA strand (reverse orientation).  Can indicate the presence of an inversion.
- *both\_positive*: Both reads in the pair are aligned to the positive DNA strand.  Can indicate the presence of an inversion.
- *outward\_facing*: The first read in the pair is aligned to the negative strand while the second is aligned to the positive strand.  Can indicate a duplication.
- *long*: The two reads are in a standard orientation however the distance between them is outside the 99.9th percentile of all mapped read pairs in the sample BAM file.
- *normal*: The two reads are in a standard orientation and are within the 99.9th percentile of all mapped read pairs in the sample BAM file.
