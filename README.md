# astreaimpute2

(c) Astrea Forensics, Ed Green - 2024

v14

Multithreaded software for generating genotype files or VCFs from mpileup data and haplotype records.

```
OPTIONS:

-H [required] Input file of haplotypes
-m [required] Input file of mpileup data with base observations
-M [optional] maximum coverage cutoff;            default = 24
   Ignores data at sites with observed coverage higher
   than this value, if specified.
-d [optional] downsample coverge to this;         default = 32
   At sites whose observed coverage is under the cutoff
   but higher than this number, down-sample the observed
   data to this amount.
-p [optional] percentile coverage cutoff;         default = none
   If specified, calculate the empirical coverage of the
   input mpileup data. Then, use all sites whose observed
   coverage is less than this percentile. For example:
   -p 95 would use all sites that are <= the observed
   coverage at 95 percent of all sites in the mpileup data
   *Note that this option overrides the value of -M   
-e [optional] rate of sequencing error;           default = 0.00500
-E [optional] rate of C->T sequencing error;      default = 0.00600
   *Note that -E must be >= -e   
-b set to true if sample is from male; default = 1
-g [optional] genome version (19 or 38); default=19
-V [optional] make VCF output
-I [optional] sampleID to assign for VCF file
-B [optional] Bayes Factor ratio cutoff           default = 100
-G [optional] Genotype likelihood ratio cutoff    default = 100
Sites with Bayes Factor and Genotype likelihood ratios
better than the cutoff will be included in the output
VCF file with PASS in the FILTER field.
```

###OUTPUT: 
Generates  tab-delimited output of this format:
`SNP_ID  chr     pos     Allele1 Allele2 call_LR Alt_fr  BayesFac`

where SNP_ID, chr, and position are read from the input haplotype table file;
Allele1 and Allele2 are the genotype call determined by analysis of the input
mpileup data;
and call_LR is the likelihood ratio of the data from the most likely genotype
call versus the next most-likely genotype call, the Alt_fr is the observed
alternate (non-reference) allele-frequency from the input haplotype table,
and BayesFac is the Bayes Factor of the most-likely genotype call versus the
next most-likely genotype call.

Alternatively, VCF 4.1 output if called with the -V option. In this case,
genotype calls passing the Genotype likelihood cutoff and Bayes Factor cutoff
specified by the -G and -B options will be included in the output with `PASS`
in the filter field.

### Recommended invocation for hair data
```
${ASTREA}/astrea-impute2.14 
   -H ${HT}/GRC38-v2-chr${CHR}-ht.txt.gz 
   -m mpileup/${SAMP}.${CHR}.mp.txt.gz 
   -g 38 -M 10 -I ${SAMP} -b ${BOY} -V 
   -B 20 -G 1 
   | bgzip > vcfs/${SAMP}.hg38.${CHR}.ai2.B20G1.vcf.gz
```

### Description of haplotype table format
astrea-imput2.14 attempts to call genotypes at target sites listed in the
input haplotype table. This table contains one or more records for each
target site. These records describe the target site and linked sites and
the haplotypes that contain them. A simple example of a haplotype table
record is shown here:
```
## TARGET 1000 A C rsTEST1
1  900    G    T rsTEST1A1
1  1000   A    C rsTEST1
1  1100   G    T rsTEST1B1
## Haplotypes
## _*_
2  000
2  111
## E_R
```
Each record has two sections: The Sites section, beginning with:
```## TARGET POS REF ALT ID```
and a haplotypes section, beginning with:
```## Haplotypes```

The Sites section contains one line per site. Each site line
should have a tab-delimited list of

`CHR	   POS	 REF  ALT  ID`.
The Target site should be one of these lines.

The Haplotypes section is a condensed description of the haplotypes
as described by the sites listed in the Sites section. It contains
one line per haplotype. Each line is of the form:

`N  BIT-STRING`

where N is the number of observed haplotypes of this type and
BIT-STRING is a description of the allelic configuration of the
sites listed in the Sites section where 0 is the reference allele
and 1 is the alternate allele.

Preceding the haplotypes line is a line of the form:

`##	  _*_`

This line marks where the Target site is in each haplotype with an
asterisk (*) and where linked sites are with underscore (_). The
purpose of this line is as a convenient visual indication of where
the target site lies amongst the linked sites within the haplotypes.

`test_hap` is a validation program that can be used to test for
properly formatted Haplotype files.

### Description of input mpileup format
astrea-impute2.14 requires mpileup format input data as generated
by `samtools Version 1.10` using this syntax:
`samtools mpileup -q 20 -Q 20 -s -a -f ${HG} -l ${SITE_POSITIONS} ${BAM_FILE}`
where `${HG}` is the reference genome, `${SITE_POSITIONS}` is a file
listing the genome positions at which to make mpileup output,
i.e., all the positions listed in any Haplotype record, and
`${BAM_FILE}` is the filename of the input bam file with the
aligned input sequence data.

The `-q` and `-Q` arguments can be set to appropriate map quality
and base quality cutoffs.


`test-pileup` is a validation program that can be used to test for
properly formated mpileup files.