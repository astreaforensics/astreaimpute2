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

OUTPUT: 

VCF or genotype file with format:

SNP_ID	chr	pos	Allele1	Allele2


### Recommended invocation for hair data
${ASTREA}/astrea-impute2.14 \
   -H ${HT}/GRC38-v2-chr${CHR}-ht.txt.gz \
   -m mpileup/${SAMP}.${CHR}.mp.txt.gz \
   -g 38 -M 10 -I ${SAMP} -b ${BOY} -V \
   -B 20 -G 1 \
   | bgzip > vcfs/${SAMP}.hg38.${CHR}.ai2.B20G1.vcf.gz
