# diplocheck
This repository provides analysis details for the following publication:

Nürnberger, Baird, Čížková, Bryjová, Mudd, Blaxter and Szymura. 2021. A dense linkage map for a large repetitive genome: Discovery of the sex-determining region in hybridising fire-bellied toads (*Bombina bombina* and *B. variegata*)

The linkage map was constructed in a F2 population of a *B. bombina* x *B. variegata* intercross. Diplotype expectations for the two F0 (grandparental, GP) individuals are *B. bombina* homozygote (BbHOM) and *B. variegata* homozygote (BvHOM), respectively, and F1 individuals are expected to be HET. All loci for which the inferred diplotypes for these individuals deviated for this expectation, were rescored based on manually selected variants as a check on the inference scheme. This repository provides the Perl scripts and data to rescore one example locus (locus 5568).

1. **Visualising the raw read data**

Raw read counts (from bam files) are summarised in `*.sync` files produced by PoPoolation2 (Kofler et al. 2011) mpileup2sync (see publication for details). Starting at column 4, each column summarises the data for one sample. Of interest here are the F0 and F1 individuals. The counts of reads supporting each of six sequence states are separated by colons (A:T:C:G:N:D), where N is IUPAC code for ‘any base’ and D is a deletion. 

locus |	position | reference_state |	BvGP |	BbGP | F1F6 | F1M |	F1F7
----- | -------- | --------------- |  ---- |  ---- | ---- | --- | ----
Bv_contig5568_asmbl_rev | 1 | C | 0:0:2:0:0:0 | 0:0:16:0:0:0 | 0:0:14:0:0:0 | 0:0:6:0:0:0 | 0:0:9:0:0:0
Bv_contig5568_asmbl_rev | 2 | T | 0:2:0:0:0:0 | 0:16:0:0:0:0 | 0:14:0:0:0:0 | 0:6:0:0:0:0 | 0:9:0:0:0:0
... | | | | | | |

Script `raw_read_cover_from_sync_GP_F1.pl` reformats these data so that they can be visualised with R/ggplot2. Input files: `gp_f1_locus5568.sync`, `locus_list.txt`, `gp_f1_inds_annotated.txt` Output file: `gp_f1_for_plotting.txt`. The output file is imported into R and plotted as follows:


```
df$state_f = factor(df$state, levels = c("R","A","C","G","T","D"))
df$symbolsize = as.factor(df$symbolsize)
library(ggplot2)
ggplot (df[which(df$locus == 5568),],aes(x=position,y=coverage,colour=state_f)) + geom_point(aes(size=symbolsize)) +
scale_size_manual(values=c(0.5,2.5)) + scale_colour_manual(values = c("#666666","#FF0033","#0033FF","#33CCFF","#FF9900","#CC33CC")) +
theme_bw(base_size = 18) + guides (size=FALSE,colour=guide_legend(title = 'state')) + facet_wrap (~ sample, ncol = 2)
```

This produces the following plot (analogous to Figure 8 in the publication). 

2. **Selecting variants**

The plot shows a homozygous variant (`A`) in the *B. bombina* GP with coverage just over 300. As expected, it is found in heterozygous state in each of the F1s. It can thus be used to identify BbHOM, HET and BvHOM diplotypes in all samples. 

Script `raw_read_cover_from_sync_GP_comparison.pl` parses the same `*.sync` file and produces a tab-delimited file of variant states in which this selection can be recorded. Input files: as above. Output file: `gp_variants.txt`. A subset of that file is shown here:

locus | position | sample | topdaown_order | state | coverage | type
----- | -------- | ------ | -------------- | ----- | -------- | ----
5568 | 187 | BvGP | 2 | G | 21 | 
5568 | 217 | BvGP | 2 | C | 27 | 
5568 | 241 | BbGP | 1 | A | 301 | 
5568 | 247 | BvGP | 2 | T | 30 | 
5568 | 248 | BvGP | 2 | G | 30 | 

topdown order = 1: the most highly covered state for a given sample and position, topdown order = 2: the second most highly covered state for a given sample and position. 

The `A` variant is at position 241 and is selected by typing `Bb` in the `type` column of that line.

3. **Rescoring samples**

For the PoPoolation2 mpileup2sync analysis, the dataset was subdivided into sets of up to 10 samples. The samples of a given set are rescored with `diplotyping_nonstandard_loci.pl`, with a set name supplied on the command line. The script then analyses `<setname>.sync` using the annotated variants stored in `gp_variants.txt`. For F2 individuals, a key file (`<setname>.key`) is also needed that identifies the samples in `<setname>.sync`. 

Because highly covered, unambiguous variants are manually selected, diplotyping is straightforward. If the two most highly covered sequence states are supported by more than `$MINCOVER` (= 5) reads, a heterozygous diplotype is inferred. If only one sequence state exceeds `$MINCOVER` reads, the sample is deemed homozygous. 

Here, the process is illustrated with set `gp_f1_locus5568`. The reference state at position 241 is `C`. Therefore AA = BbHOM, AC = HET and CC = BvHOM. The rescored diplotypes are written to`gp_f1_locus5568_rescored.tx`. The last column lists for a given individual the number of annotated positions that could not be scored, because the maximum coverage across all sequence states was smaller than or equal to `$MINCOVER`.

locus | sample | diplotype | missing
----- | ------ | --------- | -------
5568 | BvGP | v/v | 0
5568 | BbGP | b/b | 0
5568 | F1F6 | b/v | 0
5568 | F1M | b/v | 0
5568 | F1F7 | b/v | 0

Multiple variants at a given locus may allow additional haplotypes to be identified (up to four in the *B. bombina* x *B. variegata* cross). In fact, locus 5568 is a case in point. These variant positions can also be annotated (*e.g.* with b1 or b2 for different *B. bombina* haplotypes) in `gp_variants.txt`. `diplotyping_nonstandard_loci.pl` keeps track of these and reports diplotypes such as b2/v1 in that case.

4. **Converting diplotypes to Lep-MAP3 format**

`recode_diplotypes_for_lepmap.pl` reads in the output of `diplotyping_nonstandard_loci.pl` and coverts each diplotype into Lep-MAP3 format. These diplotype codes are then written to the input file for Lep-MAP3.

locus | sample | diplotype | missing | Lep-MAP3 coding
----- | ------ | --------- | ------- | ---------------
5568 | BvGP | v/v | 0 | 0 0 0 0 1 0 0 0 0 0
5568 | BbGP | b/b | 0 |  1 0 0 0 0 0 0 0 0 0
5568 | F1F6 | b/v | 0 | 0 1 0 0 0 0 0 0 0 0
5568 | F1M | b/v | 0 | 0 1 0 0 0 0 0 0 0 0 
5568 | F1F7 | b/v | 0 | 0 1 0 0 0 0 0 0 0 0 
