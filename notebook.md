# Investigating the capacity of the black soldier fly *Hermetia illucens* (Diptera: Stratiomyidae) to minimize pathogens in livestock waste

## Background 
Heightened demands for animal products worldwide has resulted in billions of tons of livestock waste per year1. Excessive animal waste has devastating impacts on the environment, including soil contamination, rising greenhouse gas emissions, and pathogen transmission. To mitigate these issues, the concept of a circular bioeconomy has been proposed, which advocates for the use of renewable biological resources to reduce waste and promote environmental sustainability. The emerging sector of insect farming has the potential to enhance this initiative through specific insects, namely *Hermetia illucens*, the black soldier fly (BSF). While BSF is primarily utilized as feed for livestock and poultry, it also exhibits remarkable waste conversion capabilities. This insect can consume virtually any biowaste, reduce waste dry mass by up to 50%, and efficiently convert waste to biomass, all of which contribute to a significant reduction in greenhouse gas emissions compared to meat production. Although the use of BSF for agricultural sustainability is promising, numerous crucial aspects of its biology remain unknown, including its microbial composition before and after waste consumption, as well as the potential pathogens they may spread into the food chain. 

**Objective**: Using metagenomics approaches, we aim to analyze changes in microbial community composition in dairy manure before and after black soldier fly digestion.

## Data Overview
**_Experimental Design_**

We had a total of three sampling time points (T0, T3, and T7) and two treatment types. The first treatment type, the control (DM_C), had three replicates and contained **only** 1 kg of dairy manure. The next treatment type was dairy manure **with** BSF (DM_BSF), which had four replicates and contained 1000 BSF larvae along with 1 kg of dairy manure. At each sampling timepoint (T0, T3, and T7), three swabs of dairy manure and ~10 BSF larvae were sampled for DNA extraction. 

**_DNA Extraction_**

After completing sample collection, we extracted DNA from **all** fecal and larval samples using the Qiagen PowerFecal Pro DNA kit, following the manufacturer's standard protocol. DNA concentration and purity were subsequently assessed using Qubit and Nanodrop measurements, respectively.

**_Library Preparation and Sequencing_**

Due to funding constraints, we sequencing only nine fecal samples representing all three timepoints. Before submission to sequencing, we pooled 5uL from the three highest-quality samples per replicate - based on Nanodrop and Qubit assessments - into a single 1.5mL tube, yielding a total volume of 15uL of total sample. Library preparation and sequencing were conducted at the University of Tennessee Genomics Core usng half a lane on an Illumina NovaSeq S4 flow cell.

## Data Analyses


![BSF Metagenomics Pipeline](https://github.com/AllysonDekovich/BSF_Metagenomics/blob/main/figures/BSF%20Metagenomics%20Schematic%20.png)

The above schematic summarizes the metagenomics workflow to analyze the data. Detailed descriptions of each software tool, along with code snippets, are provided below.

**_FastQC_**

```
#!/bin/bash
#SBATCH --job-name=fastqc_pretrimmed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu

for file in *fastq.gz
do
	fastqc -t 6 -o QC_output_pretrim $file
done
```

Initial assessment was pretty good! Reads were of high quality and there were no red flags, other than removing adapters.

**_Fastp_**

Adapters and poor quality reads were trimmed with `fastp`

```
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-9

# Only select lines with _R1_001.fastq.gz and pick the correct one by SLURM_ARRAY_TASK_ID
infile1=$(grep "_R1_001.fastq.gz" filenames.txt | sed -n "${SLURM_ARRAY_TASK_ID}p")
#echo "R1: $infile1"

# Replace _R1_001.fastq.gz with _R2_001.fastq.gz
infile2=$(basename "$infile1" | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
#echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_001.fastq.gz//')
#echo "outfile: $outfile"

# run fastp
fastp \
--in1 "$infile1" \
--in2 "$infile2" \
--out1 ./output/"${outfile}_R1_trimmed.fastq.gz" \
--out2 ./output/"${outfile}_R2_trimmed.fastq.gz" \
--adapter_sequence TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
--adapter_sequence_r2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
--detect_adapter_for_pe \
--qualified_quality_phred 25 \
--html ./output/"${outfile}".html \
--thread 4
```

The `multiqc` report revealed the "nextera transposase sequence" was the most prominent adapter, so I incorporated the forward and reverse sequence using the `--adapter-sequence` flag. I also turned on the `--detect_adapter_for_pe` setting to get any sequences I may have missed, as `fastp` has the capability to automatically detect common adapters. Post-trim results look great!

**_Nonpareil_**

Ultimately, I would like to assemble the metagenomic reads into MAGs to evaluate species composition and diversity comparisons across timepoints, however, you need pretty high coverage to make assembling MAGs even worth it. To check this, I am running a program called _Nonpareil_, which estimates how well our sequencing effort was able to capture the DNA in microbial samples. It first performs a redundancy estimation, where it starts by selecting a subset of sequencing reads and performs ungapped alignments of these reads against the whole dataset to identify redundancy (i.e., how often the same or similar sequences show up). It then applies statistical testing (Turing-Good) to infer how many new, unseen sequences you would get if you sequenced more, which helps you estimate how complete your sampling is (i.e., how close you are to sequencing everything in the sample).

_Nonpareil_ utilizes **abundance-weighted coverage**, which better reflects the reality of shotgun metagenomics. In contrast to species-based coverage - where each species is treated equally and rare species can drastically lower your estimated coverage - abundance-weighted coverage emphasizes how much of the **total DNA** (thus, the dominant organisms) is represented. Since rare species are often missed in environmental samples unless sequencing is extremely deep, species-based methods can underestimate coverage and complicate tasks such as MAG assembly. Abundance-weighted coverage helps overcome this by focusing on the most abundant organisms, which are typically of greater ecological importance and easier to recover with less sequencing effort.

Prior to running the program,  `fastq` files need to be converted to `fasta` format, which I did using BioPython:

```
#!/usr/bin/python3
from Bio import SeqIO
import gzip
from io import StringIO
import glob
import os

path_to_files = "/lustre/isaac24/proj/UTK0032/adekovic/BSF_metagenomics/04_nonpareil/*.fastq.gz"

for file in glob.glob(path_to_files):
	filename = os.path.basename(file)

	if filename.endswith(".fastq.gz"):
		prefix = filename[:-9]
		output_file = os.path.join(os.path.dirname(file), prefix + ".fasta")

		with gzip.open(file, "rt") as in_handle, open(output_file, "w") as out_handle:
			count = SeqIO.convert(in_handle, "fastq", out_handle, "fasta")
			print("Converted " + str(count) + " records from " + filename + " to " + prefix + ".fasta")
```
**Note**: only R1 reads were used in this analysis because the program assumes independence of events between sequencing reads.

SLURM script to utilize HPC resources:

```
#!/bin/bash
#SBATCH --job-name=fasta_conversion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu



python3 /lustre/isaac24/proj/UTK0032/adekovic/BSF_metagenomics/04_nonpareil/fasta_conversion.py
```

Perform the redundancy estimation:

```
#!/bin/bash
#SBATCH --job-name=nonpareil_redundancy_alignment
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-9

infile=$(sed -n -e "${SLURM_ARRAY_TASK_ID} p" filenames.txt)
#echo $infile

outfile=$(basename "$infile" | sed 's/.fasta//')
#echo $outfile

nonpareil -s "$infile" -T alignment -f fasta -b "$outfile"
```
*`-s`: path to the input file (in fasta format)<br/>
*`-T`: _nonpareil_ algorithm. I choose `alignment` as it is more accurate, however, it takes much longer to run/more computational resources<br/>
*`-f`: format of the input file<br/>
*`-b`: prefix for the output files

**Update: kmer estimation - alignment method was taking too long**

```
#!/bin/bash
#SBATCH --job-name=nonpareil_redundancy_kmer
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-9

infile=$(sed -n -e "${SLURM_ARRAY_TASK_ID} p" filenames.txt)
#echo $infile

outfile=$(basename "$infile" | sed 's/.fastq.gz//')
#echo $outfile

nonpareil -s "$infile" -T kmer -k 15 -f fastq -b "$outfile"
```

Nothing has changed, except for the `kmer` mode and fastq instead of fasta (the `kmer` mode uses alignment scores to help with estimations). Originally, I kept getting this error:
```
Nonpareil v3.5.5
[ 0.0] WARNING: The kmer kernel implements an error correction function only compatible with FastQ
[ 0.0] Reading DM_C_T0_r1_Owings_S13_L001_R1_trimmed.fasta
[ 0.0] Picking 10000 random sequences
[ 0.0] Counting Kmers
Fatal error:
Reads are required to have a minimum length of kmer size
[ 0.3] Fatal error: Reads are required to have a minimum length of kmer size
```
The default kmer size is 24bp, and it turns out I have some reads that were 15 bp, so the program was failing. I was able to get it working by specifying a different kmer size with `-k`.

**Creating and Interpreting the Nonpareil curves**

To plot the curves for samples on a single plot, I needed to create a text file containing the sample names, an alias for labeling, and a color hex code.

```
File	Name	Color
DM_C_T0_r1_Owings_S13_L001_R1_trimmed.npo	T0_r1	"#D27D2D"
DM_C_T0_r2_Owings_S13_L001_R1_trimmed.npo	T0_r2	"#7B3F00"
DM_C_T0_r3_Owings_S13_L001_R1_trimmed.npo	T0_r3	"#6E260E"
DM_T3_r1_Owings_S16_L001_R1_trimmed.npo		T3_r1	"#FF69B4"
DM_T3_r3_Owings_S17_L001_R1_trimmed.npo		T3_r3	"#FF00FF"
DM_T3_r4_Owings_S18_L001_R1_trimmed.npo		T3_r4	"#F8C8DC"
DM_T7_r1_Owings_S19_L001_R1_trimmed.npo		T7_r1	"#0047AB"
DM_T7_r2_Owings_S20_L001_R1_trimmed.npo		T7_r2	"#6495ED"
DM_T7_r4_Owings_S21_L001_R1_trimmed.npo		T7_r4	"#CCCCFF"
```
Run an `Rscript` to plot the curves:
```
library(Nonpareil)

# read in the sample file, which gives each sample a color for easy identification.
samples <- read.table('npo_samples.txt', sep='\t', header=TRUE, as.is=TRUE)

# Create a single plot that visualizes the metagenomic coverage for all samples at once.
np_curves <- Nonpareil.set(as.vector(samples$File), col=samples$Color, labels=samples$Name, plot.opts=list(plot.observed=FALSE))
```

**Coverage curves from `nonpareil`**
![Nonpareil output: coverage curves per sample](https://github.com/AllysonDekovich/BSF_Metagenomics/blob/main/figures/nonpareil_coverage_curve_DM_T0_T3_T7_pre_host_trim.png)

* The x-axis plots the log-transformed sequencing effort, which is the number of reads or base pairs sequenced.<br/>
* The y-axis plots the estimated coverage of the community diversity or complexity on a 0 to 1 scale.

Key interpretative elements of the plots:
* Curve shape: the steeper the intial rise, the lower the community diversity.<br>
* Saturation point: where the curve plateaus indicates the maximum achievable coverage for the sequencing effort.<br>
* Diversity arrows: indicate sequence diversity (`Nd`) values for each sample.<br>
	* higher `Nd` values (skewed right) indicate more diverse samples, which require more sequencing effort than less diverse (skewed left) samples.

According to this plot, it seems like all of my samples are good candidates for MAG assembly! All samples seem to plateau around the upper limit of the y-axis. **Our sequencing effort has captured about 95-100% of the abundance-weighted DNA in each community (90% coverage is recommended by standard benchmarking).** 

Additionally, all samples seem to be quite diverse and experience similar diversity levels -- the `Nd` values all cluster together and are skewed more toward the right, indicating deeper sequencing is needed to fully capture all of the diversity present in the samples.

So far so good! To double check that our sequencing effort wasn't inflated because of host contamination, I will also re-run this analysis on host trimmed reads.

**bbduk**

I used `bbduk` (a function in the `bbtools` package) to remove host contamination from the samples. Excess host DNA can distort estimates of microbial species abundance and decrease the proportion of microbial reads recovered, which in turn diminishes the sensitivity of microbial community detection and obscures the identification of rare taxa.

``bduk` uses k-mer matching, where it breaks both the sequencing reads and the reference into small sequences. Any reads that contain enough k-mers (see below) in common with the host reference are discarded. It is fast and memory efficient, as it does not require full all-to-all alignments.

_Bos taurus_ reference genome can be found [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/). **The code shown below is for a single timepoint (DM_C_T0) for clarity; however, I applied the same analysis across all timepoints and samples.**

```
#!/bin/bash
#SBATCH --job-name=bbduk_DM_C_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-3

# Only select lines with _R1_001.fastq.gz and pick the correct one by SLURM_ARRAY_TASK_ID
infile1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" DM_C_T0_filenames.txt)
echo "R1: $infile1"

# Replace _R1_001.fastq.gz with _R2_001.fastq.gz
infile2=$(basename "$infile1" | sed 's/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz/')
echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_trimmed.fastq.gz//')
echo "outfile: $outfile"

# run bbduk

java -ea -Xmx80g -Xms80g -cp /nfs/home/adekovic/mambaforge/envs/bbmap/opt/bbmap-39.33-0/current/ jgi.BBDuk \
in1=$infile1 \
in2=$infile2 \
out1=${outfile}_R1_host_trimmed.fastq.gz \
out2=${outfile}_R2_host_trimmed.fastq.gz \
ref=GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz \
k=25 \
hdist=0 \
ktrim=r \
outm1=bos_taurus_R1.fastq.gz \
outm2=bos_taurus_R2.fastq.gz \
tpe \
stats=DM_C_T0_alignment_stats.txt
```
* `k=25`: size of k-mers used to scan for matches between reads and host reference. I set mine to 25 (which I also believe is default).<br>
* `hdist=0`: Hamming distance allowed for k-mer matching is zero, meaning only exact matches are considered indicative of contamination.<br>
* `ktrim=r`: Right trim mode; once a contaminant k-mer is found, the read is trimmed from that k-mer to the end (right side) of the read.<br>
* `tpe`: Trim pairs equally; if one read pair is trimmed, the mate is trimmed to the same length to keep them in sync.



**_MEGAHIT_**

Since the samples have more than enough coverage, I will go ahead with MAG construction. First, I used `MEGAHIT`, a de novo assembler, to assemble metagenomic regions into contigs for each biological replicate (DM_C_T0, DM_T3, DM_T7). `MEGAHIT` builds multiple de Bruijn graphs using a range of k-mer sizes -- smaller k-mer sizes help fill gaps and cover low-depth regions, while larger k-mers help resolve repeats and improve contiguity. The final output is a single fasta file for each time point (see 'co-assembly' below) that contains **contigs**, or continuous sequences from overlapping reads.

Since we have three technical replicates for each timepoint, I can combine them all in a co-assembly. The benefits of a co-assembly are:

* **Improved Genome Recovery**: recovers a higher fraction of the genome more than a single assembly alone.<br>
* **Enhanced Continuity**: Produces less fragmented assemblies and longer contigs.<br>
* **Efficiency**: By pooling replicates into a single MAG assembly, the computational effort decreases and it is simpler to combine now than to wait until after MAG assembly.

**The code shown below is for a single timepoint (DM_C_T0) for clarity; however, I applied the same analysis across all timepoints and samples.**

```
#!/bin/bash
#SBATCH --job-name=megahit_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long-bigmem
#SBATCH --qos=long-bigmem
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


megahit \
-1 DM_C_T0_r1_Owings_S13_L001_R1_host_trimmed.fastq.gz,DM_C_T0_r2_Owings_S14_L001_R1_host_trimmed.fastq.gz,DM_C_T0_r3_Owings_S15_L001_R1_host_trimmed.fastq.gz \
-2 DM_C_T0_r1_Owings_S13_L001_R2_host_trimmed.fastq.gz,DM_C_T0_r2_Owings_S14_L001_R2_host_trimmed.fastq.gz,DM_C_T0_r3_Owings_S15_L001_R2_host_trimmed.fastq.gz \
-t 12 \
-m 260000000000 \
--k-list 21,29,39,49,59,69,79,89,99,109,119,129,139,149,159,169,179,189,199 \
--min-contig-len 1000 \
-o DM_C_T0
```
* `--k-list`: The multiple kmer sizes `MEGAHIT` uses iteratively to build the assembly graph.
* `--min-contig-length`: sets a minimum threshold for contig length, to avoid keeping short, unreliable contigs in the assembly.

**_Fairy_**

Prior to binning reads, I generate coverage lists with `fairy`. Coverage lists will be generated for each timepoint; each list will contain contigs (identified by `megahit`), contig length, average depth across all samples, and depth in each sample. Contigs from the same microbial genome tend to have similar coverage patterns across multiple samples or sequencing runs. Coverage lists enable more accurate binning downstream, especially distinguishing closely related species and resolving complex communities. 

`fairy` does not require alignments to estimate coverage, but also uses k-mer libraries for speed. It builds hash tables of subsampled k-mers from reads from each sample. It then queries these subsampled k-mers against the contig sequences to estimate how frequently contigs appear in each sample and approximates coverage.

**The code shown below is for a single timepoint (DM_C_T0) for clarity; however, I applied the same analysis across all timepoints and samples.**

First, each sample needs to be indexed: 
```
#!/bin/bash
#SBATCH --job-name=fairy_index_DM_C_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-3

# Only select lines with _R1_001.fastq.gz and pick the correct one by SLURM_ARRAY_TASK_ID
infile1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" DM_C_T0_filenames.txt)
echo "R1: $infile1"

# create variable name for the second paired end read set
infile2=$(basename "$infile1" | sed 's/_R1_host_trimmed.fastq.gz/_R2_host_trimmed.fastq.gz/')
echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_host_trimmed.fastq.gz//')
echo "outfile: $outfile"

# run fairy to index short reads
fairy sketch -1 $infile1 -2 $infile2 -d DM_C_T0_indexed_files
```

Then coverage can be estimated. I used two binning softwares (`maxbin2` and `metabat2`), which require different coverage list formats. `fairy` has functionality to make this happen:

```
#!/bin/bash
#SBATCH --job-name=fairy_coverage_DM_C_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


# run fairy to calculate coverage on indexed files
fairy coverage DM_C_T0_final.contigs.fa ./DM_C_T0_indexed_files/*.bcsp -t 8 --maxbin-format -o DM_C_T0_coverage_maxbin2.tsv
fairy coverage DM_C_T0_final.contigs.fa ./DM_C_T0_indexed_files/*.bcsp -t 8 -o DM_T0_coverage_metabat2.tsv
```
I would run one at a time, however, you may be able to run both lines together in one submission. I don't know if this will generate any errors though, so proceed with caution.


**_Maxbin2_**
After coverage lists were generated, I then started to bin the metagenomic reads into metagenome assembled genomes (MAGs) with two softwares: `maxbin2` and `metabat2`. Using multiple binning softwares is beneficial because different tools use distinct algorithms and have different assumptions, which can result in biases. In other words, one binning software may miss certain bins that another catches. Combining multiple softwares can result in:

* Increased recovery of high-quality MAGs by capturing diverse genome signatures missed by any single method.<br>
* Improve accuracy in assigning contigs into bins by leveraging different approaches.
* Enhance strain-level resolution and recovery of rare or low-abundance species by cross-validating bins across methods.

`maxbin2` uses an Expectation-Maximization (EM) algorithm that integrates multiple features to assign contigs into MAGs. It identifies conserved, single-copy marker genes within contigs that are frequently present and unique to specific prokaryotic genomes. By combining contig coverage across multiple metagenomic samples (generated with `fairy`) and nucleotide sequence composition, the EM algorithm iterativelty estimates the probability that each contig belongs to a given genomic bin and adjusts membership to maximize the overall likelihood. The output bins are **draft** MAGs, or sets of contigs clustered into bins with predicted completeness and contamination scores.

**The code shown below is for a single timepoint (DM_C_T0) for clarity; however, I applied the same analysis across all timepoints and samples.**

```
#!/bin/bash
#SBATCH --job-name=DM_C_T0_MAG_binning
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


run_MaxBin.pl -contig DM_C_T0_final.contigs.fa -abund DM_C_T0_coverage_maxbin2.tsv -min_contig_length 1500 -thread 16 -out DM_C_T0_bins
```
* `-min-contig-length`: each bin needs to have contigs greater than 1500 base pairs.


**_Metabat2_**

Similar to `maxbin2`, `metabat2` utilizes nucleotide sequence composition (tetranucleotide frequencies) and `fairy`-generated coverage lists to determine MAG binning. However, instead of EM likelihood, `metabat2` will use those factors to calculate normalized, composite scores to build a graph structure (scores are edge weights, contigs are nodes). Similar contigs will then be connected, and the software will partition similar contigs into bins. Less manual tuning is needed because `metabat2` is an adaptive algorithm.

**The code shown below is for a single timepoint (DM_C_T0) for clarity; however, I applied the same analysis across all timepoints and samples.**

```
#!/bin/bash
#SBATCH --job-name=DM_C_T0_metabat2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


metabat2 \
-i DM_C_T0_final.contigs.fa \
-o ./DM_C_T0bins \
-a DM_T0_coverage_metabat2.tsv \
-t 16
```

**_DAS Tool_**

After I obtained metagenomic bins from `maxbin2` and `metabat2`, I used `DAS Tool` (Dereplication, Aggregation, and Scoring Tool) to refine the bin results in order to generate a non-redundant, high-quality set of MAGs to use from downstream analyses. Some advantages of using bin refinement tools (also discussed above):

* Different binning tools have unique strengths and biases, producing varying bin sets with overlaps and conflicts.
* Single binning tools may miss some genomes or produce bins with contamination and completeness.
* Refinement tools can leverage complimentary bins from all tools, improving the number, quality, and accuracy of recovered MAGs.

`DAS Tool` works by first identifying redundant bins across multiple binners using single-copy marker genes. It then combines candidate bins from all input tools into a comprehensive set and scores each bin based on completeness and contamination metrics from single-copy gene presence:
* Completeness is the proportion of expected single-copy genes found in the bin.
* Contamination is indicated by multiple copies of single-copy genes, suggesting mixed genomes/poor assembly.
* The score also penalizes "megabins", which are overly large or contain many single copy genes, to avoid chimeric bins.


**The code shown below is for a single timepoint (DM_C_T0) for clarity; however, I applied the same analysis across all timepoints and samples.**

`DAS Tool` requires a specific file `contig2bin`, which is a `.tsv` file that contains the contig name is one column and its associated bin in another. `DAS Tool` has a script to do this for you, so just run:

```
Fasta_to_Contigs2Bin.sh -i /lustre/isaac24/proj/UTK0032/adekovic/BSF_metagenomics/08_maxbin2_MAG_binning -e fa > DM_C_T0_maxbin2_contig2bin.tsv
```
* `-i`: the path to the bins
* `-e`: the extension of the resulting fasta files (either fa or fasta)

Do this with bins from both programs.

```
#!/bin/bash
#SBATCH --job-name=DM_C_T0_DAS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


DAS_Tool \
-i ./DM_C_T0_maxbin2_contig2bin.tsv,./DM_C_T0_metabat2_contig2bin.tsv \
-l maxbin2,metabat \
-c ./DM_C_T0_final.contigs.fa \
-o ./DM_C_T0_refined_bins \
-t 30
```

* `-i`: provide the output files from each binning program
* `-l`: provide the name of the binning softwares, in order that matches the `contig2bin` files

**_CheckM_**

After bins are refined and selected with `DAS Tool`, I used `CheckM` to assess the quality of these MAGs for further filtering. `CheckM` identifies lineage-specific, single-copy marker genes that are expected to be present once per genome within a given phylogenetic group using a reference genome tree, [downloaded here](https://data.ace.uq.edu.au/public/CheckM_databases). It then assesses the fraction of these markers genes found in each MAG to estimate completeness (how much of the genome is recovered). Contamination is assessed based on the presence of multiple copies of single-copy marker genes, which can indicate mixed genomes or chimeric sequences. Contamination can also be reliably distinguished from closely related strains using amino acid identity of multicopy genes.

**The code shown below is for a single timepoint (DM_C_T0) for clarity; however, I applied the same analysis across all timepoints and samples.**

Before running, the database of marker genes must be downloaded. Set the appropriate path within your linux environment:

```
export CHECKM_DATA_PATH=/lustre/isaac24/proj/UTK0032/adekovic/BSF_metagenomics/11_checkM_contamination/CheckM2_database
```

```
#!/bin/bash
#SBATCH --job-name=DM_C_T0_checkM
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu



checkm2 predict --threads 30 --input ../09_metabat2_MAG_binning/DM_C_T0/bins --output-directory ./DM_C_T0_metabat2 -x fa --database_path ./CheckM2_database/uniref100.KO.1.dmnd
```

After reviewing the reports, I only kept bins with **>90% completeness and <10% contamination**.

**_GTDBTK_**

Once I finalized bins for every timepoint, I then taxonomically classified them using `gtdbtk`.

First, the database must be downloaded and the correct path needs to be set:

```
#!/bin/bash
#SBATCH --job-name=gtdbtk_database
#SBATCH --nodes=1
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu

# change into scratch directory
cd /lustre/isaac24/scratch/adekovic/gtdbtk_2.5.1_database

# Download the database in this directory
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz

# Export path to new data location
conda env config vars set GTDBTK_DATA_PATH="/lustre/isaac24/scratch/adekovic/gtdbtk_2.5.2_database/release226/"
```

Run `gtdbtk`:

```
#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long-bigmem
#SBATCH --qos=long-bigmem
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

gtdbtk classify_wf \
--genome_dir /lustre/isaac24/scratch/adekovic/gtdbtk_analyses_BSF/DM_C_T0 \
--extension fa \
--skip_ani_screen \
--out_dir /lustre/isaac24/scratch/adekovic/gtdbtk_analyses_BSF/DM_C_T0/output \
--cpus 8
```
