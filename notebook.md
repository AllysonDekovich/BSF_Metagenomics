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

**kmer estimation - alignment method was taking too long**

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

outfile=$(basename "$infile" | sed 's/.fasta//')
#echo $outfile

nonpareil -s "$infile" -T kmer -k 20 -f fastq -b "$outfile"
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
The default kmer size is 24bp, and it turns out I have some reads that were 20 bp, so the program was failing. I was able to get it working by specifying a different kmer size with `-k`.

**Creating and Interpreting the Nonpareil curves**

To plot the curves for samples on a single plot, I needed to create a text file containing the sample names, an alias for labeling, and a color hex code.

```
File	Name	Color
DM_C_T0_r1_Owings_S13_L001_R1_001.fastq.gz.npo	T0_r1	"#D27D2D"
DM_C_T0_r2_Owings_S14_L001_R1_001.fastq.gz.npo	T0_r2	"#7B3F00"
DM_C_T0_r3_Owings_S15_L001_R1_001.fastq.gz.npo	T0_r3	"#6E260E"
DM_T3_r1_Owings_S16_L001_R1_001.fastq.gz.npo	T3_r1	"#FF69B4"
DM_T3_r3_Owings_S17_L001_R1_001.fastq.gz.npo	T3_r3	"#FF00FF"
DM_T3_r4_Owings_S18_L001_R1_001.fastq.gz.npo	T3_r4	"#F8C8DC"
DM_T7_r1_Owings_S19_L001_R1_001.fastq.gz.npo	T7_r1	"#0047AB"
DM_T7_r2_Owings_S20_L001_R1_001.fastq.gz.npo	T7_r2	"#6495ED"
DM_T7_r4_Owings_S21_L001_R1_001.fastq.gz.npo	T7_r4	"#CCCCFF"
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
![Nonpareil output: coverage curves per sample](https://github.com/AllysonDekovich/BSF_Metagenomics/blob/main/DM_T0_T3_T7_coverage_curves-1.png)

* The x-axis plots the log-transformed sequencing effort, which is the number of reads or base pairs sequenced.<br/>
* The y-axis plots the estimated coverage of the community diversity or complexity on a 0 to 1 scale.

Key interpretative elements of the plots:
* Curve shape: the steeper the intial rise, the lower the community diversity.<br>
* Saturation point: where the curve plateaus indicates the maximum achievable coverage for the sequencing effort.<br>
* Diversity arrows: indicate sequence diversity (`Nd`) values for each sample.<br>
	* higher `Nd` values (skewed right) indicate more diverse samples, which require more sequencing effort than less diverse (skewed left) samples.

My results indicate moderately complex communities (steep rise around 1 Mbp and high Nd values) in each sample, and the community diversity seems to be similar in all samples (clusters of `Nd` values and curves). However, the curves did not reach a plateau (or just barely did at the end) - if I were to do MAG assembly, I would likely be able to assemble decent MAGs for abundant/common species, but I would struggle with rare species. 

I want to do both analyses (taxonomic profiling/annotation on raw reads vs. MAGs). I think the former will lead to better characterization of less abundant taxa, so I will try that first.

**_KneadData_**

Use `bowtie` to index the _Bos taurus_ reference genome prior to removing host contaminants.
```
#!/bin/bash
#SBATCH --job-name=bowtie_index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


# run bowtie-build to index the Bos genome
bowtie2-build ./bos_reference_genome/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz \
bos_taurus
```


Use `KneadData` to remove host (_Bos taurus_) contamination for more accurate microbial profiling.<br>

```
#!/bin/bash
#SBATCH --job-name=kneaddata
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
infile1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filenames.txt)
echo "R1: $infile1"

# Replace _R1_001.fastq.gz with _R2_001.fastq.gz
infile2=$(basename "$infile1" | sed 's/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz/')
echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_trimmed.fastq.gz//')
echo "outfile: $outfile"

# run kneaddata

kneaddata \
--input1 $infile1 \
--input2 $infile2 \
-db ./bos_reference_genome \
-o ${outfile}_host_trimmed
```

