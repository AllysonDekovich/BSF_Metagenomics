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

So far so good! I am now going to remove host contamination with `kneaddata` (see below). To double check that our sequencing effort wasn't inflated because of host contamination, I will also re-run this analysis on host trimmed reads.

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
#SBATCH --mem=70G
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

After running this analysis, host reads were successfull trimmed from all samples EXCEPT: DM_T3_r3 and DM_T7_r1. I got these error messages:
```
Error message returned from bowtie2 :
521367425 reads; of these:
  521367425 (100.00%) were unpaired; of these:
    520840270 (99.90%) aligned 0 times
    219704 (0.04%) aligned exactly 1 time
    307451 (0.06%) aligned >1 times
0.10% overall alignment rate
```
This error indicates that **none** of the reads aligned to the reference genome. There are three possible reasons that this could fail. The first being that I used an out-of-date or incorrect reference genome -- however, I do not believe this is the case, as I used the most up-to-date _Bos taurus_ reference genome and it worked well on all of the other samples. The second reason could be that the R1/R2 pair files were not matched properly. I looked back in the log file, and it looks like the proper R1/R2 files were used as input for both of these runs. The third and final reason could be biological: that these samples had minimal host DNA in general / was mostly bacterial, which would make sense as to why they did not map properly. 

I think this is the likely reason, so I will continue on with the initial files for these, but the trimmed ones for the others. I will be able to see later down the line if there are any issues. 

**`Nonpareil` coverage curves for the host trimmed samples:**
![Nonpareil output: coverage curves per sample on host trimmed samples](https://github.com/AllysonDekovich/BSF_Metagenomics/blob/main/figures/nonpareil_coverage_curve_DM_T0_T3_T7_host_trim.png)

**_MEGAHIT_**

Since the samples have more than enough coverage, I will go ahead with MAG construction. I will be using `MEGAHIT` and I will assemble a single MAG for each time point. Since we have three biological replicates for each timepoint, I can combine them all in a co-assembly. The benefits of a co-assembly are:

* **Improved Genome Recovery**: recovers a higher fraction of the genome more than a single assembly alone.<br>
* **Enhanced Continuity**: Produces less fragmented assemblies and longer contigs.<br>
* **Efficiency**: By pooling replicates into a single MAG assembly, the computational effort decreases and it is simpler to combine now than to wait until after MAG assembly.

`SLURM` scripts to run `MEGAHIT` for `T0`, `T3`, and `T7`:
```
#!/bin/bash
#SBATCH --job-name=megahit_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=260G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long-bigmem
#SBATCH --qos=long-bigmem
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


megahit \
-1 DM_C_T0_r1_Owings_S13_L001_R1_trimmed_kneaddata_paired_1.fastq,DM_C_T0_r2_Ow$
-2 DM_C_T0_r1_Owings_S13_L001_R1_trimmed_kneaddata_paired_2.fastq,DM_C_T0_r2_Ow$
-t 12 \
-m 260000000000 \
-o DM_C_T0
```

```
#!/bin/bash
#SBATCH --job-name=megahit_T3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=260G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long-bigmem
#SBATCH --qos=long-bigmem
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu



megahit \
-1 DM_T3_r1_Owings_S16_L001_R1_trimmed_kneaddata_paired_1.fastq,DM_T3_r3_Owings_S17_L001_R1_trimmed_kneaddata.trimmed.1.fastq,DM_T3_r4_Owings_S18_L001_R1_trimmed_kneaddata_paired_1.fastq \
-2 DM_T3_r1_Owings_S16_L001_R1_trimmed_kneaddata_paired_2.fastq,DM_T3_r3_Owings_S17_L001_R1_trimmed_kneaddata.trimmed.2.fastq,DM_T3_r4_Owings_S18_L001_R1_trimmed_kneaddata_paired_2.fastq \
-t 12 \
-m 260000000000 \
-o DM_T3
```

```
#!/bin/bash
#SBATCH --job-name=megahit_T7
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=260G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long-bigmem
#SBATCH --qos=long-bigmem
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


megahit \
-1 DM_T7_r1_Owings_S19_L001_R1_trimmed_kneaddata.trimmed.1.fastq,DM_T7_r2_Owings_S20_L001_R1_trimmed_kneaddata_paired_1.fastq,DM_T7_r4_Owings_S21_L001_R1_trimmed_kneaddata_paired_1.fastq \
-2 DM_T7_r1_Owings_S19_L001_R1_trimmed_kneaddata.trimmed.2.fastq,DM_T7_r2_Owings_S20_L001_R1_trimmed_kneaddata_paired_2.fastq,DM_T7_r4_Owings_S21_L001_R1_trimmed_kneaddata_paired_2.fastq \
-t 12 \
-m 260000000000 \
-o DM_T7
```

Genome fraction (1.323%) is very small, has over 3.3 million contigs, and 148 misassemblies, just to name a few. To try and improve things, I re-ran `MEGAHIT` on `DM_C_T0` with advanced parameters, such as a minimum contig length and a more diverse set of kmers. 

```
#!/bin/bash
#SBATCH --job-name=megahit_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=260G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long-bigmem
#SBATCH --qos=long-bigmem
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


megahit \
-1 DM_C_T0_r1_Owings_S13_L001_R1_trimmed_kneaddata_paired_1.fastq,DM_C_T0_r2_Owings_S14_L001_R1_trimmed_kneaddata_paired_1.fastq,DM_C_T0_r3_Owings_S15_L001_R1_trimmed_kneaddata_paired_1.fastq \
-2 DM_C_T0_r1_Owings_S13_L001_R1_trimmed_kneaddata_paired_2.fastq,DM_C_T0_r2_Owings_S14_L001_R1_trimmed_kneaddata_paired_2.fastq,DM_C_T0_r3_Owings_S15_L001_R1_trimmed_kneaddata_paired_2.fastq \
-t 12 \
-m 260000000000 \
--k-list 21,29,39,49,59,69,79,89,99,109,119,129,139,149,159,169,179,189,199 \
--min-contig-len 1000 \
-o DM_C_T0_custom
```

However, this did not improve the assembly statistics by much. I was reading online and apparently highly diverse samples can complicate these statistics. I also did a co-assembly (combined the reps within each timepoint), which may further complicate the assembler. It was recommended that I just move forward and assemble contigs into MAGs and **then** perform QUAST/BUSCO/CheckM/etc on each individual bin.

**_Maxbin2_**
I will be using `Maxbin2` to assemble contigs into separate MAGs. This program will utilize the contigs from the `MEGAHIT` assembly and the coverage information to cluster contigs into bins, each representing a putative genome. 

I need to calculate the coverage profiles of each individual sample per timepoint. I will use `bowtie2` to map back the raw reads to the new MAG assembly.

First, the 'reference' genome (DM_C_t0_final.contigs.fa) for this timepoint was indexed with `bowtie2-build`:
```
#!/bin/bash
#SBATCH --job-name=bowtie2_index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu

bowtie2-build DM_C_t0_final.contigs.fa DM_C_t0_final.contigs
```
Next, each `fastq` file was aligned to the reference genome with `bowtie2` and each bam was sorted and indexed with `samtools`:

```
#!/bin/bash
#SBATCH --job-name=bowtie2_coverage
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-3

# Only select lines with _R1_001.fastq.gz and pick the correct one by SLURM_ARRAY_TASK_ID
infile1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filenames.txt)
echo "R1: $infile1"

# create variable name for the second paired end read set
infile2=$(basename "$infile1" | sed 's/_R1_trimmed_kneaddata_paired_1.fastq/_R1_trimmed_kneaddata_paired_2.fastq/')
echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_trimmed_kneaddata_paired_1.fastq//')
echo "outfile: $outfile"

# run bowtie to calculate coverage and pipe output to bam format using samtools
bowtie2 \
-p 8 \
-x DM_C_t0_final.contigs \
-1 $infile1 \
-2 $infile2 | samtools sort -o ${outfile}.sorted.bam
```
```
#!/bin/bash
#SBATCH --job-name=index_bam_DM_C_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu

samtools index -M *.bam
```
Coverage for `DM_C_T0` timepoint was calculated with `BamM`:

```
#!/bin/bash
#SBATCH --job-name=DM_C_T0_bamm_coverage_list
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


bamm parse -c DM_C_T0_coverage.tsv -m pmean -b DM_C_T0_r1_Owings_S13_L001.sorted.bam DM_C_T0_r2_Owings_S14_L001.sorted.bam DM_C_T0_r3_Owings_S15_L001.sorted.bam
```
`bamm parse`: `BamM` command to analyze BAM files for coverage calculation
`-c`: Combines output coverage from all three replicates into a single `.tsv` file
`-m pmean`: coverage mode set to `pmean`; calculated using a particular mean method (trimmed or pooled)
`-b`: indicates the `bam` files that contain the mapped sequences to be analyzed.

Once a coverage list is generated, I will run `maxbin2`:

```
#!/bin/bash
#SBATCH --job-name=DM_C_T0_MAG_binning
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu



run_MaxBin.pl -contig DM_C_t0_final.contigs.fa \
-abund DM_C_T0_coverage_updated.tsv -min_contig_length 1500 \
-thread 8 -out DM_C_T0_bins
```

Along with providing the newly generated coverage file and the index contig assembly, I set the minimum contig length to be 1500 to ensure high quality contigs for binning later on.

**_Metabat2_**

A common recommendation is for users to run several binning programs and then apply refinement tools to retain the highest-quality bins from both for downstream analysis. So, I will also run `metabat2`:

```
#!/bin/bash
#SBATCH --job-name=DM_C_T0_metabat2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


metabat2 \
-i DM_C_t0_final.contigs.fa \
-o ./bins \
-a DM_C_T0_coverage.tsv \
-t 8
```

**_DAS Tool_**

I used `DAS` to identify the highest-quality bins from both binning softwares:

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
-i ./DM_C_T0_maxbin2_contigs2bin.tsv,./DM_C_T0_metabat2_contigs2bin.tsv \
-l maxbin2,metabat \
-c ./DM_C_t0_final.contigs.fa \
-o ./DM_C_T0_refine/DM_C_T0_DAS_refine \
-t 30
```

*`-i`: provide the output files from each binning program
*`-l`: provide the name of the binning softwares, in order that matches the `contigs2bin` files

**_CheckM_**

Once I identify the highest-quality bins for each timepoint, `CheckM` is used to identify contamination and completeness of each bin.

Before running, the `CheckM` database must be downloaded. It can be downloaded from [here](https://data.ace.uq.edu.au/public/CheckM_databases) and the path must be set within the linux environment:

```
export CHECKM_DATA_PATH=/path/to/my_checkm_data
```

Run `CheckM` on the bins from each binning software (below is code from DM_C_T0 and metabat2)
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

After reviewing the reports, I only kept bins with >90% completeness and <10% contamination.

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
export GTDBTK_DATA_PATH=/lustre/isaac24/scratch/adekovic/gtdbtk_2.5.1_database
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
