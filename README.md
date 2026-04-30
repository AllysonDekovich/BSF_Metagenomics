## README
This repository contains all scripts associated with the following publication (in preparation):
* Dekovich, A. & Owings, C.G. **Metagenomic evidence for microbial community structuring and pathogen reduction in livestock waste by the black soldier fly _Hermetia illucens_ (Linnaeus)(Diptera: Stratiomyidae)**

Raw metagenomic data from this study can be found [here]().

----
## Project Structure
All relevant analyses and scripts are located within the `analyses/` directory. To maintain a reproducible and organized workflow, the repository is structured as follows:
  * **Organization**: Each analysis type is contained within its own subdirectory, with labeling following the standard organizational format described in [Noble 2009](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)
  * **Script Location**: All scripts associated with a specific step (e.g., assembly, binning, trimming, etc.) are located directly within their respective analysis subdirectory.

----
## Usage & Computing Environment
The computational work in this study is divided into two primary categories:

### Command-line (SBATCH)
  Most command-line analyses were executed on ISAAC, the University of Tennessee's HPC cluster, using SLURM job scheduling/task arrays. SLURM scripts have been **anonymized**. Sensitive or system-specific directives, including partitions, email addresses, Quality of Service (QOS), and project accounts, have been removed and replaced with placeholders.
  * To replicate these analyses:
    * Input your own account information and job directives
    * Update file paths to match your local environment
    * *Note*: Original memory requirements have been retained to serve as a guide for resource allocation.
   
  ### R Scripts
  Post-metagenomic binning analyses and visualizations were performed in R. As above, file paths have been anonymized and must be replaced by the user to avoid errors.
