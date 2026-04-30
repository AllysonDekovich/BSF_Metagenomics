library(Nonpareil)

# read in the sample file, which gives each sample a color for easy identification.
samples <- read.table('npo_samples.txt', sep='\t', header=TRUE, as.is=TRUE)

# Create a single plot that visualizes the metagenomic coverage for all samples at once.
np_curves <- Nonpareil.set(as.vector(samples$File), col=samples$Color, 
labels=samples$Name, plot.opts=list(plot.observed=FALSE));
