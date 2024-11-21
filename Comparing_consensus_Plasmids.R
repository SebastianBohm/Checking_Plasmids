###Before starting - copy the data of your plasmid on a text document. 
#Important: Start copying at a consensus site (e.g. both times at EGFP). 
#Just add "<" and rename to fasta file. 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

library(DECIPHER)
library(ggplot2)
library(reshape2)
library(Biostrings)  

# Load sequences from FASTA files
file1 <- "C:/Users/Sebastian Böhm/Downloads/new - Copy.fa"
file2 <- "C:/Users/Sebastian Böhm/Downloads/old - Copy.fa"

seq1 <- readDNAStringSet(file1)
seq2 <- readDNAStringSet(file2)

sequences <- c(seq1, seq2)

# Perform the alignment using AlignSeqs from DECIPHER package
alignment <- DECIPHER::AlignSeqs(sequences, iterations = 100, refinements = 60)

# Convert the alignment to a character matrix
alignment_matrix <- as.matrix(alignment)

# Create a data frame for plotting
alignment_df <- as.data.frame(alignment_matrix)
alignment_df$Sequence <- rownames(alignment_df)

# Melt the data frame for ggplot2
alignment_melt <- melt(alignment_df, id.vars = "Sequence")

# Create the plot
plot <- ggplot(alignment_melt, aes(x = variable, y = Sequence, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("A" = "green", "C" = "blue", "G" = "yellow", "T" = "red", "-" = "white")) +
  theme_minimal() +
  labs(title = "Sequence Alignment", x = "Position", y = "Sequence")

# Save the plot as a JPG file
output_plot_filepath <- "C:/Users/Sebastian Böhm/Downloads/alignment_plot_NEW.jpg"
ggsave(output_plot_filepath, plot = plot, width = 10, height = 5, dpi = 300)

# Highlight mismatches in the alignment plot
alignment_mismatches <- alignment_df
alignment_mismatches$value <- ifelse(alignment_df[,1] == alignment_df[,2], "match", "mismatch")
mismatch_melt <- melt(alignment_mismatches, id.vars = "Sequence")

# Calculate pairwise distances
distance_matrix <- DECIPHER::DistanceMatrix(alignment, correction = "none")

# Extract the identity percentage
identity_percentage <- (1 - distance_matrix[1, 2]) * 100

# Print the identity percentage
cat("Identity Percentage based on pairwise distance:", identity_percentage, "%\n")

# View the alignment
DECIPHER::BrowseSeqs(alignment)
