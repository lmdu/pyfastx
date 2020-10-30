library(ggplot2)
library(patchwork)

df <- read.table("matrix_fasta_random_access.tsv", header = TRUE, sep = "\t")
#df1 <- df[which(df$tool != 'pyfasta' & df$tool != 'pyfastx_gzip' & df$tool != 'biopython'),]
#df2 <- df[which(df$tool == 'pyfasta' | df$tool == 'pyfastx_gzip' | df$tool == 'biopython'),]

p1 <-
  ggplot(df, aes(
    x = reorder(genome, size),
    time,
    colour = tool,
    group = tool
  )) +
  geom_point(aes(size = memory)) +
  geom_line() +
  ylim(0, 1000) +
  facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Elapsed time (s)") +
  guides(colour = guide_legend(title = "Tools"),
         size = guide_legend(title = "Memory (GB)"))

#p2 <-
#  ggplot(df2, aes(
#    x = reorder(genome, size),
#    time,
#    colour = tool,
#    group = tool
#  )) +
#  geom_point(aes(size = memory)) +
#  geom_line() +
  #facet_grid(rows = vars(tool), scales="free_y") +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 30, hjust = 1),
#        axis.title.x = element_blank()) +
#  ylab("Elapsed time (s)") +
#  guides(colour = guide_legend(title = "Tools"),
#         size = guide_legend(title = "Memory (GB)"))

df2 <- read.table("matrix_fasta_extract_subsequences.tsv", header=TRUE, sep="\t")
p2 <-
  ggplot(df2, aes(
    x = reorder(genome, size),
    time,
    colour = tool,
    group = tool
  )) +
  geom_point(aes(size = memory)) +
  geom_line() +
  ylim(0, 50) + 
  facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Elapsed time (s)") +
  guides(colour = guide_legend(title = "Tools"),
         size = guide_legend(title = "Memory (GB)"))

df3 <- read.table("matrix_fastq_random_access.tsv", header = TRUE, sep = "\t")

p3 <-
  ggplot(df3, aes(
    x = reorder(file, fsize),
    time,
    colour = tool,
    group = tool
  )) +
  geom_point(aes(size = memory)) +
  geom_line() +
  ylim(0, 500) + 
  facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Elapsed time (s)") +
  guides(colour = guide_legend(title = "Tools"),
         size = guide_legend(title = "Memory (GB)"))

pp <- p1 | p2 | p3
pp <- pp + plot_annotation(tag_levels = 'A')
show(pp)
