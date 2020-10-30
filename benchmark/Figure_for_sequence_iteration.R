library(ggplot2)
library(patchwork)

df1 <- read.table("matrix_fasta_sequence_iteration.tsv", header = TRUE, sep = "\t")

#df2 <- df1[which(df1$tool != 'pyfasta'),]
#df3 <- df1[which(df1$tool == 'pyfasta'),]

p1 <-
  ggplot(df1, aes(
    x = reorder(genome, size),
    time,
    colour = tool,
    group = tool
  )) +
  geom_point(aes(size = memory)) +
  geom_line() +
  ylim(0,500) +
  facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Elapsed time (s)") +
  guides(colour = guide_legend(title = "Tools"),
         size = guide_legend(title = "Memory (GB)"))

df2 <- read.table("matrix_fastq_sequence_iteration.tsv", header = TRUE, sep = "\t")
p2 <-
  ggplot(df2, aes(
    x = reorder(file, fsize),
    time,
    colour = tool,
    group = tool
  )) +
  geom_point(aes(size = memory)) +
  geom_line() +
  ylim(0,2500) +
  facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Elapsed time (s)") +
  guides(colour = guide_legend(title = "Tools"),
         size = guide_legend(title = "Memory (GB)"))

pp <- p1 | p2
pp <- pp + plot_annotation(tag_levels = 'A')
show(pp)