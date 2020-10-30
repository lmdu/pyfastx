library(ggplot2)
library(patchwork)

df1 <- read.table("matrix_fasta_build_index.tsv", header = TRUE, sep = "\t")

df2 <- df1[c("genome", "size", "count")]
df2 <- unique(df2)
#df2['label'] <- "Genome size"

p1 <-
  ggplot(df1, aes(
    x = reorder(genome, size),
    time,
    colour = tool,
    group = tool
  )) +
  geom_point(aes(size = memory)) +
  geom_line() +
  facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Elapsed time (s)") +
  guides(colour = guide_legend(title = "Tools"),
         size = guide_legend(title = "Memory (GB)"))

p2 <- ggplot(df2, aes(reorder(genome, size), size)) +
  geom_point(aes(size = count)) +
  scale_size(
    trans = "log10",
    breaks = c(10 ^ 2, 10 ^ 4, 10 ^ 6),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ylab("Genome size (GB)") +
  #facet_grid(rows = vars(label)) +
  guides(size = guide_legend(title = "Sequence count"), nrow = 3)

p3 <- ggplot(df1, aes(x = reorder(genome, size), index, group = tool)) +
  geom_bar(aes(fill = tool), stat = "identity", position="dodge") +
  #facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Index file size (MB)") +
  guides(fill = guide_legend(title = "Tools"))

df3 <- read.table("matrix_fastq_build_index.tsv", header = TRUE, sep = "\t")

df4 <- df3[c("file", "fsize", "gsize", "count")]
df4 <- unique(df4)
#df4['label'] <- "File size"

p4 <- ggplot(df3, aes(
    x = reorder(file, fsize),
    time,
    colour = tool,
    group = tool
  )) +
  geom_point(aes(size = memory)) +
  geom_line() +
  ylim(0, 3500) + 
  #facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Elapsed time (s)") +
  guides(colour = guide_legend(title = "Tools"),
         size = guide_legend(title = "Memory (GB)"))

p5 <- ggplot(df4, aes(reorder(file, fsize), fsize)) +
  geom_point(aes(size = count)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=30,hjust=1)) +
  ylab("File size (GB)") +
  #facet_grid(rows = vars(label)) +
  guides(size = guide_legend(title = "Read count"), nrow = 3)

p6 <- ggplot(df3, aes(x = reorder(file, fsize), index, group = tool)) +
  geom_bar(aes(fill = tool), stat = "identity", position="dodge") +
  #facet_grid(rows = vars(tool)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Index file size (GB)") +
  guides(fill = guide_legend(title = "Tools"))

pp <- (p1 | (p2 / p3))/(p4+p5+p6) + plot_layout(ncol = 1, heights = c(3, 1))
pp <- pp + plot_annotation(tag_levels = 'A')
show(pp)
