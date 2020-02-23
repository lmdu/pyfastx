library(ggplot2)
library(patchwork)

df <- read.table('reverse_complement_matrix.tsv', header=T, sep="\t")

#benchmark for index building and random access
p1 <- ggplot(df, aes(x=reorder(genome,size),time,colour=tool, group=tool)) + geom_point(aes(size=memory)) + geom_line() + facet_grid(rows=vars(tool)) + theme_bw() + theme(axis.text.x=element_text(angle=30,hjust=1), axis.title.x=element_blank()) + ylab("Elapsed time (s)") + guides(colour=guide_legend(title="Tools"), size=guide_legend(title="Memory (Gb)"))

gf <- df[c("genome","size","count")]
gf <- unique(gf)
gf['label'] <- "Genome"
p2 <- ggplot(gf, aes(reorder(genome, size), size)) + geom_point(aes(size=count)) + scale_size(trans="log10",breaks=c(10^2, 10^4, 10^6), labels=scales::trans_format("log10",scales::math_format(10^.x))) + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + ylab("Size (Gb)") + facet_grid(rows=vars(label)) + guides(size=guide_legend(title="Sequence count"), nrow=3)

p2 + p1 + plot_layout(ncol=1, heights=c(1.3,7))

#benchmark for reverse complement

sf <- df[c("genome","size","longlen")]
sf <- unique(sf)
sf['label'] <- "Longest seq"

p3 <- ggplot(sf, aes(reorder(genome,size), longlen, group=1)) + geom_point() + geom_line() + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + ylab("Length (Mb)") + facet_grid(rows=vars(label))

p2 + p3 + p1 + plot_layout(ncol=1, heights=c(1.3,1,6))
