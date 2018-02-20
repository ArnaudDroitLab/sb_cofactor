#setwd("C:/Dev/Projects/sb_cofactor/A549")

# Load libraries
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
library(ChIAnalysis)
library(ENCODExplorer)
library(R.utils)

output.dir="output/HiC"
dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output.dir, "dex"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output.dir, "ct"), recursive=TRUE, showWarnings=FALSE)


# Load DB regions
csaw_regions = GRanges(read.table("output/analyses/ENCODE-csaw-results.txt", sep="\t", header=TRUE))

# Load and annotate HiC contacts
hic_contacts_dex = read.table("input/ENCFF385DHX.tsv", sep="\t", header=TRUE)[,c(1:6, 8)]
hic_contacts_ctl = read.table("input/ENCFF803ZOW.tsv", sep="\t", header=TRUE)[,c(1:6, 8)]
 
dex_chia = load_chia(input.df=hic_contacts_dex, excluded.chr="chrM")
ctl_chia = load_chia(input.df=hic_contacts_ctl, excluded.chr="chrM")

de_results = read.table("results/a549_dex_time_points/6h", header=TRUE, sep=",")
de_results$padj[is.na(de_results$padj)] = 1
de_indices = de_results$padj <= 0.05 & abs(de_results$log2FoldChange) >= log2(1.5)
de_6h = de_results$log2FoldChange[de_indices]
names(de_6h) = de_results$symbol[de_indices]

chia.params = build_chia_params(biosample = "A549",
                                genome.build = "GRCh38",
                                centrality.measures=c("Degree", "Closeness"),
                                weight.attr="Reads",
                                tssRegion = c(-1000, 1000),
                                gene.annotations = list(DE_6h=de_6h))
#cache(chia.params <- add_encode_data(chia.params), dir=output.dir)
cache(chia.params <- add_encode_data(chia.params, treatment=NA), dir=output.dir)

# Read GR binding sites.
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
gr_dex = import("output/ENCODE-chip/peak_call/GR_1hr/GR_1hr_peaks.narrowPeak", format = "BED", extraCols = extraCols_narrowPeak)
gr_ctl = import("output/ENCODE-chip/peak_call/GR_0hr/GR_0hr_peaks.narrowPeak", format = "BED", extraCols = extraCols_narrowPeak)
mcols(gr_dex) = data.frame(signalValue=mcols(gr_dex)$signalValue)
mcols(gr_ctl) = data.frame(signalValue=mcols(gr_ctl)$signalValue)
chia.params$tf.regions = c(chia.params$tf.regions, GRangesList(GR_1hr=gr_dex, GR_0hr=gr_ctl))

# Import POLR2A binding under dexamethasone condition.
dex_polr2a_bed = downloadEncode(file_acc="ENCFF825AZS", dir="input/")
dex_polr2a_bed = gunzip(dex_polr2a_bed)
polr2a_dex = import(dex_polr2a_bed, format = "BED", extraCols = extraCols_narrowPeak)
mcols(polr2a_dex) = data.frame(signalValue=mcols(polr2a_dex)$signalValue)
chia.params$pol.regions = c(chia.params$pol.regions, GRangesList(POLR2A_Dex=polr2a_dex))

# Import loftOver'ed chromatin states
chia.params$input.chrom.state = import("input/E114_18_core_K27ac_mnemonics_hg38_liftOver.bed", format="BED")
chia.params$input.chrom.state$name = gsub("^(\\d)_", "0\\1_", as.character(chia.params$input.chrom.state$name))

dir.create(file.path(output.dir, "dex"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output.dir, "ctl"), recursive=TRUE, showWarnings=FALSE)
cache(dex_chia <- annotate_chia(dex_chia, chia.params, output.dir=file.path(output.dir, "dex")), dir=file.path(output.dir, "dex"))
cache(ctl_chia <- annotate_chia(ctl_chia, chia.params, output.dir=file.path(output.dir, "ctl")), dir=file.path(output.dir, "ctl"))

#normed.states = gsub("^(\\d)_", "0\\1_", as.character(dex_chia$Regions$Chrom.State))
dex_chia$Regions$Chrom.State = normed.states

# Get degree by gene.
dex_genes = select_gene_representative(dex_chia)
dex_degree = data.frame(SYMBOL=dex_genes$Regions$SYMBOL, DEX=dex_genes$Regions$Degree)

ctl_genes = select_gene_representative(ctl_chia)
ctl_degree = data.frame(SYMBOL=ctl_genes$Regions$SYMBOL, CTL=ctl_genes$Regions$Degree)

merged_degree = merge(dex_degree, ctl_degree, by="SYMBOL", all=TRUE)
merged_degree$DEX[is.na(merged_degree$DEX)] = 0
 merged_degree$CTL[is.na(merged_degree$CTL)] = 0

ggplot(merged_degree, aes(x=CTL, y=DEX)) + 
  geom_point() +
  stat_density2d(aes(fill=..density..), geom = "tile", contour = FALSE, n=15) +
  scale_fill_gradient2(low = "white", high = "red")

merged_degree$Diff = merged_degree$DEX - merged_degree$CTL
merged_degree$DE = merged_degree$SYMBOL %in% names(de_6h)
merged_degree$FC = de_6h[merged_degree$SYMBOL]

merged_degree$DE_Category = ifelse(!(merged_degree$SYMBOL %in% names(de_6h)), "STABLE", 
                                ifelse(merged_degree$SYMBOL %in% names(de_6h) & de_6h[as.character(merged_degree$SYMBOL)] > 0, "UPREGULATED", "DOWNREGULATED"))
  
ggplot(merged_degree, aes(x=DE_Category, y=Diff)) + 
  geom_violin()
  
t.test(merged_degree$Diff[merged_degree$DE_Category == "STABLE"], merged_degree$Diff[merged_degree$DE_Category != "STABLE"])
  
bah = rbind(hic_contacts_dex, hic_contacts_ctl)[,c(1:6, 8)]
test.chia = load_chia(bah, "chrM")

hic_left_dex = GRanges(data.frame(seqnames=hic_contacts_dex$chr1, start=hic_contacts_dex$x1, end=hic_contacts_dex$x2))
hic_right_dex = GRanges(data.frame(seqnames=hic_contacts_dex$chr2, start=hic_contacts_dex$y1, end=hic_contacts_dex$y2)) 

hic_left_ctl = GRanges(data.frame(seqnames=hic_contacts_ctl$chr1, start=hic_contacts_ctl$x1, end=hic_contacts_ctl$x2))
hic_right_ctl = GRanges(data.frame(seqnames=hic_contacts_ctl$chr2, start=hic_contacts_ctl$y1, end=hic_contacts_ctl$y2)) 

hic_left_dex_olap = findOverlaps(hic_left_dex, GRanges(test.chia$Regions))
hic_right_dex_olap = findOverlaps(hic_right_dex, GRanges(test.chia$Regions))
hic_left_ctl_olap = findOverlaps(hic_left_ctl, GRanges(test.chia$Regions))
hic_right_ctl_olap = findOverlaps(hic_right_ctl, GRanges(test.chia$Regions))


dex_contacts = data.frame(Left=subjectHits(hic_left_dex_olap), Right=subjectHits(hic_right_dex_olap))
ctl_contacts = data.frame(Left=subjectHits(hic_left_ctl_olap), Right=subjectHits(hic_right_ctl_olap))

dex_contacts_str = paste(c(dex_contacts$Left, dex_contacts$Right), c(dex_contacts$Right, dex_contacts$Left), sep="-")
ctl_contacts_str = paste(ctl_contacts$Left, ctl_contacts$Right, sep="-")


# Naive comparison.
load_chia <- function(chia.raw, excluded.chr=c()) {
    chia.raw = chia.raw[,1:7]
    colnames(chia.raw) = c("L.chr", "L.start", "L.end", "R.chr", "R.start", "R.end", "Reads")
    
    chia.raw = chia.raw[!(chia.raw$L.chr %in% excluded.chr) & !(chia.raw$R.chr %in% excluded.chr),]
    
    # Separate and extend on both sides
    split.raw.chia <- function(chia.raw, columns, flank.size=0) {
       result = chia.raw[,columns]
       colnames(result) = c("chr", "start", "end")
       result$start = pmax(0, result$start - flank.size)
       result$end = result$end + flank.size

       return(GRanges(result))
    }

    chia_left.ranges = split.raw.chia(chia.raw, 1:3)
    chia_right.ranges = split.raw.chia(chia.raw, 4:6)

    # Reduce to a single set of coordinates
    single.set = reduce(unlist(GRangesList(chia_left.ranges, chia_right.ranges)))

    # Build a graph.
    # Map back to original contact points
    chia_left.indices = findOverlaps(chia_left.ranges, single.set, select="first")
    chia_right.indices = findOverlaps(chia_right.ranges, single.set, select="first")

    # Find and remove self loops.
    mapped.df = cbind(chia.raw, Left=chia_left.indices, Right=chia_right.indices)
    mapped.df = mapped.df[mapped.df$Left != mapped.df$Right,]

    # Summarize multiple edges.
    max.df = ddply(mapped.df, ~Left*Right, summarize, L.chr=head(L.chr, n=1), L.start=min(L.start), L.end=max(L.end),
                                                      R.chr=head(R.chr, n=1), R.start=min(R.start), R.end=max(R.end), Reads=sum(Reads))

    # Remap IDs: igraph will create as many nodes as max(ids), which will cause problems for us.
    only.unique = sort(unique(c(max.df$Left, max.df$Right)))
    remapped.ids = data.frame(Original=only.unique, Remapped=1:length(only.unique))
    new.left = remapped.ids$Remapped[match(max.df$Left, remapped.ids$Original)]
    new.right = remapped.ids$Remapped[match(max.df$Right, remapped.ids$Original)]
                                                      
    # Create iGraph object and set the original coordinates and the number of supporting reads as edge attributes.
    chia.graph = make_graph(c(rbind(new.left, new.right)), directed=FALSE)
    edge_attr(chia.graph) <- max.df

    # Regions which were only part of a self-loop will have been filtered above
    # and will not be part of the graph object. This will cause a difference between
    # the lengths of single.set and chia.graph.
    # To fix this, we remove all regions which are not in max.df.
    single.set = single.set[only.unique]
    
    chia.obj = list(Regions=as.data.frame(single.set), Graph=chia.graph)
    class(chia.obj) = "ChIA"
    
    return(chia.obj)
}



library(Homo.sapiens)
library(ggbio)
TxDb(Homo.sapiens) = TxDb.Hsapiens.UCSC.hg38.knownGene

plot_range <- function(chia.obj, plot_range, chia.params, tf_vector, pol_vector, file.out) {
    # Plot a track for the genes.
    #vegfa.genes = ggbio::autoplot(TxDb.Hsapiens.UCSC.hg19.knownGene, wh=vegfa.range)
    component.genes = ggbio::autoplot(Homo.sapiens, wh=plot_range)
    
    # Plot interactions.
    left.middle = start(GRanges(chia_left(chia.obj))) + width(GRanges(chia_left(chia.obj))) / 2
    right.middle = start(GRanges(chia_right(chia.obj))) + width(GRanges(chia_right(chia.obj))) / 2
    seq.middle = as.character(seqnames(GRanges(chia_left(chia.obj))))
    component.arch.ranges = GRanges(data.frame(seqnames=seq.middle, start=left.middle, end=right.middle))
    
    component.arches = ggplot() + geom_arch(data=component.arch.ranges)
    
    # Plot flat regions
    component.chia.regions = ggbio::autoplot(get_granges(chia.obj), wh=plot_range)
    
    # Plot transcription factors
    tf.plots = list()
    for(tf in tf_vector) {
        tf.intersect = intersect(chia.params$tf.regions[[tf]], plot_range)
        if(length(tf.intersect) > 0) {
            tf.plots[[tf]] = ggbio::autoplot(intersect(chia.params$tf.regions[[tf]], plot_range), wh=plot_range)
        }
    }
 
    # Plot polymerases
    pol.plots = list()
    for(pol in pol_vector) {
        pol.intersect = intersect(chia.params$pol.regions[[pol]], plot_range)
        if(length(pol.intersect) > 0) {
            pol.plots[[pol]] = ggbio::autoplot(intersect(chia.params$pol.regions[[pol]], plot_range), wh=plot_range)
        }
    } 
    cs.plot = list()
    if(!is.null(chia.params$input.chrom.state)) {
        cs_subset = chia.params$input.chrom.state[subjectHits(findOverlaps(plot_range, chia.params$input.chrom.state))]
        cs.plot[[1]] = autoplot(cs_subset, aes(color=name, fill=name))
    }
 
    pdf(file.out, width=14, height=7)
    track.list = c(list(Interactions=component.arches, 
                        Regions=component.chia.regions, 
                        Genes=component.genes),
                   tf.plots,
                   pol.plots,
                   cs.plot)
    track.heights = c(1, 1, 3, rep(1, length(tf.plots)), rep(1, length(pol.plots)), 1)
    
    print(do.call(tracks, c(track.list, list(heights=track.heights))) + theme_bw() + theme(panel.border = element_blank()))
    dev.off()
}

plot_component <- function(component.id, chia.obj, chia.params, tf_vector, pol_vector, file.out="Components.pdf") {
    # Select component.
    chia.obj = select_by_components(chia.obj, component.id)
    
    # Determine the range.
    component.range = range(get_granges(chia.obj), ignore.strand=TRUE)
    plot_range(chia.obj, component.range, chia.params, tf_vector, pol_vector, file.out)
}

plot_gene_component <- function(gene_name, chia.obj, file.out) {
    component.id = unique(chia.obj$Regions$Component.Id[chia.obj$Regions$SYMBOL == gene_name])
    plot_component(component.id, chia.obj, chia.params, "NR3C1", "POLR2A", file.out)
}

chia_subset_range <- function(chia.obj, gr) {
    return(chia_vertex_subset(chia.obj, subjectHits(findOverlaps(gr, get_granges(chia.obj)))))
}


dex_FKBP5 = unique(dex_chia$Regions$Component.Id[grepl("FKBP5", dex_chia$Regions$SYMBOL)])
plot_component(dex_FKBP5, dex_chia, chia.params, "GR_1hr", "POLR2A_Dex", file.out="A549_DEX_FKBP5.pdf")
#ctl_FKBP5 = unique(ctl_chia$Regions$Component.Id[grepl("FKBP5", ctl_chia$Regions$SYMBOL)])
#plot_component(ctl_FKBP5, ctl_chia, chia.params, "GR_0hr", "POLR2A", file.out="A549_CTL_FKBP5.pdf")

dex_FKBP5_range = range(get_granges(select_by_components(dex_chia, dex_FKBP5)), ignore.strand=TRUE)
ctl_subset = chia_subset_range(ctl_chia, dex_FKBP5_range)
plot_range(ctl_subset, dex_FKBP5_range, chia.params, "GR_0hr", "POLR2A", file.out="A549_CTL_FKBP5.pdf")

dex_TGFB2 = unique(dex_chia$Regions$Component.Id[grepl("TGFB2", dex_chia$Regions$SYMBOL)])
plot_component(dex_TGFB2, dex_chia, chia.params, "GR_1hr", "POLR2A_Dex", file.out="A549_DEX_TGFB2.pdf")
#ctl_TGFB2 = unique(ctl_chia$Regions$Component.Id[grepl("TGFB2", ctl_chia$Regions$SYMBOL)])
#plot_component(ctl_TGFB2, ctl_chia, chia.params, "GR_0hr", "POLR2A", file.out="A549_CTL_TGFB2.pdf")

dex_TGFB2_range = range(get_granges(select_by_components(dex_chia, dex_TGFB2)), ignore.strand=TRUE)
ctl_subset = chia_subset_range(ctl_chia, dex_TGFB2_range)
plot_range(ctl_subset, dex_TGFB2_range, chia.params, "GR_0hr", "POLR2A", file.out="A549_CTL_TGFB2.pdf")

autoplot(chia.params$input.chrom.state[subjectHits(findOverlaps(dex_FKBP5_range, chia.params$input.chrom.state))])


de_results$GR_Bound_Dex = dex_genes$Regions$TF.overlap.NR3C1[match(de_results$symbol, dex_genes$Regions$SYMBOL)] > 0
de_results$GR_Bound_Ctl = ctl_genes$Regions$TF.overlap.NR3C1[match(de_results$symbol, ctl_genes$Regions$SYMBOL)] > 0