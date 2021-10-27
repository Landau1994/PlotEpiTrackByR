###PLotEpiTrackByR
###Author:wlt
###Date:20211025
###Replicate Figure 1b. in nature 2016 The landscape of accessible chromatin in
###mammalian preimplantation embryos
###Using GEO data GSE66581 for demonstration
###Extra: using mESC from GSM1014187


###----0. Basic configuration----
library(ggplot2)
library(patchwork)
library(cowplot)
library(tidyverse)
library(GEOquery)
library(GenomicRanges)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(clusterProfiler)
library(future.apply)
library(rtracklayer)
library(GENOVA)
library(sigminer)
library(RColorBrewer)
source("myfunction_lib.R")



####define colorbar
solarExtra <- colorRampPalette(c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
hic.rdbu <- colorRampPalette(c( "#0065b2","#2399de","#67cdfe", "#f5f5f5","#fcb09d",
                                "#ed6855","#be2a21"))
hic.sadle.rdbu <- colorRampPalette(colors = c("#1f67d1","#5299f2","#69b2fd",
                                              "lightgrey",
                                              "#ffc3aa","#ff876b","red1"))
hic.contact.hot <- colorRampPalette(colors = GENOVA:::bezier_corrected_hot) 
hic.red <- colorRampPalette(c("white", "red"))
hcl.pals(type = NULL)
hcl.colors(200,palette = "YlGnBu")
blueteal <- ggthemes::tableau_color_pal(palette = "Blue-Teal",type = "ordered-sequential")
blueteal <- colorRampPalette(colors = blueteal(20))

scales::show_col(colours = hcl.colors(10,"YlGnBu"))

###---1. download data-----------
###don't run it again

####you can donwload GEO supp data by getGEOSuppFiles
tmp.id <- "GSE66581"
tmp.dir <- "data"
dir.create(path = tmp.dir)
getGEOSuppFiles(GEO = tmp.id,makeDirectory = T,baseDir = "data")
#### unzip data
untar(tarfile="data/GSE66581/GSE66581_RAW.tar",exdir="data/GSE66581/byReplicate/bedgraph")

#### or directly donwload url
tmp.input <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1014nnn/GSM1014187/suppl/GSM1014187_mm9_wgEncodeUwDnaseEscj7S129ME0SigRep2.bigWig"
tmp.output <- "data/ENCODE_mESC_DNaseI.bw" 
download.file(url=tmp.input,destfile=tmp.output)

###----2. Plot gene Track-------

###----2.1 Prepare gene track data------
mm9 <- TxDb.Mmusculus.UCSC.mm9.knownGene
mm9.tx <- myGetGRangesFromsTxDb(txdb = mm9,
                                standard.chromosomes = T,
                                verbose = T)
tmp <- bitr(mm9.tx$gene_id,
            fromType = "ENTREZID",
            toType = "SYMBOL",
            OrgDb = org.Mm.eg.db)
tmp <- tmp[!duplicated(tmp$SYMBOL),]
ttt.mm9.tx <- mm9.tx[mm9.tx$gene_id %in% tmp$ENTREZID]
ttt.mm9.tx$gene_id <- plyr::mapvalues(ttt.mm9.tx$gene_id,
                                      from = tmp$ENTREZID,
                                      to = tmp$SYMBOL)
mm9.tx <- ttt.mm9.tx

tmp.region <- "chr10:128222630-128407404"
p_gene <- myGenePlot(annotation = mm9.tx,
                     region = tmp.region,
                     arrow_sbreaks = 600,
                     font_size = 18,
                     label_size = 3)+
  ylab("Gene")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.y = element_text(angle = 0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
p_gene

###----2. plot DNase-seq track from bigwig---------

tmp.files <- "data/ENCODE_mESC_DNaseI.bw"
tmp.colors <- "#446195"
tmp.levels <- "DNase-seq"
p_DNase <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                         bigwig = tmp.files,
                         smooth = 100,
                         lognorm = F,
                         type = "coverage",
                         y_label = tmp.levels,
                         fontsize=18,
                         track.color=tmp.colors,
                         tmp.ylimits=c(0,100),
                         max.downsample = 3000,
                         downsample.rate = 0.1,
                         tmp.seed=42)
p_DNase
ggsave(filename = "fig/DNase_coverage_plot.png",
       width = 12,height = 8,dpi = 350,bg = "white")


p_DNase <- myBigwigTrack(region = as(tmp.region,"GRanges"),
              bigwig = tmp.files,
              smooth = 100,
              lognorm = F,
              type = "bar",
              y_label = tmp.levels,
              fontsize=18,
              track.color=tmp.colors,
              tmp.ylimits=c(0,100),
              max.downsample = 3000,
              downsample.rate = 0.1,
              tmp.seed=42)
p_DNase
ggsave(filename = "fig/DNase_bar_plot.png",
       width = 12,height = 8,dpi = 350,bg = "white")


###----3. Plot ATAC-seq track------------

###----3.1. Convert bedgraph file to bigwig--------

###unzip file
tmp.dir <- "data/GSE66581/byReplicate/bedgraph/"
tmp.files <- list.files(path = tmp.dir,full.names = T)
plan(multisession, workers = 8)
future_lapply(seq_along(tmp.files),FUN = function(x){
  gunzip(filename = tmp.files[x])
})

###rename file, remove GSMXXXXX
tmp.dir <- "data/GSE66581/byReplicate/bedgraph/"
tmp.files <- list.files(path = tmp.dir,full.names = T)
future_lapply(seq_along(tmp.files),FUN = function(x){
  file.rename(from = tmp.files[x],
              to = gsub(pattern = "GSM[0-9]{7}_",
                        replacement = "",x = tmp.files[x]))
})


###create out dir
dir.create("data/GSE66581/byReplicate/bigwig")
dir.create("data/GSE66581/bigwig")

###bedgraph to bigwig
###write bw
tmp.dir <- "data/GSE66581/byReplicate/bedgraph/"
tmp.files <- list.files(path = tmp.dir,full.names = T)
future_lapply(seq_along(tmp.files),FUN = function(x){
  tmp.sample <- gsub(pattern = "data/GSE66581/byReplicate/bedgraph/",
                     replacement = "",
                     x = tmp.files[x])
  tmp.sample <- gsub(pattern = ".bedGraph",
                     replacement = ".bw",
                     x = tmp.sample)
  tmp.output <- paste0("data/GSE66581/byReplicate/bigwig/",tmp.sample)
  myBegraphToBigwig(tmp.input = tmp.files[x],
                    tmp.output = tmp.output,
                    region = NULL,
                    tmp.seqlength = seqlengths(mm9.tx),
                    tmp.levels = seqlevels(mm9.tx))
})

###non replicate
tmp.dir <- "data/GSE66581/"
tmp.files <- list.files(path = tmp.dir,
                        pattern = ".bedGraph.gz",
                        full.names = T)
future_lapply(seq_along(tmp.files),FUN = function(x){
  gunzip(filename = tmp.files[x])
})

tmp.dir <- "data/GSE66581/"
tmp.files <- list.files(path = tmp.dir,
                        pattern = ".bedGraph",
                        full.names = T)
future_lapply(seq_along(tmp.files),FUN = function(x){
  file.rename(from = tmp.files[x],
              to = gsub(pattern = "GSE66581_atac_",
                        replacement = "",
                        x = tmp.files[x]))
})


tmp.dir <- "data/GSE66581"
tmp.files <- list.files(path = tmp.dir,
                        pattern = ".bedGraph",
                        full.names = T)
future_lapply(seq_along(tmp.files),FUN = function(x){
  tmp.sample <- gsub(pattern = "data/GSE66581/",
                     replacement = "",
                     x = tmp.files[x])
  tmp.sample <- gsub(pattern = ".bedGraph",
                     replacement = ".bw",
                     x = tmp.sample)
  tmp.output <- paste0("data/GSE66581/bigwig/",
                       tmp.sample)
  myBegraphToBigwig(tmp.input = tmp.files[x],
                    tmp.output = tmp.output,
                    region = NULL,
                    tmp.seqlength = seqlengths(mm9.tx),
                    tmp.levels = seqlevels(mm9.tx))
})


plan(strategy = sequential)

####---------3.2 plot ATAC-seq track--------------

tmp.levels <- c("2cell","4cell","8cell","icm")
tmp.colors <- c("#3d7b42","#6a4181","#cf8e67","#a9477e")
scales::show_col(tmp.colors)

tmp.list.1 <- lapply(seq_along(tmp.levels),function(ii){
  cat(ii,sep = "\n")
  p1 <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                   bigwig = paste0("data/GSE66581/byReplicate/bigwig/",tmp.levels[ii],"_re1",".bw"),
                   smooth = 100,
                   lognorm = F,
                   type = "bar",
                   y_label = paste0(tmp.levels[ii]," rep1"),
                   fontsize=18,
                   track.color=tmp.colors[ii],
                   tmp.ylimits=c(0,10),
                   max.downsample = 3000,
                   downsample.rate = 0.1,
                   tmp.seed=42)
  p2 <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                    bigwig = paste0("data/GSE66581/byReplicate/bigwig/",tmp.levels[ii],"_re2",".bw"),
                    smooth = 100,
                    lognorm = F,
                    type = "bar",
                    y_label =paste0(tmp.levels[ii]," rep2"),
                    fontsize=18,
                    track.color=tmp.colors[ii],
                    tmp.ylimits=c(0,10),
                    max.downsample = 3000,
                    downsample.rate = 0.1,
                    tmp.seed=42)
  p <- list(p1,p2)
  return(p)
})
tmp.list.1 <- Reduce(c,tmp.list.1)

####tmp.list.2
tmp.levels <- c("200_mESC","1k_mESC","50k_mESC")
tmp.colors <- c("#76492a","#76492a","#76492a")
scales::show_col(tmp.colors)

tmp.list.2 <- lapply(seq_along(tmp.levels),function(ii){
  cat(ii,sep = "\n")
  p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                      bigwig = paste0("data/GSE66581/bigwig/",tmp.levels[ii],".bw"),
                      smooth = 100,
                      lognorm = F,
                      type = "bar",
                      y_label = tmp.levels[ii],
                      fontsize=18,
                      track.color=tmp.colors[ii],
                      tmp.ylimits=c(0,10),
                      max.downsample = 3000,
                      downsample.rate = 0.1,
                      tmp.seed=42)
  return(p)
})

p_ATAC_track_list <- c(tmp.list.1,tmp.list.2)
plot.list <- c(p_ATAC_track_list,list(p_DNase))


(wrap_plots(plot.list,ncol = 1)  + 
    p_gene  + plot_layout(ncol = 1,heights = c(rep(1,12),3))) + 
  plot_annotation(title = "chr10:128,222,630–128,407,404") & 
  theme(plot.title = element_text(hjust = 1),
        axis.line.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank())

###create output directory
dir.create(path = "fig")
ggsave(filename = "fig/ATAC-seq_track_by_R.png",
       width = 12,height = 8,
       dpi = 350)


####-------4.add pyramid track-------

####-------4.1 load hic data -------
####The ICN hic data downloaded and converted from GSE82185
####load centeromere file
data("centromeres.mm9")
hic_ICM_10k <- load_contacts(signal_path = "data/ICM.hic", 
                               resolution = 10e3, 
                               balancing = T,
                               z_norm = F, 
                               sample_name = "ICM",
                               centromeres = centromeres.mm9,
                               colour = "black")

###-------4.2 plot pyramid ------------------
tmp.crop <- 1e6
p_hic <- pyramid(exp = hic_ICM_10k,
        chrom = "chr10",
        start = 94.08e6,
        end = 97.08e6,
        crop_y=c(0,tmp.crop),
        display_yaxis=T,
        colour = scale_fill_gradientn(name ="Contacts",
                                      colours = rev(hcl.colors(100,palette = "YlGnBu")),
                                      limits=c(0,80)))+
  theme(plot.title = element_text(face="bold",),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="bold"),
        axis.text = element_blank(),
        axis.title.y = element_text(angle = 0,size=18),
        axis.title.x = element_blank(),
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))
p_hic$facet$y_scales$.pyramid$name <- "ICM"
p_hic
ggsave(filename = "fig/signle_pyramid.png",width = 8,height = 8,dpi = 350,bg = "white")





###---- 4.3 plot mutiple track-------

tmp.region <- "chr10:94080000-97080000"
p_gene <- myGenePlot(annotation = mm9.tx,
                     region = tmp.region,
                     arrow_sbreaks = 80000,
                     font_size = 18,showlabel = F,
                     label_size = 3)+
  ylab("Gene")+
  scale_x_continuous(breaks = c(9.5e7,9.6e7,9.7e7),labels=c("95","96","97"))+
  xlab("chr10 position (Mb)")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.y = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
p_gene


tmp.levels <- c("2cell","4cell","8cell","ICM")
tmp.colors <- c("#3d7b42","#6a4181","#cf8e67","#a9477e")
scales::show_col(tmp.colors)


tmp.list <- lapply(seq_along(tmp.levels),function(ii){
  cat(ii,sep = "\n")
  p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                     bigwig = paste0("data/GSE66581/bigwig/",tmp.levels[ii],".bw"),
                     smooth = 500,
                     lognorm = T,
                     type = "coverage",
                     y_label = tmp.levels[ii],
                     fontsize=18,
                     track.color=tmp.colors[ii],
                     tmp.ylimits=c(0,1),
                     max.downsample = 3000,
                     downsample.rate = 0.1,
                     tmp.seed=42)+
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank())
  return(p)
})



p_hic + wrap_plots(tmp.list,ncol = 1) + 
  p_gene + 
  plot_layout(ncol = 1) + 
  plot_annotation(title = tmp.region) &
  theme(plot.title = element_text(face = "bold",
                                  hjust = 0.5))

ggsave(filename = "fig/pyramid_hic.png",width = 12,height = 6,dpi = 350)
