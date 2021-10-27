###-----Data process------
myTrunclog <- function(x){
  res <- ifelse(x>1,log10(x),0)
}

###-----IO function-------

myBegraphToBigwig <- function(tmp.input,tmp.output,region=NULL,tmp.seqlength=NULL,tmp.levels=NULL){
  cat("load bedgraph",sep = "\n")
  if(!is.null(region)){
    region_data <- rtracklayer::import.bedGraph(
      con = tmp.input,
      which = as(region,"GRanges")
    )
  }else{
    region_data <- rtracklayer::import.bedGraph(
      con = tmp.input,
    )
  }
  #rtracklayer::close(con = tmp.input)
  seqlevels(region_data) <- tmp.levels
  seqlengths(region_data) <- tmp.seqlength
  cat("export bedgraph",sep = "\n")
  rtracklayer::export.bw(object = region_data,con = tmp.output)
  #rtracklayer::close(con = tmp.output)
}

###-----Visualization function---------------
myGetGRangesFromsTxDb <- function(
	txdb = mm9KG_txdb,
	standard.chromosomes = TRUE,
	verbose = TRUE
) {
	if (!requireNamespace("biovizBase", quietly = TRUE)) {
		stop("Please install biovizBase\n",
			 "https://www.bioconductor.org/packages/biovizBase/")
	}
	# convert seqinfo to granges
	whole.genome <-  as(object = seqinfo(x = txdb), Class = "GRanges")
	if (standard.chromosomes) {
		whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
	}
	# extract genes from each chromosome
	if (verbose) {
		####crunch featching GRanges from various data source
		tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
			cat(seqlevels(whole.genome)[x],sep = "\n")
			biovizBase::crunch(
				obj = txdb,
				which = whole.genome[x],
				columns = c("tx_id", "tx_name", "gene_id"))
		})
	} else {
		tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
			suppressMessages(expr = biovizBase::crunch(
				obj = txdb,
				which = whole.genome[x],
				columns = c("tx_id", "tx_name","gene_id")))
		})
	}
	# combine
	tx <- do.call(what = c, args = tx)
	return(tx)
}

myBigwigTrack <- function(
	region = tmp.region,
	bigwig = K4me3,
	smooth = 200,
	lognorm = F,
	type = "coverage",
	y_label = "Score",
	fontsize=18,
	track.color="blue",
	tmp.ylimits=c(0,16),
	max.downsample = 3000,
	downsample.rate = 0.1,
	tmp.seed=42
) {
	possible_types <- c("line", "heatmap", "coverage","bar")
	if (!(type %in% possible_types)) {
		stop(
			"Invalid type requested. Choose ",
			paste(possible_types, collapse = ", ")
		)
	}
	if (!requireNamespace("rtracklayer", quietly = TRUE)) {
		message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
		return(NULL)
	}
	if (!inherits(x = region, what = "GRanges")) {
		region <- as(region,"GRanges")
	}
	### to smooth track data
	region_data <- rtracklayer::import(
		con = bigwig,
		which = region,
		as = "NumericList"
	)[[1]]
	if (!is.null(x = smooth)) {
		region_data <- RcppRoll::roll_mean(x = region_data, n = smooth, fill = 0L)
	}
	region_data <- data.frame(
		position = start(x = region):end(x = region),
		score = region_data,
		stringsAsFactors = FALSE
	)
	window.size = width(x = region)
	sampling <- min(max.downsample, window.size * downsample.rate)
	set.seed(tmp.seed)
	coverages <- slice_sample(.data = region_data, n = sampling)
	if(lognorm){
		coverages$score <- myTrunclog(coverages$score)
	}
	p <- ggplot(
		data = coverages,
		mapping = aes_string(x = "position", y = "score")
	)
	if (type == "line") {
		p <- p + geom_line()
	} else if (type == "heatmap") {
		# different downsampling needed for heatmap
		# cut into n bins and average within each bin
		region_data$bin <- floor(x = region_data$position / smooth)
		region_data <- group_by(region_data, bin)
		region_data <- mutate(region_data, score = mean(x = score))
		region_data <- ungroup(region_data)
		region_data <- unique(x = region_data[, c("bin", "score")])
		p <- ggplot(
			data = region_data,
			mapping = aes_string(x = "bin", y = 1, fill = "score")
		) + geom_tile() + scale_fill_viridis_c()
	} else if (type == "coverage") {
		p <- p + 
			geom_area(fill=track.color,color=track.color)
	} else if (type == "bar") {
		p <- p + 
			geom_bar(stat="identity",fill=track.color,color=track.color)
	} 
	chromosome <- as.character(x = seqnames(x = region))
	p <- p  +
		xlab(label = paste0(chromosome, " position (bp)")) +
		ylab(label = y_label) +
		scale_x_continuous(expand = c(0,0),limits = c(start(region),end(region)))+
		scale_y_continuous(expand = c(0,0),limits = tmp.ylimits,
						   breaks = c(tmp.ylimits[2]),
						   labels = paste0("[",0,"-",tmp.ylimits[2],"]"))+
		theme_cowplot(font_size = fontsize)+
		theme(axis.line = element_line(size = 0.5),
			  axis.ticks = element_line(size = 0.5),
			  axis.title.y = element_text(angle = 0),
			  axis.line.y = element_blank(),
			  axis.ticks.y = element_blank())
	return(p)
}



split_body <- function(df, width = 1000) {
	wd <- df$end - df$start
	nbreak <- wd / width
	if (nbreak > 1) {
		steps <- 0:(nbreak)
		starts <- (width * steps) + df$start
		starts[starts > df$end] <- NULL
	} else {
		starts <- df$end
	}
	breaks <- data.frame(
		seqnames = df$seqnames[[1]],
		start = starts,
		end = starts + 1,
		strand = df$strand[[1]],
		gene_name = df$gene_name[[1]],
		type = "arrow"
	)
	return(breaks)
}

record_overlapping <- function(annotation = gene_bodies, min.gapwidth = 1000) {
	# convert back to granges
	annotation.stash <- annotation
	annotation$strand <- "*"
	gr <- makeGRangesFromDataFrame(
		df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
	)
	# work out which ranges overlap
	collapsed <- GenomicRanges::reduce(
		x = gr, with.revmap = TRUE, min.gapwidth = min.gapwidth
	)$revmap
	idx <- seq_along(gr)
	for (i in seq_along(collapsed)) {
		mrg <- collapsed[[i]]
		for (j in seq_along(mrg)) {
			idx[[mrg[[j]]]] <- j
		}
	}
	names(x = idx) <- gr$gene_name
	return(idx)
}

reformat_annotations <- function(
	annotation = annotation.subset,
	start.pos = start.pos,
	end.pos = end.pos,
	arrow_sbreaks=2000
) {
	annotation <- annotation[annotation$type == "exon"]
	exons <- as.data.frame(x = annotation)
	annotation <- split(
		x = annotation,
		f = as.character(annotation$gene_id)
	)
	annotation <- lapply(X = annotation, FUN = as.data.frame)
	
	# add gene total start / end
	gene_bodies <- lapply(seq_along(annotation), function(i){
		df <- data.frame(
			seqnames = annotation[[i]]$seqnames[[1]],
			start = min(annotation[[i]]$start),
			end = max(annotation[[i]]$end),
			strand = annotation[[i]]$strand[[1]],
			gene_name = annotation[[i]]$gene_id[[1]],
			type = "body"
		)
		# trim any that extend beyond region
		df$start <- ifelse(
			test = df$start < start.pos,
			yes = start.pos,
			no = df$start
		)
		df$end <- ifelse(
			test = df$end > end.pos,
			yes = end.pos,
			no = df$end
		)
		breaks <- split_body(df = df,width = arrow_sbreaks)
		df <- rbind(df, breaks)
		return(df)
	})
	
	gene_bodies <- do.call(what = rbind, args = gene_bodies)
	
	# record if genes overlap
	overlap_idx <- record_overlapping(annotation = gene_bodies, min.gapwidth = 1000)
	gene_bodies$dodge <- overlap_idx[as.character(gene_bodies$gene_name)]
	exons$dodge <- overlap_idx[as.character(exons$gene_id)]
	
	label_df <- gene_bodies[gene_bodies$type == "body", ]
	label_df$width <- label_df$end - label_df$start
	label_df$position <- label_df$start + (label_df$width / 2)
	
	onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
	onminus <- gene_bodies[gene_bodies$strand == "-", ]
	
	return(
		list(
			"labels" = label_df,
			"exons" = exons,
			"plus" = onplus,
			"minus" = onminus
		)
	)
}



myGenePlot <- function(annotation=tx, 
					   region="chr7:152042291-152056148",showlabel=T,
					   arrow_sbreaks=2000,font_size=18,label_size=8) {
	
	if (is.null(x = annotation)) {
		return(NULL)
	}
	if (!inherits(x = region, what = "GRanges")) {
		region <- as(region,"GRanges")
	}
	start.pos <- start(x = region)
	end.pos <- end(x = region)
	chromosome <- seqnames(x = region)
	
	# get names of genes that overlap region, then subset to include only those
	# genes. This avoids truncating the gene if it runs outside the region
	annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
	genes.keep <- unique(x = annotation.subset$gene_id)
	annotation.subset <- annotation[
		fastmatch::fmatch(x = annotation$gene_id, table = genes.keep, nomatch = 0L) > 0L
	]
	
	if (length(x = annotation.subset) == 0) {
		# make empty plot
		p <- ggplot(data = data.frame())
		y_limit <- c(0, 1)
	} else {
		annotation_df_list <- reformat_annotations(
			annotation = annotation.subset,
			start.pos = start.pos,
			end.pos = end.pos,
			arrow_sbreaks=arrow_sbreaks
		)
		p <- ggplot() +
			# exons
			geom_segment(
				data = annotation_df_list$exons,
				mapping = aes_string(
					x = "start",
					y = annotation_df_list$exons$dodge,
					xend = "end",
					yend = annotation_df_list$exons$dodge,
					color = "strand"
				),
				show.legend = FALSE,
				size = 5
			) +
			# gene body
			geom_segment(
				data = annotation_df_list$labels,
				mapping = aes_string(
					x = "start",
					y = annotation_df_list$labels$dodge,
					xend = "end",
					yend = annotation_df_list$labels$dodge,
					color = "strand"
				),
				show.legend = FALSE,
				size = 1/2
			)
		if (nrow(x = annotation_df_list$plus) > 0) {
			# forward strand arrows
			p <- p + geom_segment(
				data = annotation_df_list$plus,
				mapping = aes_string(
					x = "start",
					y = annotation_df_list$plus$dodge,
					xend = "end",
					yend = annotation_df_list$plus$dodge,
					color = "strand"
				),
				arrow = arrow(
					ends = "last",
					type = "open",
					angle = 45,
					length = unit(x = 0.05, units = "inches")
				),
				show.legend = FALSE,
				size = 1/2
			)
		}
		if (nrow(x = annotation_df_list$minus) > 0) {
			# reverse strand arrows
			p <- p + geom_segment(
				data = annotation_df_list$minus,
				mapping = aes_string(
					x = "start",
					y = annotation_df_list$minus$dodge,
					xend = "end",
					yend = annotation_df_list$minus$dodge,
					color = "strand"
				),
				arrow = arrow(
					ends = "first",
					type = "open",
					angle = 45,
					length = unit(x = 0.05, units = "inches")
				),
				show.legend = FALSE,
				size = 1/2
			)
		}
		# label genes
		n_stack <- max(annotation_df_list$labels$dodge)
		annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + (n_stack * 0.2)
		if(showlabel){
			p <- p + geom_text(
				data = annotation_df_list$labels,
				mapping = aes_string(x = "position", y = "dodge", label = "gene_name"),
				size = label_size,fontface="bold"
			)
		}
		y_limit <- c(0.9, n_stack + (n_stack * 0.5))
	}
	p <- p +
		theme_cowplot(font_size = font_size) +
		ylab("Genes") +
		xlab(label = paste0(chromosome, " position (bp)")) +
		# scale_x_continuous(expand = c(0,0),
		#                    limits = c(start.pos,end.pos),
		#                    breaks = seq(start.pos,end.pos,length.out=5))+
		scale_x_continuous(expand = c(0,0),
						   limits = c(start.pos,end.pos))+
		ylim(y_limit) +
		theme(
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.line.y = element_blank()
		) +
		scale_color_manual(values = c("darkblue", "darkgreen"))
	return(p)
}
