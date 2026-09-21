#Part 1
#cBioPortal data extraction
if(FALSE)
{
	cbio <- cBioPortal()
	studies <- getStudies(cbio, buildReport = TRUE)
	print(head(studies))
	acc <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "acc_tcga",
		    #genePanelId = "IMPACT341",
			#      geneID = 'TP53',
		    molecularProfileIds = c("acc_tcga_linear_CNA")
		)
	print(acc)
}

#Study igc get data
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	print(head(mat))
	study <- mat[,'studyId']
	for( i in study[1])
	{
		print(i)
		file <- downloadStudy(i, ask = FALSE)
		print(file)
		file_dir <- untarStudy(file, tempdir())
		print(file_dir)
		tmp <- loadStudy(file_dir)
		print(tmp)
		if(FALSE)
		{
		cna <- tmp[['cna']]
		print(head(cna))
		saveRDS(cna, file = paste0('R_data/',i,'_cna.rds'))
		mRNA <- tmp[['mrna_illumina_microarray_zscores_ref_diploid_samples']]
		print(head(mRNA))
		saveRDS(mRNA, file = paste0('R_data/',i,'_mRNA.rds'))
		}
	}
}

#Study igc save data
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	print(head(mat))
	study <- mat[,'studyId']
	for( i in study[1])
	{
		if(TRUE)
		{
		cna <- readRDS(file = paste0('R_data/',i,'_cna.rds'))
		cna <-assays(cna)
		print(cna)
		write.table( cna , file = paste0('R_data/',i,'_cna.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		}
		mRNA <- readRDS(file = paste0('R_data/',i,'_mRNA.rds'))
		mRNA <-assays(mRNA)
		print(mRNA)
		write.table( mRNA , file = paste0('R_data/',i,'_mRNA.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
	}
}

#Study ctpac read data and save data
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	print(head(mat))
	study <- mat[,'studyId']
	for( i in study[2])
	{
		print(i)
		file <- downloadStudy(i, ask = FALSE)
		print(file)
		file_dir <- untarStudy(file, tempdir())
		print(file_dir)
		tmp <- loadStudy(file_dir)
		print(tmp)
		if(FALSE)
		{
		cna <- tmp[['cna']]
		print(head(cna))
		cna <-assays(cna)
		print(cna)
		write.table( cna , file = paste0('R_data/',i,'_cna.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		mRNA <- tmp[['mrna_seq_fpkm_zscores_ref_all_samples']]
		print(head(mRNA))
		mRNA <-assays(mRNA)
		print(mRNA)
		write.table( mRNA , file = paste0('R_data/',i,'_mRNA.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		}
	}
}


#Study mbc read data and save data
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	#print(head(mat))
	study <- mat[,'studyId']
	for( i in study[3])
	{
		print(i)
		file <- downloadStudy(i, ask = FALSE)
		print(file)
		file_dir <- untarStudy(file, tempdir())
		print(file_dir)
		tmp <- loadStudy(file_dir)
		print(tmp)
		if(FALSE)
		{
		cna <- tmp[['cna']]
		print(head(cna))
		cna <-assays(cna)
		print(cna)
		write.table( cna , file = paste0('R_data/',i,'_cna.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		mRNA <- tmp[['mrna_seq_v2_rsem_zscores_ref_all_samples']]
		print(head(mRNA))
		mRNA <-assays(mRNA)
		print(mRNA)
		write.table( mRNA , file = paste0('R_data/',i,'_mRNA.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		}
	}
}


#Study metabric read data
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	#print(head(mat))
	study <- mat[,'studyId']
	for( i in study[5])
	{
		print(i)
		file <- downloadStudy(i, ask = FALSE)
		print(file)
		file_dir <- untarStudy(file, tempdir())
		print(file_dir)
		tmp <- loadStudy(file_dir)
		print(tmp)
		if(FALSE)
		{
		cna <- tmp[['cna']]
		print(head(cna))
		saveRDS(cna, file = paste0('R_data/',i,'_cna.rds'))
		mRNA <- tmp[['mrna_illumina_microarray_zscores_ref_diploid_samples']]
		print(head(mRNA))
		saveRDS(mRNA, file = paste0('R_data/',i,'_mRNA.rds'))
		}
	}
}

#Study metabric output csv
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	print(head(mat))
	study <- mat[,'studyId']
	for( i in study[4])
	{
		if(FALSE)
		{
		cna <- readRDS(file = paste0('R_data/',i,'_cna.rds'))
		cna <-assays(cna)
		print(cna)
		write.table( cna , file = paste0('R_data/',i,'_cna.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		}
		mRNA <- readRDS(file = paste0('R_data/',i,'_mRNA.rds'))
		mRNA <-assays(mRNA)
		print(mRNA)
		write.table( mRNA , file = paste0('R_data/',i,'_mRNA.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
	}
}

#Study TCGA read and save data
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	print(head(mat))
	study <- mat[,'studyId']
	for( i in study[4])
	{
		print(i)
		file <- downloadStudy(i, ask = FALSE)
		print(file)
		file_dir <- untarStudy(file, tempdir())
		print(file_dir)
		tmp <- loadStudy(file_dir)
		print(tmp)
		if(FALSE)
		{
		cna <- tmp[['cna']]
		print(head(cna))
		saveRDS(cna, file = paste0('R_data/',i,'_cna.rds'))
		mRNA <- tmp[['mrna_seq_v2_rsem_zscores_ref_diploid_samples']]
		print(head(mRNA))
		saveRDS(mRNA, file = paste0('R_data/',i,'_mRNA.rds'))
		mRNA <- tmp[['mrna_seq_v2_rsem']]
		print(head(mRNA))
		saveRDS(mRNA, file = paste0('R_data/',i,'_mRNA_RSEM.rds'))
		}
	}
}

#Study TCGA output csv
if(FALSE)
{
	library(cBioPortalData)
	library(AnVIL)
	mat <- read.csv( paste('data/study_breast.tsv'),sep='\t',row.names = 1);
	mat <- mat[mat[,'select']==1,]
	print(head(mat))
	study <- mat[,'studyId']
	for( i in study[4])
	{
		if(FALSE)
		{
		cna <- readRDS(file = paste0('R_data/',i,'_cna.rds'))
		cna <-assays(cna)
		print(cna)
		write.table( cna , file = paste0('R_data/',i,'_cna.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		mRNA <- readRDS(file = paste0('R_data/',i,'_mRNA.rds'))
		mRNA <-assays(mRNA)
		print(mRNA)
		write.table( mRNA , file = paste0('R_data/',i,'_mRNA.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		}
		mRNA <- readRDS( file = paste0('R_data/',i,'_mRNA_RSEM.rds'))
		mRNA <-assays(mRNA)
		print(mRNA)
		write.table( mRNA , file = paste0('R_data/',i,'_mRNA_rsem.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
	}
}

#Study TCGA GDC for differential expression analysis
if(FALSE)
{
	library(TCGAbiolinks)
	library(SummarizedExperiment)
# Gene expression aligned against hg19.
query <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
			   sample.type = "Primary Tumor"
                           #platform = "Illumina", 
                           #file.type  = "results",
                           #experimental.strategy = "RNA-Seq",
                           #legacy = TRUE)
			   )
#GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query, save = FALSE)
BRCAMatrix <- assay(expdat,"unstranded") 
saveRDS(BRCAMatrix,'R_data/gdc_brca_rsem_ec.rds') 
#expdat <- GDCprepare(query = query, save = TRUE, save.filename = "brca.hg38.rda")
}

if(FALSE)
{
TCGA <- readRDS('R_data/gdc_brca_rsem_ec.rds') 
print(head(TCGA))
write.table( TCGA , file = paste0('R_data/gdc_tcga_brca_rsem.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
}

if(FALSE)
{
	library(TCGAbiolinks)
# Gene expression aligned against hg19.
query <- GDCquery(project = "TCGA-BRCA",
		  legacy = TRUE,
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq",
                           file.type = "results",
                           experimental.strategy = "RNA-Seq")
#GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query, save = TRUE, save.filename = "brca.hg19.rda")
}

#Part 2 CNA analysis
#Karyotype plot
if(FALSE)
{
	library(karyoploteR)
	library(biomaRt)
	library(GenomicRanges)
	library(BSgenome.Hsapiens.UCSC.hg19)
	Hsapiens <- BSgenome.Hsapiens.UCSC.hg19

	cna_data <- read.csv('R_data/CNA_freq_subtype.csv', check.names = FALSE)

	#Genome coordinate
	mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)
	gene_coords <- getBM(
	  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
	  filters = "hgnc_symbol",
	  values = cna_data$gene,
	  mart = mart
	)

	#create genome position
	cna_annotated <- merge(cna_data, gene_coords, by.x = "gene", by.y = "hgnc_symbol")
	cna_annotated <- subset(cna_annotated, chromosome_name %in% c(1:22, "X"))
	cna_annotated$chr <- paste0("chr", cna_annotated$chromosome_name)
	cna_annotated$midpoint <- round((cna_annotated$start_position + cna_annotated$end_position) / 2)

	#pdf name
	pdf("figure/test_karyoplot_bar_subtype_all.pdf", width = 20, height = 13)

	kp <- plotKaryotype(
	  genome = "hg19",
	  chromosomes = paste0("chr", c(1:22, "X")),
	  plot.type = 4,
	  cex = 1.5,
	  labels.plotter = NULL,
	  main = "Breast Cancer CNA Frequency by Subtype"
	)
	kpAddCytobands(kp)
	kpAddChromosomeNames(kp,cex=0.8)

	# ---- Shared Scale ----
	baseline <- 0.5
	track_height <- 0.22       # height per CNA track
	track_spacing <- 0.03      # space between tracks
	max_freq <- 0.4

	# ---- Plot Tracks ----
	tracks <- c("All", "ER+", "HER2+", "TNBC")
	patient_num <- c('3970','2474','712','632')
	for (i in seq_along(tracks)) {
	  sub <- tracks[i]

	  if (sub == "All") {
	    amp_y <- cna_annotated$All_amp_freq
	    del_y <- -cna_annotated$All_del_freq
	    amp_colors <- "red"
	    del_colors <- "blue"
	    label_data <- cna_annotated
	  } else {
	    amp_freq_col <- paste0(sub, "_amp_freq")
	    del_freq_col <- paste0(sub, "_del_freq")
	    amp_fdr_col <- paste0(sub, "_amp_FDR")
	    del_fdr_col <- paste0(sub, "_del_FDR")

	    label_data <- cna_annotated[!is.na(cna_annotated[[amp_freq_col]]) & !is.na(cna_annotated[[del_freq_col]]), ]

	    amp_y <- label_data[[amp_freq_col]]
	    del_y <- -label_data[[del_freq_col]]
	    amp_sig <- label_data[[amp_fdr_col]] < 0.0001
	    del_sig <- label_data[[del_fdr_col]] < 0.0001

	    amp_colors <- ifelse(amp_sig, "darkred", "red")
	    del_colors <- ifelse(del_sig, "black", "blue")
	  }

	  # Define vertical track boundaries

	  total_track <- track_height + track_spacing
	  r1 <- 1 - (i - 1) * total_track
	  r0 <- r1 - track_height
	  y_center <- (r0 + r1) / 2

	  # Plot amplifications
	  kpSegments(kp,
	    chr = label_data$chr,
	    x0 = label_data$midpoint,
	    x1 = label_data$midpoint,
	    y0 = y_center,
	    y1 = y_center + 0.5 * amp_y / max_freq / 4/0.25 *0.22,
	    col = amp_colors
	  )

	  # Plot deletions
	  kpSegments(kp,
	    chr = label_data$chr,
	    x0 = label_data$midpoint,
	    x1 = label_data$midpoint,
	    y0 = y_center,
	    y1 = y_center + 0.5 * del_y / max_freq /4/0.25 *0.22,
	    col = del_colors
	  )
	  #kpAbline(kp, y0 = y_center+0.5*0.15/max_freq/4, y1 = y_center+0.5*0.16/max_freq/4, lty = 2, col = "black")
	  kpAbline(kp, h = y_center+0.5*0.15/max_freq/4/0.25 *0.22, lty = 2, col = "black")

	  # Add label on the left for the track
	  kpAddLabels(kp, labels = paste0(sub,'\n','n = ',patient_num[i]), r0 = r0+0.01, r1 = r1, pos = 1, srt = 90, label.margin = 0.02, cex = 1.2)
	  kpAxis(kp,
	    side = "right",
	    #ymin = -max_freq,
	    #ymax = max_freq,
	    numticks = 5,
	    r0 = r0,
	    r1 = r1,
	    labels = c(max_freq, max_freq/2, 0, max_freq/2, max_freq)
	  )
	}

	# ---- Add Unified Y-Axis ----
	kpAddLabels(kp, side = "right", labels = "Frequency", srt = 90, pos = 1, label.margin = 0.04)

	# ---- Global Legend ----
	legend(0.25,0.15, inset = -0.02, xpd = TRUE,
	       legend = c("Amp (sig)", "Amp", "Del (sig)", "Del"),
	       #col = c("darkred", "red", "darkblue", "blue"),
	       col = c("darkred", "red", "black", "blue"),
	       lty = 1, horiz = TRUE, bty = "n", cex = 1.5)

	dev.off()
}

#plot ADC target of interest on chr1, 8, 17
if(FALSE)
{
		# Load libraries
	library(karyoploteR)
	library(GenomicRanges)
	library(biomaRt)


	fc_mat <- read.csv( paste('R_data/adc_target_log2.csv'),sep=',',row.names = 1);
	gene_list <- rownames(fc_mat)

	pdf("figure/gene_locations_chr1_8_17_hg19.pdf", width = 10, height = 6)

	mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)

	gene_coords <- getBM(
	  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
	  filters = "hgnc_symbol",
	  values = gene_list,
	  mart = mart
	)

	# Keep only standard chromosomes (1–22, X, Y)
	gene_coords <- subset(gene_coords, chromosome_name %in% c(1:22, "X", "Y"))

	cat("Detected genes and their chromosome locations (hg38):\n")
	print(gene_coords)

	gene_gr <- GRanges(
	  seqnames = paste0("chr", gene_coords$chromosome_name),
	  ranges = IRanges(
	    start = gene_coords$start_position,
	    end = gene_coords$end_position
	  ),
	  gene = gene_coords$hgnc_symbol
	)


	chrom_order <- paste0("chr", c(1:22, "X", "Y"))
	chroms_to_plot <- intersect(chrom_order, unique(seqnames(gene_gr)))

	cat("Chromosomes to plot (sorted):", paste(chroms_to_plot, collapse = ", "), "\n")

	kp <- plotKaryotype(
	  genome = "hg38",
	  chromosomes = chroms_to_plot,
	  plot.type = 2  # Side-by-side layout
	)

	kpPlotMarkers(
	  kp,
	  data = gene_gr,
	  labels = gene_gr$gene,
	  #text.orientation = "horizontal",
	  adjust.label.position = TRUE,   # avoid overlaps
	  label.margin = 0.005,
	  cex = 0.7,
	  r0 = 0, r1 = 0.9,
	  label.col = "darkred",
	  line.color = "gray30"
	)
	dev.off()
}

if(FALSE)
{
	library(karyoploteR)
	library(biomaRt)
	library(GenomicRanges)
	library(BSgenome.Hsapiens.UCSC.hg19)
	Hsapiens <- BSgenome.Hsapiens.UCSC.hg19

	# ---- Load Data ----
	cna_data <- read.csv('R_data/CNA_freq_subtype.csv', check.names = FALSE)

	# ---- Get gene coordinates using biomaRt ----
	mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)
	gene_coords <- getBM(
	  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","band"),
	  filters = "hgnc_symbol",
	  values = cna_data$gene,
	  mart = mart
	)

	# ---- Merge & Clean ----
	cna_annotated <- merge(cna_data, gene_coords, by.x = "gene", by.y = "hgnc_symbol")
	cna_annotated <- subset(cna_annotated, chromosome_name %in% c(1:22, "X"))
	cna_annotated$chr <- paste0("chr", cna_annotated$chromosome_name)
	#cna_annotated$midpoint <- round((cna_annotated$start_position + cna_annotated$end_position) / 2)
	print(head(cna_annotated))
	write.table( cna_annotated , file = paste0('R_data/cna_freq_subtype_chromosome.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
}

#part 3 differential expression
#deseq2
if(FALSE)
{
	mat <- read.csv( paste('R_data/normal_tissue_gdc.tsv'),sep='\t',row.names = 1);
	saveRDS(mat, 'R_data/normal_tissue_gdc.rds')
}

#differential expression of each tissue
if(FALSE)
{
library(DESeq2)
library(readr)

expr <-	readRDS('R_data/normal_tissue_gdc.rds')
sample_info <- read.csv("R_data/normal_tissue_info_gdc.tsv", sep = '\t', row.names = 1)  # Or use readRDS if already serialized
#rownames(sample_info) <- sample_info$sample_id
# Set tissue as a factor with reference = "breast_cancer"
sample_info$tissue <- factor(sample_info$tissue)
sample_info$tissue <- relevel(sample_info$tissue, ref = "Breast_Cancer")


	# Get all tissue types except breast_cancer
tissue_levels <- unique(sample_info$tissue)
tissue_levels <- setdiff(tissue_levels, "Breast_Cancer")

# Loop over tissues, compare each to breast_cancer
for (tissue in tissue_levels) {
  cat("🔄 Running DE for:", tissue, "vs breast_cancer\n")

  # Subset samples for this comparison
  samples_subset <- sample_info$tissue %in% c("Breast_Cancer", tissue)
  sub_sample_info <- sample_info[samples_subset, ]
  sub_expr <- expr[, rownames(sub_sample_info)]

  # Set tissue as factor with breast_cancer as reference
  sub_sample_info$tissue <- factor(sub_sample_info$tissue, levels = c("Breast_Cancer", tissue))

  # Construct DESeq2 object
  dds_sub <- DESeqDataSetFromMatrix(countData = sub_expr,
                                    colData = sub_sample_info,
                                    design = ~ tissue)

  # Filter low-count genes
  keep <- rowSums(counts(dds_sub)) >= 10
  dds_sub <- dds_sub[keep, ]

  # Run DESeq only on this subset
  dds_sub <- DESeq(dds_sub)

  # Get results: tissue vs breast_cancer
  res <- results(dds_sub, contrast = c("tissue", tissue, "Breast_Cancer"))
  res <- res[order(res$padj), ]

  # Save results
  out_file <- paste0("R_data/deseq/DE_", tissue, "_vs_Breast_Cancer_gdc.csv")
  write.csv(as.data.frame(res), file = out_file)

  cat("✅ Saved:", out_file, "\n")
}
}

if(FALSE)
{
library(ggplot2)
library(dplyr)

expr_matrix <-	readRDS('R_data/normal_tissue.rds')
sample_info <- read.csv("R_data/normal_tissue_info.tsv", sep = '\t', row.names = 1)  # Or use readRDS if already serialized
# Step 1: Transpose expression matrix to sample x gene (for merging)
expr_t <- t(expr_matrix)
expr_df <- as.data.frame(expr_t)

# Step 2: Add sample IDs if not already rownames
expr_df$sample_id <- rownames(expr_df)

# Step 3: Merge with sample metadata
sample_info$sample_id <- rownames(sample_info)
merged_df <- left_join(expr_df, sample_info, by = "sample_id")

# Step 4: Plot ERBB2 expression by tissue
	pdf("figure/test_ERBB2.pdf", onefile = TRUE, width = 25, height = 13)
p <- ggplot(merged_df, aes(x = tissue, y = ERBB2)) +
  geom_boxplot(aes(fill = tissue)) +
  coord_cartesian(ylim = c(0, 50000)) +  # limits y-axis display
  theme_bw() +
  labs(
    title = "ERBB2 Expression Across Tissues",
    x = "Tissue Type",
    y = "ERBB2 Expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()
}

#differential expression and average expression plot
if(FALSE)
{
# Load necessary libraries
library(ggplot2)
library(ggrepel) # For better label placement

# Read the data generated by Python
data <- read.csv("R_Data/gene_data_for_r.txt", sep = "\t", row.names=1)

data$gene_symbol <- rownames(data)
print(head(data))
# Ensure Finalist is treated as a factor for mapping
data$is_finalist <- factor(data$Finalist)

# Create the ggplot
p <- ggplot(data, aes(y = mean, x = avg_log2FoldChange, color = annotation)) +
  # Map size to Finalist status: Finalists (1) are larger than others (0)
  #geom_point(aes(size = is_finalist), alpha = 0.7) +
  geom_point(aes(alpha = is_finalist, size = is_finalist)) +
  scale_alpha_manual(values = c("0" = 0.2, "1" = 1.0), guide = "none") +
  # Manually control sizes (e.g., 4 for finalists, 1.5 for others)
  scale_size_manual(values = c("0" = 5, "1" = 10), guide = "none") +
        scale_x_continuous(limits = c(-5,7))+
  
  # Label gene names only for Finalist == 1
  geom_text_repel(
    aes(label = ifelse(Finalist == 1, as.character(gene_symbol), "")),
    show.legend = FALSE,
    color = "black",
    size = 10,
    max.overlaps = 50,
    box.padding = 0.5
  ) +
  
  # Aesthetic improvements
  labs(
    title = "Gene Expression vs Log Fold Change",
    y = "Mean Expression",
    x = "Average Log2 Fold Change",
    color = "Subtype Annotation"
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
   theme_classic(base_size = 35)+
   theme(legend.position = c(0.2,0.9))

# Display the plot
pdf("figure/exp_avg_log_latest.pdf", onefile = TRUE, width = 25, height = 17)
print(p)
dev.off()

# Save the plot
# ggsave("gene_expression_plot.png", plot = p, width = 10, height = 8)
}

#Tissue differential expression specific plot.
#complex heatmap to show log2 change
if(FALSE)
{
library(ComplexHeatmap)
library(circlize)

set.seed(123)

### 1. Generate random data

# Fold change matrix
fc_mat <- read.csv( paste('R_data/adc_target_log2.csv'),sep=',',row.names = 1);

# FDR matrix
fdr_mat <- read.csv( paste('R_data/adc_target_p.csv'),sep=',',row.names = 1);

# Significance matrix (TRUE if FDR < 0.05)
sig_mat <- fdr_mat < 0.0001

### 2. Color scale for fold change
col_fun <- colorRamp2(
    c(-5, 0, +5),
    c("blue", "white", "red")
)

### 3. Create a heatmap annotation layer for significance
# Mark significant cells with a star
sig_anno <- function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j]) {
        grid.text("*", x = x, y = y, gp = gpar(fontsize = 12, col = "black"))
    }
}

### 4. Build the heatmap
ht <- Heatmap(
    fc_mat,
    name = "Fold Change",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    cell_fun = sig_anno,
    row_title = "Genes",
    column_title = "Tissues",
    heatmap_legend_param = list(title = "Fold Change")
)

### 5. Draw heatmap
pdf("figure/Fig_tissue.pdf", onefile = TRUE, width = 25, height = 17)
print(ht)
dev.off()
}
