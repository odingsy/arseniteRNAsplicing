library(tidyverse)
basepath <- '' # base path
fc_cutoff <- 2
pval_cutoff <- 0.05
red <- hcl.colors(hcl.pals("sequential")[8], n = 5)[4] #|> scales::show_col()
blue <- hcl.colors(hcl.pals("sequential")[4], n = 5)[4] #|> scales::show_col()
x <- seq(0, 1, length.out = 4)
cp <- scales::gradient_n_pal(c("#9C5D41", "#97928A", "#C1ABAD", "#D1D4D0"), x)(x)


# Table S7 overlap majiq and leafcutter output for GSE205332 -----
# rosa ma criteria: both MAJIQ (P(ΔΨ > 0.1) > 0.95) and LeafCutter (P < 0.05)
lc <- inner_join(read_delim(file.path(basepath, 'lc/3v3_effect_sizes.txt')) %>%
             separate(intron, into = c('chr', 'start', 'end', 'cid'), sep = ':',remove = TRUE, convert = TRUE) %>%
             mutate(strand = str_extract(cid, pattern = "(\\+|-)"),
                    cluster = paste(chr, cid, sep = ':')),
           read_delim(file.path(basepath, 'lc/3v3_cluster_significance.txt')) %>%
             filter(p.adjust < 0.05), by = c('cluster')) %>% 
# leafcutter_threshold	leafcuter_sites 	overlap_with_majiq_(61sites) 	SEMA3F_identified?
# adj p-value < 0.05  772	41	yes 
# adj p-value < 0.01	293	29	yes 
# adj p-value < 0.001	70	8	no
  mutate(width = end - start + 1) %>%
  group_by(cid) %>%
  slice(which.max(width)) %>%
  # filter(any(abs(deltapsi) > 0.1)) %>%
  # dplyr::slice_max(order_by = width, n = 1, with_ties = FALSE) %>%
  ungroup()
lc_gr <- GenomicRanges::makeGRangesFromDataFrame(df = lc, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')

env <- environment()
mj <- read_tsv(file.path(basepath, 'majiq/3vs3_filtedtbl.tsv'), skip = 10) %T>% assign('mj_ori',., envir = env) %>%
  separate_rows(probability_changing, junctions_coords, sep = ';', convert = TRUE) %>%
  separate(junctions_coords, into = c('start', 'end'), sep = '-', convert = TRUE) %>%
  mutate(width = end - start + 1) %>%
  dplyr::group_by(lsv_id) %>%
  slice(which.max(width)) %>%
  # dplyr::slice_max(order_by = width, n = 1, with_ties = FALSE) %>%
  ungroup()
mj_gr <- GenomicRanges::makeGRangesFromDataFrame(df = mj, seqnames.field = 'seqid', start.field = 'start', end.field = 'end', strand.field = 'strand')

ol <- GenomicRanges::findOverlaps(mj_gr, lc_gr)
setdiff(IRanges::from(ol), which(mj$gene_name %in% unique(lc$genes)))

ij <- inner_join(mj[IRanges::from(ol),] %>%
                   select(lsv_id),
  mj_ori, by = 'lsv_id')

if (!dir.exists(file.path(basepath, 'rds'))) dir.create(file.path(basepath, 'rds'))
# saveRDS(object = ij, file = file.path(basepath, 'rds/overlap_majiq_leafcutter'))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, 'S7.MAJIQ_LeafCutter_OL')
lsty1 <- openxlsx::createStyle(bgFill = "#d6d2d2")
openxlsx::writeData(wb, 'S7.MAJIQ_LeafCutter_OL', ij)
openxlsx::conditionalFormatting(wb, 'S7.MAJIQ_LeafCutter_OL', rows = 1:(nrow(ij)+1), cols = 1:ncol(ij), rule = paste0("$", LETTERS[which(colnames(ij) == "shade")], "1==", 1), style = lsty1)
openxlsx::openXL(wb)

# Comparing with Riedmann splicing list  -----
# https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-015-1295-9/MediaObjects/12864_2015_1295_MOESM12_ESM.pdf
drs <- scan(text = 'UCA1 MFAP5 ARHGDIB MGP
CD36 LPL SAGE1 UST
GJB5 SLC43A3 NCAM2 COL5A2
AKT3 ITPR1 GNAT3 CYP4F22
ITGB4 PLAUR KCTD12 SCN9A
COL15A1 LGALS3BP SEMA6A NPAS2
EDNRA L1CAM ETV4 MRPL17
SEMA4B ALOX15B RNF122 FRMD3
HLA-DMA LINGO2 ALS2CR8 SDC2
TTLL1 SIDT1 RIBC2 COL25A1
SASH1 ARL6IP6 ISM2 IDH2
SUSD4 C16orf58 STAT5A DIRA3
ANK2 CPM LUZP2 DOTL1L
PLXND1 RPS6KA2 TSPAN18 SHC4
ELANE MDGA2 KCNV1 CDY2B
HLA-F PTH1R FHL2 P2RX4
GATA4 IFT140 KLHL31 SHC3
UQCRC1 NUDT8 CDH5 NPFFR2
SLITRK4 PDZD4 MKX FCER1G
PTK6 PRKCG SEC22C NAT6
NEGR1-IT1', what = character(), quiet = TRUE)

urs <- scan(text = 'DOK2 THAP3 PDCD6 NOSTRIN
PRR23A PYCARD GZMH STEAP2
T1CAM1 SLC45A4 COL24A1 C1S
FAM212B RYR2 LIMD1 SQSTM1
ZNF469 FABP3 TRIM16 ZNF323
DHRS2 UGAT1A5 SEMA3A ALDH3A1
SH3GL2 ABCG2 HMOX1', what = character(), quiet = TRUE)

drs_gn <- AnnotationDbi::select(Homo.sapiens::Homo.sapiens, keys = drs, columns = c('ENSEMBL'), keytype = 'SYMBOL')
drs_gn$ENSEMBL[drs_gn$SYMBOL == 'NAT6'] <- 'ENSG00000243477';drs_gn$ENSEMBL[drs_gn$SYMBOL == 'DOTL1L'] <- 'ENSG00000104885'; drs_gn$ENSEMBL[drs_gn$SYMBOL == 'DIRA3'] <- 'ENSG00000162595'; drs_gn$ENSEMBL[drs_gn$SYMBOL == 'C16orf58'] <- 'ENSG00000140688';drs_gn$ENSEMBL[drs_gn$SYMBOL == 'ALS2CR8'] <- 'ENSG00000138380'
urs_gn <- AnnotationDbi::select(Homo.sapiens::Homo.sapiens, keys = urs, columns = c('ENSEMBL'), keytype = 'SYMBOL')
urs_gn <- urs_gn[urs_gn$SYMBOL != 'T1CAM1', ]; urs_gn <- urs_gn[urs_gn$SYMBOL != 'UGAT1A5', ]; urs_gn$ENSEMBL[urs_gn$SYMBOL == 'FAM212B'] <- 'ENSG00000197852'; urs_gn$ENSEMBL[urs_gn$SYMBOL == 'ZNF323'] <- 'ENSG00000235109'
sum(str_extract(ij$gene_id, pattern = '^[^\\.]*') %in% drs_gn$ENSEMBL)
sum(str_extract(ij$gene_id, pattern = '^[^\\.]*') %in% urs_gn$ENSEMBL)
sum(str_extract(mj$gene_id, pattern = '^[^\\.]*') %in% drs_gn$ENSEMBL)
sum(str_extract(mj$gene_id, pattern = '^[^\\.]*') %in% urs_gn$ENSEMBL)


