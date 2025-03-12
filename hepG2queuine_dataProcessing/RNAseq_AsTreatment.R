library(tidyverse)
basepath <- # base path
fc_cutoff <- 2
pval_cutoff <- 0.05
red <- hcl.colors(hcl.pals("sequential")[8], n = 5)[4] #|> scales::show_col()
blue <- hcl.colors(hcl.pals("sequential")[4], n = 5)[4] #|> scales::show_col()
x <- seq(0, 1, length.out = 4)
cp <- scales::gradient_n_pal(c("#9C5D41", "#97928A", "#C1ABAD", "#D1D4D0"), x)(x)


# 1. overlap majiq and leafcutter output for GSE205332 -----
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

# 2. gene ontology analysis ------
shortList <- readRDS(file.path(basepath, 'rds/overlap_majiq_leafcutter'))
allquant <- read_tsv(file.path(basepath, 'majiq/3vs3_alltbl.tsv'), skip = 10)

sum(shortList$gene_id %in% allquant$gene_id)

shortList <- str_extract(shortList$gene_id, pattern = "^.+(?=\\.)") |> unique()
fullList <- str_extract(allquant$gene_id, pattern = "^.+(?=\\.)") |> unique()
sum(shortList %in% fullList)

bp <- clusterProfiler::groupGO(gene = shortList, keyType = 'ENSEMBL', OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "BP", level = 6, readable = TRUE); bp@result <- bp@result %>% arrange(desc(Count))
enrichplot:::barplot.enrichResult(bp, font.size = 10, showCategory=10, label_format = 50)


mf <- clusterProfiler::groupGO(gene = shortList, keyType = 'ENSEMBL', OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "MF", level = 5, readable = TRUE); mf@result <- mf@result %>% arrange(desc(Count))
enrichplot:::barplot.enrichResult(mf, font.size = 10, showCategory=10, label_format = 50)+ scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10))
mf@result$geneID[1]

cc <- clusterProfiler::groupGO(gene = shortList, keyType = 'ENSEMBL', OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "CC", level = 3, readable = TRUE); cc@result <- cc@result %>% arrange(desc(Count))
enrichplot:::barplot.enrichResult(cc, font.size = 10, showCategory=10, label_format = 50)


View(bp4@result)
View(mf4@result)
View(cc4@result)
enrichplot:::barplot.enrichResult(bp4, showCategory=20)

# ego3 <- clusterProfiler::gseGO(geneList = sort(shortList), keyType = 'ENSEMBL', OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "ALL", minGSSize = 100,maxGSSize = 500,pvalueCutoff = 0.05, verbose = FALSE)

# cp <- clusterProfiler::enrichGO(gene = shortList,
#                           universe = fullList,
#                           OrgDb = org.Hs.eg.db::org.Hs.eg.db,
#                           keyType = 'ENSEMBL',
#                           ont           = "ALL",
#                           pAdjustMethod = "BH",
#                           pvalueCutoff  = 0.05,
#                           readable      = TRUE)

AnnotationDbi::keytypes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)


# 3. spliceosome (https://www.genenames.org/data/genegroup/#!/group/1518)----
spliceosomeProteins <- read_delim('https://www.genenames.org/cgi-bin/genegroup/download?id=1518&type=branch') 
hs <- UniProt.ws::UniProt.ws(9606)
t <- UniProt.ws::select(x = hs, keys = spliceosomeProteins$`Ensembl gene ID`, columns = c('ft_zn_fing'), keytype = 'Ensembl')
hgnc <- UniProt.ws::select(x = hs, keys = spliceosomeProteins$`HGNC ID`, columns = c('ft_zn_fing', 'gene_primary'), keytype = 'HGNC')
spliceosomeProteins$`HGNC ID` |> unique()
hgnc$From |> unique()
hgnc$Entry |> unique()
tb <- hgnc %>% 
  filter(!is.na(Zinc.finger)) 
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', tb); openxlsx::openXL(wb)
# AnnotationDbi::keytypes(hs)
View(spliceosomeProteins)


# 4. extract bam reads that overlapped with  -----
bamFile <- file.path('/Volumes/sgms/sequencing/hepG2_quieuine/bam', list.files('/Volumes/sgms/sequencing/hepG2_quieuine/bam', pattern = '*(1|2|3)_Homo_sapiens_RNA-Seq_Aligned.sortedByCoord.out.bam$'))
bamtbl <- tibble::tibble(path = bamFile, samp = paste0(rep(c('trt', 'ctrl'), 3), rep(3:1, each = 2))) %>% 
  dplyr::arrange(samp)

igvWeb <- function(gn, region){
  geneShortListed <- readxl::read_excel(file.path(basepath, 'joinMajiq_Leafcutter_rosaMacriteria.xlsx'), sheet = 'S7.MAJIQ_LeafCutter_OL') %>%  # readRDS(file.path(basepath, 'rds/overlap_majiq_leafcutter')) %>% 
    dplyr::filter(gene_name %in% gn) %>% 
    tidyr::separate_rows(exons_coords, sep = ';', convert = TRUE) %>% 
    tidyr::separate_rows(exons_coords, sep = '-', convert = TRUE) %>%
    dplyr::group_by(gene_name, seqid) %>% 
    dplyr::summarise(start = min(exons_coords, na.rm = TRUE),
                     end = max(exons_coords, na.rm = TRUE), 
                     strand = unique(strand)) %>% 
    dplyr::ungroup()
  gsl <- GenomicRanges::makeGRangesFromDataFrame(df = geneShortListed, seqnames.field = 'seqid', start.field = 'start', end.field = 'end',strand.field = 'strand')
  names(gsl) <- gn
  gsl_paste <- with(geneShortListed, sprintf("%s:%d-%d",seqid, start, end))
  
  igv <- igvR::igvR()
  BrowserViz::setBrowserWindowTitle(igv, gn)
  igvR::setGenome(igv, "hg38")
  
  igvR::showGenomicRegion(igv, gsl_paste)
  loc <- igvR::getGenomicRegion(igv)
  
  for (i in seq_len(nrow(bamtbl))){
    ga <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bamtbl$path[i]), use.names = TRUE, param = Rsamtools::ScanBamParam(which = GenomicRanges::makeGRangesFromDataFrame(loc)))
    track <- igvR::GenomicAlignmentTrack(trackName = bamtbl$samp[i], alignment = ga, visibilityWindow = 40000)
    igvR::displayTrack(igv, track)
  }
}

igvWeb('SEMA3F')



# genome to RNA ----
rearrangeRanges <- function(ir, wheretostart=IRanges::start(ir)[1], myspacer = rep(0, length(ir))) {
  # Dr. Girke's splice product function
  library(IRanges)
  mywidth <- width(ir)
  mywidth[1] <- mywidth[1] + (wheretostart - 1)
  mywidth <- c(mywidth[1], (mywidth[-1] + myspacer[-1]))
  myend <- cumsum(mywidth)
  mystart <- myend - width(ir) + 1
  ir_shift <- IRanges::IRanges(start=mystart, end=myend, names=names(ir))
  return(ir_shift)
}
genomeToRNA <- function(gn, writeSeq = FALSE){
  grs <- readxl::read_excel(file.path(basepath, 'joinMajiq_Leafcutter_rosaMacriteria.xlsx'), sheet = 'S7.MAJIQ_LeafCutter_OL') %>% 
    dplyr::filter(gene_name %in% gn) %>%
    tidyr::separate_rows(exons_coords, sep = ';', convert = TRUE) %>% 
    separate(exons_coords, into = c('start', 'end'), sep = '-', convert = TRUE) %>% 
    mutate(strand = '+') %>% 
    dplyr::rename(chr = seqid)
  # selectedExonCoord <- GenomicRanges::makeGRangesFromDataFrame(df = grs, seqnames.field = 'chr', start.field = 'start', end.field = 'end',strand.field = 'strand')
  grl <- c()
  grl <- GenomicRanges::makeGRangesListFromDataFrame(df = grs, split.field = "gene_name")
  bseq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, grl)
  
  if (writeSeq){
    finalSeq <- bseq[[1]][1]
    for (i in 2:length(bseq[[1]])){
      finalSeq <- Biostrings::xscat(finalSeq, bseq[[1]][i])
    }
    # BSgenome::width(finalSeq)
    Biostrings::writeXStringSet(x = finalSeq, file.path(basepath, 'fasta_regionOfInterest_2ndBatch', paste(gn, 'transctipt', sep = '_')))
  } else {
    print(gn)
    rearrangeRanges(ir =  IRanges::ranges(grl[[1]]))
  }
}


genomeToRNA('SEMA3F', writeSeq = TRUE)
