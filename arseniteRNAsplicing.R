library(tidyverse)
hs <- UniProt.ws::UniProt.ws(9606)
sepCollapse <- '||' # separator for concatenation annotation information. 

# load the data ----
base <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive'
raw <- file.path(base, 'PhD_study/conferences/23ASMS/rds') 
# split_pp <- readr::read_delim(file.path(base, 'PhD_study/arenite_competitive/data/combinedRedoSearch_Rep1Rep2/230329_AsCompetitive_combinedSplitSearch', 'peptides.txt'), na = c("", "NA","#NUM!"))
split_pg <- readr::read_delim(file.path(base, 'PhD_study/arenite_competitive/data/combinedRedoSearch_Rep1Rep2/230329_AsCompetitive_combinedSplitSearch', 'proteinGroups.txt'), na = c("", "NA","#NUM!")) %>%
  filter(!(grepl("(CON|REV)__", `Protein IDs`) | grepl("(CON|REV)__", `Majority protein IDs`)))

# preprocessing  -----
nuctbl <- split_pg %>%
  select(`Protein IDs`, `Protein names`, matches('Ratio H/L normalized nuc(F|R)')) %>%
  # mutate(across(where(is.double), ~na_if(., 'NaN'))) %>% # replace the 'NaN' (literal, need to be quoted) in the ratio col with NA
  pivot_longer(matches('Ratio H/L normalized.+$'), names_to = "samp", values_to = "ratio_normalized") %>%
  mutate(samp = str_extract(samp, "nuc(F|R)_\\d+"),
         rep = str_extract(samp, 'F|R'),
         ratio_normalized = if_else(str_detect(samp, 'F'), 1/ratio_normalized, ratio_normalized)) %>% # 1/forward = ctrl(L)/competitive(H); reverse = ctrl(H)/competitive(L)
  # group_by(`Protein IDs`) %>% filter(any(ratio_normalized> 1.5)) %>% # ungroup(rep) %>% filter(any(rep == 'F') & any(rep == 'R')) %>%
  # ungroup() %>% # filter at least one F and R > 1.5
  select(-rep) %>%
  nest(data = -c('Protein IDs', 'Protein names')) %>%
  mutate(mean = map_dbl(data, ~{mean(.$ratio_normalized, na.rm = TRUE) %>% as.numeric(.)}),
         rsd = map_dbl(data, ~{sd(.$ratio_normalized, na.rm = TRUE) %>% as.numeric(.)}) / mean,
         `num<1.5` = map_dbl(data, ~{ sum(c(is.na(.$ratio_normalized), sum(.$ratio_normalized < 1.5, na.rm = TRUE))) %>% as.numeric(.)}),
         numOfNA = map_dbl(data, ~{is.na(.$ratio_normalized) %>% sum(.) %>% as.numeric(.)})) %>%
  # filter(`num<1.5` <= 2 & (mean > 2.5 | (mean > 1.5 & rsd < 0.2))) %>%
  unnest(data) %>%
  pivot_wider(names_from = "samp", values_from = "ratio_normalized") %>%
  relocate(`Protein IDs`, `Protein names`, matches('nuc(F|R)_1'), matches('nuc(F|R)_2'))

# obtain total number of quantified protein 
nuctbl %>% 
  filter(numOfNA < 4) %>% 
  nrow()

# obtain annotation information from uniprot
upid <- nuctbl %>% separate_rows(`Protein IDs`) %>% pull(`Protein IDs`)
geneSymbol <- UniProt.ws::select(x = hs, keys =  upid, columns = c('gene_primary'), keytype = 'UniProtKB')
zf <- UniProt.ws::select(x = hs, keys = upid, columns = c('ft_zn_fing'), keytype = 'UniProtKB')
nuctbl <- nuctbl %>%
  rowwise() %>%
  mutate(gn = geneSymbol$Gene.Names..primary.[geneSymbol$From %in% unlist(strsplit(`Protein IDs`, split = ';'))] %>% unique(.) %>% .[complete.cases(.)] %>% paste0(., collapse = sepCollapse),
         zincFinger = zf$Zinc.finger[zf$From %in% unlist(strsplit(`Protein IDs`, split = ';'))] %>% unique(.) %>% .[complete.cases(.)] %>% paste0(., collapse = sepCollapse)) %>%
  ungroup()

# a table summarizing all the interactive candidates ----
tblout <- nuctbl %>% 
  filter(nchar(zincFinger) !=0 & `num<1.5` <= 2 & mean > 1.5) %>% 
  mutate(across(where(is.double), ~signif(.x, digits = 3))) %>% 
  arrange(desc(mean)) %>%
  mutate(gn = factor(gn, levels = gn)) %>%
  mutate(rnum = row_number(),
         shade = case_when(rnum %% 2 == 1 ~ 0, rnum %% 2 == 0 ~ 1)) %>%
  select(-rnum) %>% 
  rename(Gene = gn) %>% 
  select(Gene,`Protein IDs`, starts_with('nuc'), mean, rsd, zincFinger, everything())
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, 'test')
lsty1 <- openxlsx::createStyle(bgFill = "#d6d2d2")
openxlsx::writeData(wb, 'test', tblout)
openxlsx::conditionalFormatting(wb, 'test', rows = 1:(nrow(tblout)+1), cols = 1:ncol(tblout), rule = paste0("$", LETTERS[which(colnames(tblout) == "shade")], "1==", 1), style = lsty1)
openxlsx::openXL(wb)

# GO analysis of the shortlist protein ----
# AnnotationDbi::columns(org.Hs.eg.db)
gene.df <- clusterProfiler::bitr(tblout$`Protein IDs`, fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Hs.eg.db::org.Hs.eg.db)
bp7 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "BP", level = 7, readable = TRUE); bp7@result <- bp7@result %>% arrange(desc(Count))
cc6 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "CC", level = 6, readable = TRUE); cc6@result <- cc6@result %>% arrange(desc(Count))
mf6 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "MF", level = 6, readable = TRUE); mf6@result <- mf6@result %>% arrange(desc(Count))

enrichplot:::barplot.enrichResult(bp7, font.size = 10, showCategory=10, label_format = 50)
enrichplot:::barplot.enrichResult(cc6, font.size = 10, showCategory=10, label_format = 50)
enrichplot:::barplot.enrichResult(mf6, font.size = 10, showCategory=10, label_format = 50)

wb <- openxlsx::createWorkbook(); 
openxlsx::addWorksheet(wb, 'BiologicalProcess'); openxlsx::writeData(wb, 'BiologicalProcess', head(bp7, 20)); 
openxlsx::addWorksheet(wb, 'CellularComponent'); openxlsx::writeData(wb, 'CellularComponent', head(cc6, 20)); 
openxlsx::addWorksheet(wb, 'MolecularFunction'); openxlsx::writeData(wb, 'MolecularFunction', head(mf6, 20)); 
openxlsx::openXL(wb)


# a dotplot shows the Forward and Reverse from the first replicate --------
ftbl <- nuctbl %>% 
  mutate(zf_type = case_when(
    stringr::str_detect(zincFinger, 'UBZ4') ~ 'UBZ4', # type 4 UBZs are CCHC http://www.ebi.ac.uk/interpro/entry/profile/PS51908/
    stringr::str_detect(zincFinger, 'PHD') ~ 'PHD',
    stringr::str_detect(zincFinger, 'C2H2') ~ 'C2H2',
    stringr::str_detect(zincFinger, 'C6H2') ~ 'C6H2',
    stringr::str_detect(zincFinger, 'C3H1') ~ 'C3H1',
    stringr::str_detect(zincFinger, 'CCHC') ~ 'CCHC',
    stringr::str_detect(zincFinger, 'C4') ~ 'C4',
    TRUE ~ ''
  )) %>% 
  mutate(`Arsenite replaceable zinc-finger domain` = case_when(
    zincFinger != '' ~ 'prefered',
    zincFinger %in% c('C2H2') ~ 'not prefered',
    zincFinger == '' ~ ''
  )) %>% 
  mutate(nucF_1 = ifelse(gn == 'POLD1', nucF_2, nucF_1), 
         nucR_1 = ifelse(gn == 'POLD1', nucR_2, nucR_1),
         nucF_1 = ifelse(gn == 'METAP1', nucF_2, nucF_1), 
         nucR_1 = ifelse(gn == 'METAP1', nucR_2, nucR_1)) %>% 
  filter(gn != 'RAD18')

p <- ftbl %>%  
  filter((!is.na(nucF_1)) & (!is.na(nucR_1))) %>% 
  ggplot(aes(x = nucF_1, y = nucR_1))+
  geom_point(color = 'grey', size = 0.5, alpha = .5)+
  geom_point(data = . %>% filter(nchar(zincFinger) !=0 & `num<1.5` <= 2 & mean > 1.5), aes(color = `Arsenite replaceable zinc-finger domain`), size = 0.5)+
  ggrepel::geom_label_repel(data = . %>% filter(gn %in% c('SF1', 'METAP1', 'ZC3H15', 'POLD1')), aes(label = gn), label.size = NA)+
  geom_point(data = . %>% filter(nchar(zincFinger) !=0 & `num<1.5` <= 2 & mean > 1.5) %>% filter(gn %in% c('POLD1')), aes(color = `Arsenite replaceable zinc-finger domain`), size = 1)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0,6), ylim = c(0,6), xintercept = 1, yintercept = 1, ratio = 1)+
  labs(x = "1/forward = ctrl(L)/competitive(H)", y = "reverse = ctrl(H)/competitive(L)")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none")

ggsave(filename = file.path(base,'PhD_study/conferences/23ASMS/pngs/protein_dotplot.svg'), device = 'svg', plot = p, width = 5, height = 5, units = 'in')


# overlap majiq and leafcutter output for GSE205332,  MAJIQ (P(ΔΨ > 0.1) > 0.95) and LeafCutter (P < 0.05) -----
env <- environment()
mj <- readxl::read_excel(file.path(base, 'PhD_study/conferences/23ASMS','SI_tables.xlsx'), sheet = 'S5.MAJIQ_outputTable', skip = 10) %T>% assign('mj_ori',., envir = env) %>%
  separate_rows(probability_changing, junctions_coords, sep = ';', convert = TRUE) %>%
  separate(junctions_coords, into = c('start', 'end'), sep = '-', convert = TRUE) %>%
  mutate(width = end - start + 1) %>%
  dplyr::group_by(lsv_id) %>%
  slice(which.max(width)) %>%
  # dplyr::slice_max(order_by = width, n = 1, with_ties = FALSE) %>%
  ungroup()
mj_gr <- GenomicRanges::makeGRangesFromDataFrame(df = mj, seqnames.field = 'seqid', start.field = 'start', end.field = 'end', strand.field = 'strand')

lc <- readxl::read_excel(file.path(base, 'PhD_study/conferences/23ASMS','SI_tables.xlsx'), sheet = 'S6.LeafCutter_outputTable')
# lc <- inner_join(read_delim(file.path(basepath, 'lc/3v3_effect_sizes.txt')) %>%
#                    separate(intron, into = c('chr', 'start', 'end', 'cid'), sep = ':',remove = TRUE, convert = TRUE) %>%
#                    mutate(strand = str_extract(cid, pattern = "(\\+|-)"),
#                           cluster = paste(chr, cid, sep = ':')),
#                  read_delim(file.path(basepath, 'lc/3v3_cluster_significance.txt')) %>%
#                    filter(p.adjust < 0.05), by = c('cluster')) %>%
#   mutate(width = end - start + 1) %>%
#   group_by(cid) %>%
#   slice(which.max(width)) %>%
#   # filter(any(abs(deltapsi) > 0.1)) %>%
#   # dplyr::slice_max(order_by = width, n = 1, with_ties = FALSE) %>%
#   ungroup()
lc_gr <- GenomicRanges::makeGRangesFromDataFrame(df = lc, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')

ol <- GenomicRanges::findOverlaps(mj_gr, lc_gr)
ij <- inner_join(mj[IRanges::from(ol),] %>%
                   select(lsv_id),
                 mj_ori, by = 'lsv_id') %>% 
  mutate(rnum = row_number(),
         shade = case_when(rnum %% 2 == 1 ~ 0, rnum %% 2 == 0 ~ 1)) %>%
  select(-rnum)
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, 'S7.MAJIQ_LeafCutter_OL')
lsty1 <- openxlsx::createStyle(bgFill = "#d6d2d2")
openxlsx::writeData(wb, 'S7.MAJIQ_LeafCutter_OL', ij)
openxlsx::conditionalFormatting(wb, 'S7.MAJIQ_LeafCutter_OL', rows = 1:(nrow(ij)+1), cols = 1:ncol(ij), rule = paste0("$", LETTERS[which(colnames(ij) == "shade")], "1==", 1), style = lsty1)
openxlsx::openXL(wb)


# RNA seq validation----- 
plot_col <- function(sheetname = '5_SEMA3F_quantification', range = NULL){
  basepath <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/arenite_competitive/RNA-seq_analysis'
  tbl <- readxl::read_excel(file.path(basepath, 'joinMajiq_Leafcutter_rosaMacriteria.xlsx'), sheet = sheetname, range = range) %>% 
    group_by(samp) %>% 
    mutate(rep = row_number()) %>% 
    ungroup() %>% 
    pivot_longer(-c(samp,rep), 'band', 'intensity') %>% 
    mutate(samp = factor(samp, levels = c('ctrl', 'trt')),
           band = factor(band, levels = c('inclusion', 'exclusion'))) %>% 
    group_by(samp, rep) %>% 
    mutate(perc = value / sum(value[band == 'inclusion'], value[band == 'exclusion']) * 100) %>% 
    ungroup()
  testLab <- tbl %>% 
    group_by(band) %>% 
    rstatix::t_test(perc ~ samp) %>% 
    rstatix::add_significance("p")
  ggpubr::ggbarplot(tbl,x = "samp", y = "perc", fill = "band", add = "mean_se", position = position_stack(), width = .5)+
    labs(x = '', y = 'Ratio(%)')+
    # annotate(geom="text", x = 0.5, y = 1.3, hjust = "inward", vjust = "inward", label = testLab$ttestlab[1])+
    coord_cartesian(y = c(0, 110))+
    scale_fill_manual(name = '', values = c('inclusion' = '#C3DBFD', 'exclusion' = '#FFFFF9'))+
    scale_y_continuous(breaks = c(0, 50, 100))+
    theme(legend.position = 'none')+
    ggpubr::stat_pvalue_manual(testLab[1,], y.position = 110, step.increase = 0.1, label = "p.signif")
}




s_10uM_293t <- plot_col('5_SEMA3F_quantification')
ggsave(filename = file.path(base,'PhD_study/conferences/23ASMS/pngs/validation_HEK293T_10uM.svg'), device = 'svg', plot = s_10uM_293t)

s_10uM_hepG2 <- plot_col(sheetname = 'HepG2_10uM', range = 'R1C1:R7C3') #SEMA3F
ggsave(filename = file.path(base,'PhD_study/conferences/23ASMS/pngs/validation_HepG2_10uM.svg'), device = 'svg', plot = s_10uM_hepG2)



