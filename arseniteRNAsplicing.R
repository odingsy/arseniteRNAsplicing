# maxquant setting
# 230329 combine search of two nuclear proteome pulldown experiments
# - Multiplicity:2
# - max. labeled AAs:3 (default) 
# - 'Heavy label': Arg6, Lys8
# - Max number of modifications per peptide: 5 (default)
# - max missed cleavages: 2 (default)
# - min peptide length: 7 (default)
# - UP000005640_9606 (10/25/22) 
# - match between runs: default setting
# - 5 cores
# - split search 

# samples 
# As_competitive_nuc_F.raw
# As_competitive_nuc_R.raw
# nucNAP5F_20230207025301.raw
# nucNAP5R_20230207070239.raw


library(tidyverse)
hs <- UniProt.ws::UniProt.ws(9606)
sepCollapse <- '||' # separator for concatenation annotation information. 

# load the data ----
base <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/conferences/23ASMS/data'
# raw <- file.path(base, 'PhD_study/conferences/23ASMS/rds') 
split_pp <- readr::read_delim(file.path(base, '230329_AsCompetitive_combinedSplitSearch', 'peptides.txt'), na = c("", "NA","#NUM!"))
split_pg <- readr::read_delim(file.path(base, '230329_AsCompetitive_combinedSplitSearch', 'proteinGroups.txt'), na = c("", "NA","#NUM!")) %>%
  filter(!(grepl("(CON|REV)__", `Protein IDs`) | grepl("(CON|REV)__", `Majority protein IDs`)))
# supplientary table from Dong et al. (https://pubs.acs.org/doi/suppl/10.1021/acs.chemrestox.2c00244/suppl_file/tx2c00244_si_002.xlsx)
dxjList <- readxl::read_xlsx(file.path(base, 'tx2c00244_si_002.xlsx'), sheet = 'Table S1') %>%
  filter(!grepl(pattern = '^REV__', x = `Uniprot ID`)) %>%
  separate_rows(`Uniprot ID`, sep = ';') %>%
  pull(`Uniprot ID`)


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
         sd = map_dbl(data, ~{sd(.$ratio_normalized, na.rm = TRUE) %>% as.numeric(.)}),
         `num<1.5` = map_dbl(data, ~{ sum(c(is.na(.$ratio_normalized), sum(.$ratio_normalized < 1.5, na.rm = TRUE))) %>% as.numeric(.)}),
         numOfNA = map_dbl(data, ~{is.na(.$ratio_normalized) %>% sum(.) %>% as.numeric(.)})) %>%
  # filter(`num<1.5` <= 2 & (mean > 2.5 | (mean > 1.5 & rsd < 0.2))) %>%
  unnest(data) %>%
  pivot_wider(names_from = "samp", values_from = "ratio_normalized") %>%
  relocate(`Protein IDs`, `Protein names`, matches('nuc(F|R)_1'), matches('nuc(F|R)_2'))



# obtain annotation information from UniProt
upid <- nuctbl %>% separate_rows(`Protein IDs`) %>% pull(`Protein IDs`)
# geneSymbol <- UniProt.ws::select(x = hs, keys =  upid, columns = c('gene_primary'), keytype = 'UniProtKB')
# zf <- UniProt.ws::select(x = hs, keys = upid, columns = c('ft_zn_fing'), keytype = 'UniProtKB')
# saveRDS(geneSymbol, file.path(base, 'PhD_study/conferences/23ASMS/rds/geneSymbol.rds'))
# saveRDS(zf, file.path(base, 'PhD_study/conferences/23ASMS/rds/zf.rds'))
geneSymbol <- readRDS(file.path(base, 'PhD_study/conferences/23ASMS/rds/geneSymbol.rds'))
zf <- readRDS(file.path(base, 'PhD_study/conferences/23ASMS/rds/zf.rds'))

nuctbl <- nuctbl %>%
  rowwise() %>%
  mutate(gn = geneSymbol$Gene.Names..primary.[geneSymbol$From %in% unlist(strsplit(`Protein IDs`, split = ';'))] %>% unique(.) %>% .[complete.cases(.)] %>% paste0(., collapse = sepCollapse),
         `in Dong et al.(2022)?` = ifelse(any(unlist(strsplit(`Protein IDs`, split = ';')) %in% dxjList), 'yes', 'no'),
         zincFinger = zf$Zinc.finger[zf$From %in% unlist(strsplit(`Protein IDs`, split = ';'))] %>% unique(.) %>% .[complete.cases(.)] %>% paste0(., collapse = sepCollapse)) %>%
  ungroup() %>% 
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
  )) 

# Table S1: a table summarizing all the interactive candidates ----
tblout <- nuctbl %>% 
  filter(nchar(zincFinger) !=0 & `num<1.5` <= 2 & mean > 1.5) %>% 
  mutate(across(where(is.double), ~signif(.x, digits = 3))) %>% 
  arrange(desc(mean)) %>%
  mutate(gn = factor(gn, levels = gn)) %>%
  mutate(rnum = row_number(),
         shade = case_when(rnum %% 2 == 1 ~ 0, rnum %% 2 == 0 ~ 1)) %>%
  select(-rnum) %>% 
  rename(Gene = gn) %>% 
  select(Gene,`Protein IDs`, starts_with('nuc'), mean, sd, zincFinger, everything())
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, 'test')
lsty1 <- openxlsx::createStyle(bgFill = "#d6d2d2")
openxlsx::writeData(wb, 'test', tblout)
openxlsx::conditionalFormatting(wb, 'test', rows = 1:(nrow(tblout)+1), cols = 1:ncol(tblout), rule = paste0("$", LETTERS[which(colnames(tblout) == "shade")], "1==", 1), style = lsty1)
openxlsx::openXL(wb)

# Figure S2 GO analysis of the shortlist protein ----
gene.df <- clusterProfiler::bitr(tblout$`Protein IDs`, fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Hs.eg.db::org.Hs.eg.db)
bp7 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "BP", level = 7, readable = TRUE); bp7@result <- bp7@result %>% arrange(desc(Count))
cc6 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "CC", level = 6, readable = TRUE); cc6@result <- cc6@result %>% arrange(desc(Count))
mf6 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "MF", level = 6, readable = TRUE); mf6@result <- mf6@result %>% arrange(desc(Count))

bp <- enrichplot:::barplot.enrichResult(bp7, font.size = 10, showCategory=10, label_format = 50)
cc <- enrichplot:::barplot.enrichResult(cc6, font.size = 10, showCategory=10, label_format = 50)
mf <- enrichplot:::barplot.enrichResult(mf6, font.size = 10, showCategory=10, label_format = 50)

ggsave('./pngs/bp.pdf', bp)
ggsave('./pngs/cc.pdf', cc)
ggsave('./pngs/mf.pdf', mf)

wb <- openxlsx::createWorkbook();
openxlsx::addWorksheet(wb, 'BiologicalProcess'); openxlsx::writeData(wb, 'BiologicalProcess', head(bp7, 20));
openxlsx::addWorksheet(wb, 'CellularComponent'); openxlsx::writeData(wb, 'CellularComponent', head(cc6, 20));
openxlsx::addWorksheet(wb, 'MolecularFunction'); openxlsx::writeData(wb, 'MolecularFunction', head(mf6, 20));
openxlsx::openXL(wb)

# Figure 2, combine nucF_1, nucF_2 and nucR_1, nucR_2 for dotplot -----
p <- nuctbl %>%
  rowwise() %>%
  mutate(nucF = mean(c(nucF_1, nucF_2), na.rm = TRUE), 
         nucR = mean(c(nucR_1, nucR_2), na.rm = TRUE)) %>% 
  ggplot(aes(x = nucF, y = nucR))+
  geom_point(color = 'grey', size = 1, alpha = 0.5)+
  geom_point(data = . %>% filter(`Protein IDs`%in% tblout$`Protein IDs`), aes(color = `Arsenite replaceable zinc-finger domain`), size = 1)+
  ggrepel::geom_label_repel(data = . %>% filter((`Protein IDs`%in% tblout$`Protein IDs` & nucF > 1.8 & nucR > 1.8)| gn %in% c('ZRANB2')), aes(label = gn), box.padding = 0.5, max.overlaps = Inf)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0,6), ylim = c(0,6), xintercept = 1, yintercept = 1, ratio = 1)+
  labs(x = "1/forward = ctrl(L)/competitive(H)", y = "reverse = ctrl(H)/competitive(L)")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none")

ggsave(filename = file.path(base,'PhD_study/conferences/23ASMS/pngs/protein_dotplot.svg'), device = 'svg', plot = p, width = 5, height = 5, units = 'in')
# plotly::ggplotly(p)


# Figure S1, find the peptides -----
protn <- c('SF1', 'SF1',
           'PTBP1/3', 'PTBP1/3',
           'XRN2', 'XRN2',
           'ZC3H15', 'ZC3H15',
           'METAP1',
           'MRE11', 'MRE11',
           'RBM10', 'RBM10')
ppn <- c('QGIETPEDQNDLR', 'SITNTTVCTK', 
         'LSLDGQNIYNACCTLR', 'VTNLLMLK',
         'EGMEAAVEK', 'NLTVILSDASAPGEGEHK',
         'DEELEKDTMDNWDEK', 'EVFEFRPELVNDDDEEADDTR',
         'LFHTAPNVPHYAK',
         'GNDTFVTLDEILR', 'IDISPVLLQK',
         'LDQQTLPLGGR', 'MLPQAATEDDIR')

split_pp %>% 
  dplyr::filter(Sequence %in% ppn) %>% 
  mutate(Sequence = factor(Sequence, levels = ppn)) %>% 
  arrange(Sequence) %>%
  mutate(gn = protn) %>% 
  dplyr::select(Sequence, gn, dplyr::matches('Ratio H/L normalized nuc._1'), dplyr::matches('Ratio H/L normalized nuc._2')) %>% 
  mutate(across(matches('nucF'), ~{1/.x})) %>%
  mutate(peptide_gn = sprintf(paste0('%s', '%',7-nchar(gn)+10, 's'), Sequence, gn) %>% factor(., levels = rev(.))) %>% 
  rename_with(~{str_replace(.x, "Ratio H/L normalized ", "")}, matches('nuc(F|R)_\\d')) %>% 
  pivot_longer(cols = matches('nuc(F|R)_(1|2)'), names_to = "samp", values_to = "ratio_normalized") %>%
  mutate(samp = factor(samp, levels = c('nucF_1', 'nucR_1', 'nucF_2', 'nucR_2'))) %>% 
  group_by(gn) %>%
  ggplot(aes(peptide_gn, ratio_normalized))+
  stat_summary(aes(fill = reorder(gn, Sequence)), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  geom_hline(yintercept=1.5, linetype="dashed",color = "blue", linewidth = .5)+
  geom_point(aes(color = samp), size = 2)+
  coord_flip()+
  ylab('control/competitor')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,2,7,8)])+
  guides(color = FALSE, fill = FALSE)+
  theme_classic()
