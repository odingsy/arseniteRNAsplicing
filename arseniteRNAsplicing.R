# samples 
# As_competitive_nuc_F.raw
# As_competitive_nuc_R.raw
# nucNAP5F_20230207025301.raw
# nucNAP5R_20230207070239.raw


library(tidyverse)
hs <- UniProt.ws::UniProt.ws(9606)
sepCollapse <- '||' # separator for concatenation annotation information. 

# load the data ----
base <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/conferences/23ASMS'
# raw <- file.path(base, 'PhD_study/conferences/23ASMS/rds') 
split_pp <- readr::read_delim(file.path(base, 'data/230329_AsCompetitive_combinedSplitSearch', 'peptides.txt'), na = c("", "NA","#NUM!"))
split_pg <- readr::read_delim(file.path(base, 'data/230329_AsCompetitive_combinedSplitSearch', 'proteinGroups.txt'), na = c("", "NA","#NUM!")) %>%
  filter(!(grepl("(CON|REV)__", `Protein IDs`) | grepl("(CON|REV)__", `Majority protein IDs`)))
# Dong et al.(2022) CRT. SI table 2 (https://pubs.acs.org/doi/suppl/10.1021/acs.chemrestox.2c00244/suppl_file/tx2c00244_si_002.xlsx)
dxjtbl <- readxl::read_xlsx(file.path(base, 'data/tx2c00244_si_002.xlsx'), sheet = 'Table S1', col_types = c(rep("text",2), rep("numeric",8))) %>%
  filter(!grepl(pattern = '^REV__', x = `Uniprot ID`))
# Burger et al.(2025) Cell. Table S1 (https://www.cell.com/cms/10.1016/j.cell.2024.11.025/attachment/a959124f-f54e-42cc-8a81-d771183f3c10/mmc1.xlsx)
bl_ht <- readxl::read_xlsx(file.path(base, 'data/mmc1.xlsx'), sheet = 'HCT116-TPEN') %>% filter(`p.adj.signif` != 'ns')
bl_hept <- readxl::read_xlsx(file.path(base, 'data/mmc1.xlsx'), sheet = 'Hepatocytes-TPEN') %>% filter(`p.adj.signif` != 'ns')
bl_prot <- readxl::read_xlsx(file.path(base, 'data/mmc1.xlsx'), sheet = 'Prostate-TPEN') %>% filter(`p.signif` != 'ns')
bl_hz <- readxl::read_xlsx(file.path(base, 'data/mmc1.xlsx'), sheet = 'HCT116-ZnCl2') %>% filter(`p.adj.signif` != 'ns')
bl_hepz <- readxl::read_xlsx(file.path(base, 'data/mmc1.xlsx'), sheet = 'Hepatocytes-ZnCl2') %>% filter(`p.adj.signif` != 'ns')
bl_proz <- readxl::read_xlsx(file.path(base, 'data/mmc1.xlsx'), sheet = 'Prostate-ZnCl2') %>% filter(`p.signif` != 'ns')

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
geneSymbol <- readRDS(file.path(base, 'rds/geneSymbol.rds'))
zf <- readRDS(file.path(base, 'rds/zf.rds'))
tpen <- unique(c(bl_ht$uniprot_accession, bl_hept$uniprot_accession, bl_prot$uniprot_accession)) 
zn <- unique(c(bl_hz$uniprot_accession, bl_hepz$uniprot_accession, bl_proz$uniprot_accession)) 
dxjList <- dxjtbl %>%
  separate_rows(`Uniprot ID`, sep = ';') %>%
  pull(`Uniprot ID`)

nuctbl <- nuctbl %>%
  rowwise() %>%
  mutate(gn = geneSymbol$Gene.Names..primary.[geneSymbol$From %in% unlist(strsplit(`Protein IDs`, split = ';'))] %>% unique(.) %>% .[complete.cases(.)] %>% paste0(., collapse = sepCollapse),
         `in Dong et al.(2022)?` = ifelse(any(unlist(strsplit(`Protein IDs`, split = ';')) %in% dxjList), 'yes', 'no'),
         `constitutive Zinc binding?` = ifelse(str_replace_all(`Protein IDs`, pattern = ';', replacement = '|') %>% 
                                                 str_detect(tpen, .) %>% 
                                                 any(), 'yes', 'no'),
         `inducible Zinc binding?` = ifelse(str_replace_all(`Protein IDs`, pattern = ';', replacement = '|') %>% 
                                              str_detect(zn, .)%>% 
                                              any(), 'yes', 'no'),
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
  filter(`num<1.5` <= 2 & mean > 1.5) %>% 
  filter(`constitutive Zinc binding?` == 'yes' | `inducible Zinc binding?` == 'yes' | nchar(zincFinger) !=0) %>% 
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

# overlapping with Dong et al 2022
lst <- tblout %>% 
  separate_rows(`Protein IDs`, sep = ';') %>%
  pull(`Protein IDs`)

sum(lst %in% dxjList)
sum(!lst %in% dxjList)
sum(dxjList %in% lst)
sum(!dxjList %in% lst)

# Figure S2 GO analysis of the shortlist protein ----
gene.df <- clusterProfiler::bitr(tblout %>% separate_rows(`Protein IDs`, sep = ';') %>% pull(`Protein IDs`), fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Hs.eg.db::org.Hs.eg.db)
bp5 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "BP", level = 5, readable = TRUE); bp5@result <- bp5@result %>% arrange(desc(Count))
bp <- enrichplot:::barplot.enrichResult(bp5, font.size = 10, showCategory=12, label_format = 50)

cc7 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "CC", level = 7, readable = TRUE); cc7@result <- cc7@result %>% arrange(desc(Count))
cc <- enrichplot:::barplot.enrichResult(cc7, font.size = 10, showCategory=12, label_format = 50)

mf7 <- clusterProfiler::groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "MF", level = 7, readable = TRUE); mf7@result <- mf7@result %>% arrange(desc(Count))
mf <- enrichplot:::barplot.enrichResult(mf7, font.size = 10, showCategory=12, label_format = 50)

ggsave('./pngs/bp.pdf', plot = bp, width = 5, height = 5, units = 'in')
ggsave('./pngs/cc.pdf', plot = cc, width = 5, height = 5, units = 'in')
ggsave('./pngs/mf.pdf', plot = mf, width = 5, height = 5, units = 'in')

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
  mutate(gn = case_when(
    gn == 'HARS1||HARS2' ~ 'HARS1/2',
    gn == 'PTBP1||PTBP3' ~ 'PTBP1/3',
    TRUE ~ gn
  )) %>% 
  ggplot(aes(x = nucF, y = nucR), shape = 20)+
  geom_point(color = 'grey', size = 0.3, alpha = 0.2)+
  geom_point(data = . %>% filter(`Protein IDs`%in% tblout$`Protein IDs`), color = '#F8766D', size = 0.4)+
  ggrepel::geom_text_repel(data = . %>% filter((gn %in% c('SF1','ZRANB2', 'DDX39B', 'PTBP1/3', 'RAD50', 'XRCC1', 'TP53BP1', 'MRE11', 'POLDIP3'))& gn != c('PGK1||PGK2')), aes(label = gn), size= 2, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = Inf)+
  # ggrepel::geom_text_repel(data = . %>% filter(((`Protein IDs`%in% tblout$`Protein IDs` & nucF > 2.2 & nucR > 2.2) | gn %in% c('ZRANB2'))& gn != c('PGK1||PGK2')), aes(label = gn), box.padding = 0.5, min.segment.length = 0.1, max.overlaps = Inf)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0,7), ylim = c(0,7), xintercept = 1, yintercept = 1, ratio = 1)+
  labs(x = "Ratio(Control/Competition), Forward", y = "Ratio(Control/Competition), Reverse")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none")

ggsave(filename = file.path(base, 'pngs/protein_dotplot.svg'), device = 'svg', plot = p, width = 5, height = 5, units = 'in')
# plotly::ggplotly(p)


# Figure S1, find the peptides -----
pptbl <- tribble(
  ~protein, ~peptide, ~group,
  # RNA splicing
  'SF1', 'QGIETPEDQNDLR', 'splicing',
  'SF1', 'SITNTTVCTK', 'splicing',
  'PTBP1', 'LSLDGQNIYNACCTLR', 'splicing',
  'PTBP1/3', 'VTNLLMLK','splicing',
  'DDX39B', 'QVMMFSATLSK', 'splicing',
  'ZRANB2', 'LDEDEDEDDADLSK', 'splicing',
  'ZRANB2', 'AVGPASILK', 'splicing',
  # DNA repair
  'RAD50', 'VLASLIIR', 'repair',
  'RAD50', 'MSILGVR', 'repair',
  'TP53BP1', 'ADDPLRLDQELQQPQTQEK', 'repair',
  'TP53BP1', 'QDATVQTER', 'repair',
  'XRCC1', 'AQGAVTGKPR', 'repair',
  'MRE11', 'GNDTFVTLDEILR', 'repair',
  'MRE11', 'IDISPVLLQK', 'repair',
  # 'XRN2', 'EGMEAAVEK',
  # 'XRN2', 'NLTVILSDASAPGEGEHK',
  # 'ZC3H15', 'DEELEKDTMDNWDEK',
  # 'ZC3H15', 'EVFEFRPELVNDDDEEADDTR',
  # 'METAP1', 'LFHTAPNVPHYAK',
  # 'RBM10', 'LDQQTLPLGGR',
  # 'RBM10', 'MLPQAATEDDIR',
  'POLD1', 'ITVALPR', 'replication',
  'POLD1', 'EAADWVSGHFPSPIR', 'replication',
  'POLE3', 'AERPEDLNLPNAVITR', 'replication',
  'POLE3', 'EALPDGVNISK', 'replication'
)

splicing <- hcl.colors(hcl.pals("sequential")[10], n = 8)[c(4,4,6,6,5,3,3)]#|> scales::show_col()
repair <- hcl.colors(hcl.pals("sequential")[5], n = 8)[c(3,3,5,5,4,2,2)]#|> scales::show_col()
replication <- hcl.colors(hcl.pals("sequential")[7], n = 8)[c(4,4,5,5)]#|> scales::show_col()

p <- split_pp %>% 
  dplyr::filter(Sequence %in% pptbl$peptide) %>% 
  mutate(Sequence = factor(Sequence, levels = pptbl$peptide)) %>% 
  arrange(Sequence) %>%
  mutate(gn = pptbl$protein) %>% 
  dplyr::select(Sequence, gn, dplyr::matches('Ratio H/L normalized nuc._1'), dplyr::matches('Ratio H/L normalized nuc._2')) %>% 
  mutate(across(matches('nucF'), ~{1/.x})) %>%
  mutate(peptide_gn = sprintf(paste0('%s', '%',7-nchar(gn)+10, 's'), Sequence, gn) %>% factor(., levels = rev(.))) %>% 
  rename_with(~{str_replace(.x, "Ratio H/L normalized ", "")}, matches('nuc(F|R)_\\d')) %>% 
  pivot_longer(cols = matches('nuc(F|R)_(1|2)'), names_to = "samp", values_to = "ratio_normalized") %>%
  mutate(samp = factor(samp, levels = c('nucF_1', 'nucR_1', 'nucF_2', 'nucR_2'))) %>% 
  group_by(gn) %>%
  ggplot(aes(peptide_gn, ratio_normalized))+
  stat_summary(aes(fill = peptide_gn), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  geom_hline(yintercept=1.5, linetype="dashed",color = "red", linewidth = .5)+
  geom_point(aes(color = samp), size = 2)+
  coord_flip()+
  ylab('control/competitor')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,2,7,8)])+ # dot color
  scale_fill_manual(values = c(replication, repair, splicing))+
  guides(fill = FALSE, color = guide_legend(override.aes = list(size = 0.5)))+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))

ggsave(filename = file.path(base, 'pngs/peptide_barplot.svg'), device = 'svg', plot = p, width = 5, height = 5, units = 'in')

