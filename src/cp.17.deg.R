source("functions.R")
dirw = file.path(dird, '17_chipseq_deg')

fc = file.path(dirw, '../15_collect/cp7.rds')
tc = readRDS(fc)
tc %>% print(n=50)

fd = file.path(dirw, 'degs.rds')
td = readRDS(fd)
td %>% print(n=50)

tw = td %>% mutate(yid = str_replace(yid, '2$', '')) %>%
    rename(rnaseq.tissue=tissue) %>% select(-reg.gid) %>%
    inner_join(tc, by=c('yid','tf')) %>% select(-yid) %>%
    rename(chipseq.tissue=tissue, chipseq.rep = rep) %>%
    print(n=40)

sum_chip_deg <- function(ds, pk) {
    #{{{
    gids = ds %>% filter(padj < .01) %>% pull(gid)
    n_gene = nrow(ds)
    n_gene_DE = length(gids)
    n_peak = nrow(pk)
    pk1 = pk %>% filter(abs(dist_tss) < 5000)
    n_peak_5kb = nrow(pk1)
    pk2 = pk1 %>% filter(gid %in% gids)
    n_peak_5kb_DE = nrow(pk2)
    n_gene_DE_peak = length(unique(pk2$gid))
    tibble(n_gene=n_gene,n_gene_DE=n_gene_DE,n_peak,n_peak_5kb=n_peak_5kb,
        n_peak_5kb_DE=n_peak_5kb_DE,n_gene_DE_peak=n_gene_DE_peak)
    #}}}
}
to = tw %>% mutate(res = map2(ds, peaks, sum_chip_deg)) %>%
    select(-ds,-peaks) %>% unnest(res) %>%
    mutate(pct_deg = n_gene_DE_peak/n_gene_DE,
           pct_peak = n_peak_5kb_DE/n_peak_5kb) %>%
    print(n=20, width=Inf)

fo = file.path(dirw, '10.sum.tsv')
write_tsv(to, fo)
