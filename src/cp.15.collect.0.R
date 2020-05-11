source("functions.R")
diro = file.path(dird, 'out')
dirr = file.path(dird, 'raw')

read_nf_chipseq <- function(yid) {
    #{{{
    ti = read_design(yid)
    th = ti %>% select(group, replicate, antibody) %>% filter(!is.na(antibody)) %>%
        distinct(group, replicate) %>%
        mutate(sid = sprintf("%s_R%d", group, replicate))
    #
    to1 = th %>%
        mutate(fp = sprintf("%s/%s/bwa/mergedLibrary/macs/broadPeak/%s_peaks.annotatePeaks.txt", dirr, yid, sid)) %>%
        mutate(data = map(fp, read_tsv, col_names=F, skip=1)) %>%
        select(group, replicate, data) %>% unnest() %>%
        select(group=1,rep=2,peakid=3,chrom=4,start=5,end=6,srd=7,score=8,
            annotation=10,dist_tss=12, gid=13)
    to1 %>% count(group, rep) %>% print(n=50)
    #
    to2 = th %>%
        mutate(fp = sprintf("%s/%s/bwa/mergedLibrary/macs/broadPeak/qc/%s_peaks.FRiP_mqc.tsv", dirr, yid, sid)) %>%
        mutate(data = map(fp, read_tsv, col_names=F, comment='#')) %>%
        select(group, rep=replicate, data) %>% unnest() %>%
        select(group, rep, frip=4)
    to2 %>% print(n=50)
    to2f = to2 %>% filter(frip >= .0)
    to1 %>% count(group, rep) %>% left_join(to2f, by=c('group','rep')) %>% print(n=50)
    #
    frip = to2f
    peaks = to1 %>% inner_join(to2f, by=c('group','rep')) %>%
        replace_na(list(peakid='.')) %>% select(-frip)
    list(design = ti, frip = frip, peaks = peaks)
    #}}}
}

# ChIP-Seq / DAP-Seq
#{{{
yid = 'cp14g'
yid = 'cp15a'
yid = 'cp15b'
yid = 'cp18a'
yid = 'ca19a4'
yid = 'cp12a'
yid = 'cp12b'

res = read_nf_chipseq(yid)
fo = sprintf("%s/%s.rds", diro, yid)
saveRDS(res, file=fo)
#}}}


#{{{ share with Pete
dir1 = file.path(dird, 'out')
fi = sprintf("%s/%s.rds", dir1, yid)
res = readRDS(fi)

dirs = file.path(dirp, 'data2/share')
to = res$peaks %>% mutate(grp=str_c(group,rep,sep='_'), start=start-1) %>%
    select(chrom,start,end,peakid,score,srd,grp,gid,dist_tss) %>%
    arrange(chrom,start,end)
to %>% count(grp)
fo = sprintf("%s/%s.bed", dirs, yid)
write_tsv(to, fo, col_names=F)
#}}}



