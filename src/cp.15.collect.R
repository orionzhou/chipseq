source("functions.R")
diri = file.path(dird, '11_qc')
dirw = file.path(dird, '15_collect')
#tf = read_tf_info() %>% mutate(name=ifelse(is.na(name), tf, name))


#{{{ extract FRiP scores
yids = c("cp12a",'cp12b','cp14g','cp15a')#,'cp15b','cp18a','ca19a4')
ti = tibble(yid = yids) %>%
    mutate(meta = map(yid, read_chipseq_meta)) %>%
    mutate(bamstat = map(yid, read_chipseq_bamstats)) %>%
    mutate(macs2 = map(yid, read_chipseq_macs2))

tih = ti %>% select(yid, meta) %>% unnest(meta) %>%
    rename(rep = Replicate) %>%
    distinct(yid,group,rep,antibody,control)
tim = ti %>% select(yid, macs2) %>% unnest(macs2) %>%
    separate(group, c("group", 'rep'), sep='_R') %>%
    mutate(rep = as.integer(rep))

to = t_cfg %>% select(yid,author,year) %>% inner_join(tim,by='yid') %>%
    inner_join(tih, by=c('yid','group','rep'))

fo = file.path(dirw, '05_macs2.tsv')
write_tsv(to, fo, na='')
#}}}

#{{{ merge NF chipseq/dapseq results
read_chipseq_peaks <- function(yid, min_frip=.05, tf_info=tf, diri = '~/projects/nf/data/out') {
    #{{{
    fi = sprintf("%s/%s.rds", diri, yid)
    res = readRDS(fi)
    grps = res$frip %>% filter(frip >= min_frip) %>% select(-frip)
    to = res$peaks %>% inner_join(grps, by=c('group','rep')) %>%
        select(group, rep, chrom, start, end, srd, score,tgt.gid=gid,dist_tss)
    tf1 = tf %>% filter(yid==!!yid) %>% select(tf, reg.gid=gid, ctag=reference)
    if(yid == 'ca19a4')
        to1 = to %>% rename(tf = group, note = rep)
    else if (yid == 'cp12b')
        to1 = to %>% mutate(group = str_to_upper(group)) %>%
            rename(tf=group, note=rep)
    else if (yid == 'cp14g')
        to1 = to %>% mutate(tf = 'RA1', note=str_replace(group, '_control','')) %>%
            select(-group,-rep)
    else if (yid == 'cp15a')
        to1 = to %>% mutate(tf='FEA4') %>% select(-group) %>% rename(note=rep)
    else if (yid == 'cp15b' | yid == 'cp18a')
        to1 = to %>% mutate(tf='O2') %>% select(-group) %>% rename(note=rep)
    to1 %>% inner_join(tf1, by='tf') %>%
            select(ctag, tf, note, reg.gid, everything())
    #}}}
}

yids = c("cp12a",'cp12b','cp14g','cp15a','cp15b','cp18a','ca19a4')
tfmap = c("cp12a"="KN1","cp12b"='P1','cp14g'='RA1','cp15a'='FEA4','cp15b'='O2','cp18a'='O2')
hdrs = c('chrom','start','end','peakid','score','srd','grp',"gid","dist_tss")
ti = tibble(yid = yids) %>%
    mutate(fi = sprintf("%s/%s.bed", diri, yid)) %>%
    mutate(r = map(fi, read_tsv, col_names=hdrs))

to = ti %>%
    select(-fi) %>% unnest(r) %>% select(-peakid) %>%
    separate(grp, c('grp','rep'), sep=-2, extra='merge') %>%
    mutate(rep = as.integer(str_replace(rep, '_', ''))) %>%
    group_by(yid, grp, rep) %>%
    nest() %>% rename(peaks = data) %>% ungroup() %>%
    mutate(tf = ifelse(yid %in% names(tfmap), tfmap[yid], grp)) %>%
    mutate(tissue = ifelse(yid %in% c("cp12a","cp14g"), grp, '')) %>%
    mutate(tissue = str_replace(tissue, '_control', '')) %>%
    select(yid, tf, tissue, rep, peaks)
to %>% print(n=44)

fo = file.path(dirw, "cp7.rds")
saveRDS(to, fo)
#}}}


#{{{ Ricci2019 DAP-Seq 27 TFs
fi = file.path(dirw, 'Ricci2019.dapseq.txt')
ti = read_tsv(fi, col_names='fname') %>%
    mutate(fi = sprintf("%s/Ricci2019/%s", dirw, fname)) %>%
    mutate(tf = str_replace(fname,'\\.dap\\.bed.gz','')) %>%
    mutate(tf = str_replace(tf,'^GSM.*DAP_','')) %>%
    inner_join(tf, by='tf') %>% select(-tf) %>% rename(tf=gid)

to = ti %>%
    mutate(data=map(fi, read_tsv, col_names=c('chrom','start','end'))) %>%
    select(tf,data) %>% unnest(data) %>%
    select(chrom,start,end,tf)
fo = file.path(dirw, 'Ricci2019_tmp/01.raw.bed')
write_tsv(to, fo, col_names=F)
#}}}


