require(devtools)
load_all("~/git/rmaize")
require(ggforce)
lib = 'chipseq'
dirp = sprintf("~/projects/%s", lib)
dird = file.path(dirp, 'data')
t_cfg = read_gspread_master(lib = c("chipseq",'dapseq'))

read_chipseq_meta <- function(yid, diri=file.path(dird,'11_qc'))
    read_tsv(sprintf("%s/%s/00.meta.tsv", diri, yid))
read_chipseq_bamstats <- function(yid, diri=file.path(dird,'11_qc'))
    read_tsv(sprintf("%s/%s/stats/bamstats.tsv", diri, yid))
read_chipseq_macs2 <- function(yid, diri=file.path(dird,'11_qc')) {
    #{{{
    ti1 = read_tsv(sprintf("%s/%s/stats/macs2_count.tsv", diri, yid), col_names=c('group','count')) %>%
        separate(group, c('fn','group'), sep=":", fill='left') %>% select(-fn)
    ti2 = read_tsv(sprintf("%s/%s/stats/macs2_frip.tsv", diri, yid), col_names=c('group','frip')) %>%
        separate(group, c('fn','group'), sep=":", fill='left') %>% select(-fn)
    ti1 %>% inner_join(ti2, by='group')
    #}}}
}
