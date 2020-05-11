source("functions.R")
t_cfg
diro = file.path(dird, 'design')

#{{{ ChIP-Seq, DAP-Seq and ATAC-Seq
#{{{ ca19a
yid = 'ca19a'
ti = read_samplelist(yid)
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    select(group=Tissue,replicate=Replicate,fastq_1,fastq_2)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ ca19a3
yid = 'ca19a3'
ti = read_samplelist(yid)
ctl = 'input'
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(group=str_c(Tissue,Treatment,sep="_")) %>%
    mutate(antibody=Treatment, control=str_c(Tissue, ctl, sep="_")) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ ca19a4
yid = 'ca19a4'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=ctl) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    mutate(group=Treatment) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cd18a
yid = 'cd18a'
ti = read_samplelist(yid)
ctl = 'halo'
to = ti %>% mutate(Treatment=str_replace(Treatment, "GST-HALO", 'halo')) %>%
    mutate(Treatment=str_replace(Treatment, '^Zm', '')) %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=ctl) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    mutate(group=Treatment) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp12a
yid = 'cp12a'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=str_c(Tissue,ctl,sep="_")) %>%
    mutate(group=str_c(Tissue,Treatment,sep="_")) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp12b
yid = 'cp12b'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(group=ifelse(Genotype=='rr', 'p1', 'control')) %>%
    mutate(antibody=group, Treatment=group, control=ctl) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp14g
yid = 'cp14g'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>% separate(Treatment, c("Treatment", 'suf'), sep='_') %>%
    mutate(Treatment=ifelse(Treatment=='chip', 'ra1', 'control')) %>%
    mutate(Tissue=str_c(Tissue,suf, sep="_")) %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=str_c(Tissue,ctl,sep="_")) %>%
    mutate(group=str_c(Tissue,Treatment,sep="_")) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp15a
yid = 'cp15a'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(Treatment=ifelse(Treatment=='chip', 'fea4', ctl)) %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=ctl) %>%
    mutate(group=Treatment) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp16a - dunno how to do this
yid = 'cp16a'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(Treatment=ifelse(Treatment=='chip', 'fea4', ctl)) %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=ctl) %>%
    mutate(group=Treatment) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp18a
yid = 'cp18a'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(Treatment=ifelse(Genotype=='o2', 'o2', ctl)) %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=ctl) %>%
    mutate(group=Treatment) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp18b
yid = 'cp18b'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(Treatment=ifelse(Treatment=='chip', 'bzip22', ctl)) %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=ctl) %>%
    mutate(group=Treatment) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp19a
yid = 'cp19a'
ti = read_samplelist(yid)
ctl = 'control'
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control='') %>%
    mutate(group=str_c(Tissue, Treatment, sep='_')) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cp19b
yid = 'cp19b'
ti = read_samplelist(yid)
ctl = 'input'
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(antibody=Treatment, control=ctl) %>%
    mutate(group=Treatment) %>%
    mutate(antibody=ifelse(Treatment==ctl, '', antibody)) %>%
    mutate(control=ifelse(Treatment==ctl, '', control)) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2, antibody, control)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}

#{{{ cm16a
yid = 'cm16a'
ti = read_samplelist(yid)
to = ti %>%
    mutate(fastq_1 = map2_chr(SampleID, paired, get_fastq, read1=T, yid=yid)) %>%
    mutate(fastq_2 = map2_chr(SampleID, paired, get_fastq, read1=F, yid=yid)) %>%
    mutate(group=str_c(Tissue,Treatment,sep='_')) %>%
    select(group,replicate=Replicate,fastq_1,fastq_2)
to %>% select(-fastq_1) %>% print(n=50)

fo = sprintf("%s/%s.csv", diro, yid)
write_csv(to, fo)
#}}}
#}}}






