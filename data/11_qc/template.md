## Pipeline outputs:
* [sample meta](00.meta.tsv)
* [alignment BAM stats](stats/bamstats.tsv)
* [number of peaks called by MACS2](stats/macs2_count.tsv)
* [sample FRiP scores](stats/macs2_frip.tsv)
* [IGV session file](igv_session.xml)

File location on MSI:
`/home/springer/zhoux379/projects/chipseq/data/11_qc/**/`
* `bigwig/`: bigWig files
* `peaks/`: peak BED files
* `peaks_annotated/`: (HOMER) annotated peak BED files
* `peaks_annotated/`: consensus peaks (if available)

[Result file walk through](https://github.com/orionzhou/chipseq/blob/master/output.md)
