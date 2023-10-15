# agnes_scRNA
# agnes_scRNA
# Preprocessing
The data was obtained in the FASTQ format from 10X. I ran FASTQC to assess the quality first and then passed it through cellranger count. 

```{bash} 
cellranger count --id=DMSO --fastqs=BL_Chrm_007_Nephew_1 --sample=DMSO --transcriptome=refdata-gex-GRCh38-2020-A
cellranger count --id=EZH2i --fastqs=BL_Chrm_007_Nephew_2 --sample=GSK126 --transcriptome=refdata-gex-GRCh38-2020-A
cellranger count --id=RACi --fastqs=BL_Chrm_007_Nephew_3 --sample=RACi --transcriptome=refdata-gex-GRCh38-2020-A
cellranger count --id=EZH2-RACi --fastqs=BL_Chrm_007_Nephew_4 --sample=combo --transcriptome=refdata-gex-GRCh38-2020-A
```
