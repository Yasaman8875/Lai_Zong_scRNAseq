#!/bin/bash
set -e
set -u
set -o pipefail

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif velocyto run10x -@ 10 --samtools-memory 5000 -m hg38_rmsk.gtf DMSO refdata-gex-GRCh38-2020-A/genes/genes.gtf
singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif velocyto run10x -@ 10 --samtools-memory 5000 -m hg38_rmsk.gtf EZH2i refdata-gex-GRCh38-2020-A/genes/genes.gtf
singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif velocyto run10x -@ 10 --samtools-memory 5000 -m hg38_rmsk.gtf RACi refdata-gex-GRCh38-2020-A/genes/genes.gtf
singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif velocyto run10x -@ 10 --samtools-memory 5000 -m hg38_rmsk.gtf EZH2-RACi refdata-gex-GRCh38-2020-A/genes/genes.gtf
