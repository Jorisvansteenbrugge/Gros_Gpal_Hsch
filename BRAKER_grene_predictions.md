The prediction of genes in the assemblies was done with BRAKER2 (https://github.com/Gaius-Augustus/BRAKER). Gene prediction was supported by publicly available RNAseq data under accesions 

RNAseq data was mapped on the reference genomes using hisat2 :

```
$ hisat2-build -p 8 G_pallida_D383.fasta G_pallida_D383

$ mkdir rnaseq
$ for fwd_read in *_1.fastq.gz
$   do
$     read_base=$(basename $fwd_read _1.fastq.gz)
$     hisat2 -x G_pallida_D383 -1 $fwd_read -2 "$read_base"_2.fastq.gz --dta -p 8 \
        | samtools sort -@ 2 -m 30G -o rnaseq/"$read_base"_gpD383.bam
$     samtools index rnaseq/"$read_base"_gpD383.bam
$ samtools merge rnaseq_gpD383.bam *.bam
$   done
```


```
$ braker.pl --genome=G_pallida_D383.fasta --species=G_pallida_D383 --cores=10 --gff3  --bam=rnaseq_gpD383.bam
