### Programs not in this repository

- canu (https://github.com/marbl/canu)
- wtdbg2 (https://github.com/ruanjue/wtdbg2)
- BUSCO (https://busco.ezlab.org/)
- minimap2 (https://github.com/lh3/minimap2)
- purge_haplotigs (https://bitbucket.org/mroachawri/purge_haplotigs/)
- finisherSC (https://github.com/kakitone/finishingTool)
- Hydraslayer/GapFiller (https://github.com/Jorisvansteenbrugge/GapFiller)
- SSPACE-Longread (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-211)
- Arrow/GenomicConsensus (https://github.com/PacificBiosciences/GenomicConsensus)
- Pilon (https://github.com/broadinstitute/pilon)
- BioPython (https://biopython.org/)

# Trimming of reads

Since we sequenced a pool of genomes we perform a read correction step to reduce the number of haplotypes in the reads. Canu is used for the error correction of the reads. Specifically, parameters these parameters were tweaked: `corOutCoverage=200 correctedErrorRate=0.15`.

```
$ canu -correct -p gpD383 -d corrected_reads genomeSize=100m corOutCoverage=200 correctedErrorRate=0.15 \
    -pacbio-raw raw_pb_reads.fasta

```
The corrected reads produced by canu are stored as **correctedReads.fasta.gz**



# Assembly with wtdgb2

For each genome approximately 100 initial assemblies were produced with wtdbg2, optimizing the parameters minimal read length, k-mer size, and minimal read depth. This optimization step is performed with a custom script [optimize_wtdbg2.py](https://github.com/Jorisvansteenbrugge/GROS_genomes/blob/main/Assembly/optimize_wtdbg2.py)

```
$ python Assembly/optimize_wtdbg2.py --reads correctedReads.fasta.gz --wd initial_assemblies -t 8
```

 **Assembly_L5000_p18e6** was selected to continue,

# Purging of Haplotigs

```
$ minimap2 -t 4 -ax map-pb Assembly_L5000_p18e6.fasta raw_pb_reads.fasta --secondary=no | samtools sort -m 30G \
    -@ 2 -o Assembly_L5000_p18e6_minimapped_AllPB.bam
$ purge_haplotigs hist -b Assembly_L5000_p18e6_minimapped_AllPB.bam -g Assembly_L5000_p18e6.fasta -t 8
$ purge_haplotigs cov -i Assembly_L5000_p18e6_minimapped_AllPB.bam.genecov -l 4 -m 15 -h 110
$ purge_haplotigs purge -g Assembly_L5000_p18e6.fasta -c coverage_stats.csv
$ mv curated.fasta G_pallidaD383_PH.fasta
```


# Finishing Tool

```
$ mkdir finisheroutput
$ ln -s correctedReads.fasta.gz finisheroutput/raw_reads.fasta.gz
$ finisherSC.py -par 8 finisheroutput/ ~/tools/miniconda/bin/mummer
```



# Scaffolding SSPACE-Longread

**Gr-Line19**
```
$ mkdir scaffolding
$ SSPACE-LongReadunix.pl -c finisheroutput/improved3.fasta -p correctedReads.fasta.gz \
     -t 8 -b scaffolding/ -o 1000 -g 500 -k1
```

# Gap Filling

**Gr-Line19**
```
$ mkdir gapfilling
$ minimap2 -t <threads> -Hax map-pb scaffolding/scaffolds.fasta correctedReads.fasta.gz \
    | samtools view -hF 256 - | samtools sort -@ 2 -m 30G -o gapfilling/correctedreads_onScaffolds.bam - 
$ samtools index gapfilling/correctedreads_onScaffolds.bam
$ Hydraslayer -o gapfilling/scaffolds_gapfilled.fasta -c 10 -a 60 scaffolding/scaffolds.fasta \
    gapfilling/correctedreads_onScaffolds.bam
```

# Polishing with Arrow

the following is repeated three times. The first run uses the output of the gapfilling step `gapfilling/scaffolds_gapfilled.fasta`. The second and third run use the output fasta file of the first and second run, respectively.
```
$ mkdir Arrow
$ blasr raw_pb_reads.fasta gapfilling/scaffolds_gapfilled.fasta --out Arrow/R1.bam --bam --bestn 10 \
        --minMatch 12 --maxMatch 30 --nproc 8 --minSubreadLength 50 --minAlnLength 50 --minPctSimilarity 70 \
        --minPctAccuracy 70 --hitPolicy randombest --randomSeed 1
$ samtools sort Arrow/R1.bam -@ 2 -m 30G -o Arrow/R1_sort.bam
$ pbindex Arrow/R1_sort.bam
$ arrow Arrow/R1_sort.bam -j 8 --referenceFilename gapfilling/scaffolds_gapfilled.fasta -o Arrow/scaffolds_gapfilled_arrow1.fasta -o Arrow/scaffolds_gapfilled_arrow1.gff -o Arrow/scaffolds_gapfilled_arrow1.fastq
```


## Polishing with Pilon

```
$ mkdir gp_pilon
$ pilonAuto -g Arrow/scaffolds_gapfilled_arrow3.fasta -m 100 -t 8 -1 Illumina_reads_1.fastq \
        -2 Illumina_reads_2.fastq -it 5 --pilon pilon-1.23.jar -o pilon/
```


## Cleanup
After multiple rounds of polishing with arrow and pilon, the scaffold names tend to get messy; e.g. > Scaffold1|Arrow|Arrow|Arrow|Pilon|Pilon|Pilon|Pilon|Pilon so a custom script was used to generate simple and clean Scaffold names.

```
$ python swap_fasta_header.py Gp_D383_Scaffold pilon/round4/pilon_run4.fasta > G_pallida_D383_named.fasta
```
