### Programs not in this repository

- canu (https://github.com/marbl/canu)
- wtdbg2 (https://github.com/ruanjue/wtdbg2)
- BUSCO (https://busco.ezlab.org/)
- minimap2 (https://github.com/lh3/minimap2)
- purge_haplotigs (https://bitbucket.org/mroachawri/purge_haplotigs/)
- finisherSC (https://github.com/kakitone/finishingTool)
- Hydraslayer/GapFiller (https://github.com/Jorisvansteenbrugge/GapFiller)
- SSPACE-Longread (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-211)
- Medaka (https://github.com/nanoporetech/medaka)
- Pilon (https://github.com/broadinstitute/pilon)


# Trimming of reads

Since we sequenced a pool of genomes we perform a read correction step to reduce the number of haplotypes in the reads. Canu is used for the error correction of the reads. Specifically, parameters these parameters were tweaked: `corOutCoverage=200 correctedErrorRate=0.15`.
The reads that we used are nanopore reads and pacbio reads.

```
$ canu -correct -p gr19 -d HS_correctedReads_Nanopore genomeSize=180m corOutCoverage=200 correctedErrorRate=0.15 \
    -pacbio-raw HS_nanoporeReads.fastq.gz
$ canu -correct -p gr22 -d HS_correctedReads_Pacbio genomeSize=180m corOutCoverage=200 correctedErrorRate=0.15 \
    -pacbio-raw HS_pacbioReads.fasta.gz
```
The corrected reads produced by canu are stored as **HS_correctedReads_Pacbio.fasta.gz** and **HS_correctedReads_Nanopore.fasta.gz**


# Assembly with wtdgb2

For each genome approximately 100 initial assemblies were produced with wtdbg2, optimizing the parameters minimal read length, k-mer size, and minimal read depth. This optimization step is performed with a custom script [optimize_wtdbg2.py](https://github.com/Jorisvansteenbrugge/GROS_genomes/blob/main/Assembly/optimize_wtdbg2.py). For the initial assembly, only nanopore reads were used.

```
$ mkdir HS_initial_assemblies
$ python Assembly/optimize_wtdbg2.py --reads HS_correctedReads_Nanopore.fasta.gz --wd HS_initial_assemblies -t 8
```

The quality of the initial assemblies was assessed based on whether the assembly size was close to the genome size estimate of appproximately 180Mb and the completeness ([BUSCO](https://github.com/Jorisvansteenbrugge/GROS_genomes/blob/main/BUSCO.md)).**Assembly_L6000_g180_X100_e8** was selected to continue.


# Purging of Haplotigs
```
$ minimap2 -t 4 -ax map-pb Assembly_L6000_g180_X100_e8.fasta HS_nanoporeReads.fastq.gz --secondary=no | samtools sort -m 30G \
    -@ 2 -o HS_Assembly_L6000_g180_X100_e8.bam
$ purge_haplotigs hist -b HS_Assembly_L6000_g180_X100_e8.bam -g Assembly_L6000_g180_X100_e8.fasta -t 8
$ purge_haplotigs cov -i HS_Assembly_L6000_g180_X100_e8.bam.genecov -l 4 -m 15 -h 110
$ purge_haplotigs purge -g Assembly_L6000_g180_X100_e8.fasta -c coverage_stats.csv
$ mv curated.fasta H_schachtii_PH.fasta
```

# Finishing Tool
```
$ mkdir finisheroutput_HS
$ ln -s HS_nanoporeReads.fasta.gz finisheroutput_HS/raw_reads.fasta.gz
$ finisherSC.py -par 8 finisheroutput_HS/ ~/tools/miniconda/bin/mummer
```

# Scaffolding SSPACE-Longread
```
$ mkdir scaffolding_HS
$ SSPACE-LongReadunix.pl -c finisheroutput_HS/improved3.fasta -p HS_nanoporeReads.fastq.gz \
     -t 8 -b scaffolding_HS/ -o 1000 -g 500 -k1
```
# Gap Filling

**Gr-Line19**
```
$ mkdir gapfilling_HS
$ minimap2 -t <threads> -Hax map-pb scaffolding_HS/scaffolds.fasta  HS_correctedReads_Nanopore.gz \
    | samtools view -hF 256 - | samtools sort -@ 2 -m 30G -o gapfilling_HS/correctedreads_onScaffolds.bam - 
$ samtools index gapfilling_HS/correctedreads_onScaffolds.bam
$ Hydraslayer -o gapfilling_HS/HS_scaffolds_gapfilled.fasta -c 10 -a 60 gapfilling_HS/scaffolds.fasta \
    gapfilling_HS/correctedreads_onScaffolds.bam
```

# Polishing with Medaka
```
medaka_consensus -i HS_nanoporeReads.fastq.gz -o H_schachti_scaffolded_gapFilled.fasta -t 12
```
