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
The reads that we used are nanopore reads and pacbio reads.

```
$ canu -correct -p gr19 -d HS_correctedReads_Nanopore genomeSize=180m corOutCoverage=200 correctedErrorRate=0.15 \
    -pacbio-raw HS_nanoporeReads.fastq.gz
$ canu -correct -p gr22 -d HS_correctedReads_Pacbio genomeSize=180m corOutCoverage=200 correctedErrorRate=0.15 \
    -pacbio-raw HS_pacbioReads.fasta.gz
```
The corrected reads produced by canu are stored as **HS_correctedReads_Pacbio.fasta.gz** and **HS_correctedReads_Nanopore.gz**
