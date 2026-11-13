# E. Coli
```bash
curl -sLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz
curl -sLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

Sample paired read for this genome:
```bash
curl -sLO https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR33366732/SRR33366732
```

You need to use the SRA toolkit to extract the actual `.fastq` files from this.
