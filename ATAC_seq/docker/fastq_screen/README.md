# fastq-screen docker

Following installation guide [here](https://stevenwingett.github.io/FastQ-Screen/)


Once this is built, within the docker image run

```bash
fastq_screen --get_genomes
```

And save the genomes to a helpful location. Then update the `fastq_screen.conf` and keep with the docker image. Pass this using

```bash
fastq_screen --conf /path/to/fastq_screen.conf sample1.fastq sample2.fastq
```