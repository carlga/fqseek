# ![fqseek](docs/images/fqseek_logo_whitebg.png#gh-light-mode-only) ![fqseek](docs/images/fqseek_logo_darkbg.png#gh-dark-mode-only)

[![lang](https://img.shields.io/badge/Language-Python-yellow.svg?style=flat-square)](http://dx.doi.org/10.1093/bioinformatics/btw354)

---

`fqseek` is an integrated framework that aims to streamline the processing and analysis
of publicly available or newly generated high-throughput sequencing datasets. Its
modular architecture makes it applicable to a variety of omics technologies and by
using the snakemake workflow management system it allows to efficiently process files.
`fqseek` provides a powerful and reproducible execution environment which makes it 
particularly suitable for multi-omics integration and large-scale meta-analysis projects.

Relevant information for execution and reproducibility is stored in the configuration 
file. To create and clear the configuration file use `fqseek makeconfig` and 
`fqseek clearconfig`, respectively. `fqseek init` helps generate required index files
and databases. Retrieval of reference genomes and annotations is facilitated with
`fqseek getref`, and publicly available datasets with `fqseek getfq`. Quality control
steps with extended functionality can be accessed via `fqseek runqc`. Tailored
bioinformatics workflows can be executed with `fqseek run_rnaseq`, `fqseek run_atacseq`,
`fqseek run_scrnaseq` and `fqseek run_scatacseq`.

For more details and use case examples check the [documentation]().


## Installation

`fqseeq` can be installed via [conda](http://anaconda.org/) from the 
[bioconda](https://bioconda.github.io/) channel:

```bash
conda install -c bioconda -c conda-forge fqseek
```

or, alternatively, the more efficient:

```bash
mamba install -c bioconda -c conda-forge fqseek
```


## Usage



## Citation

If using `fqseek` please cite the [article]() and other publications relevant to your
analysis workflow.

> fqseek: an integrated framework for reproducible and efficient processing of 
> high-throughput sequencing data.


## Development & Contributions

