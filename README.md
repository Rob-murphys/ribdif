# RibDif2

RibDif evaluates the usefulness of a given amplicon at differentiating within a genus or species

# Installation

## Setup conda

For this to work, you need [miniconda](https://docs.conda.io/en/latest/miniconda.html) to be installed on your system (in ubuntu if you are using WSL)

```
conda create --name ribdif2 python=3.11 -y
conda activate ribdif2
```

Ribdif is most easily installed with pip inside a virtual environment. Make sure your pip version is up to date!


`pip install git+https://github.com/Rob-murphys/ribdif.git`

Next we need to install some dependencies:

`conda install -c bioconda vsearch fasttree mafft barrnap pyani -y`

Check your install:

`ribdif -h`

# Usage

## Running with default primers

`ribdif -g <some genus name>`

Use `-h` for full options list.

An example run with Ruegeria as the genus:

`ribdif -g Ruegeria`

You can also run with a species aswell:

`ribdig -g "Mycoplasma bovis"`

## Running with custom primers

`ribdif -g <some genus name> -p /path/to/my_primer_file.txt`

### Primer file format

```
v3v4	CCTACGGGNGGCNGCAG	GACTACNNGGGTATCTAATCC	550
v1v9	AGRGTTYGATYMTGGCTCAG 	RGYTACCTTGTTACGACTT	1800
```
Where the columns are: primer name | forward sequence | reverse sequence | expected amplicon size

## Output

RibDif generate a new directory within the specified output directory (<current working directory>/results by default) names after the genus in question.

```
Pseudoalteromonas
    ├── amplicons
        ├── <genus>-<primer name>-meta.tsv          # information about each leaf tip in the amplicon tree
        ├── <genus>-<primer name>.aln               # aligned amplicons from a given primer
        ├── <genus>-<primer name>.amplicons         # the amplicons generated by a given primer
        ├── <genus>-<primer name>.summary           # a summary of each amplicon generated (ID, originating sequence, position and length)
        ├── <genus>-<primer name>.tree              # tree of alligned amplicons for a given primer
        ├── <genus>-<primer name>.uc                # cluster file of amlicons for a given primer
        ├── <genus>-<primer name>_confusion.csv     # matrix of which genomes we can tell apart with these primers (a heatmap is made of this)
    │   └── <primer name>-clusters                  # directory of clusters generated from the amplicons of a givem primer
    ├── figures
    │   ├── <genus>-<primer name>_graphs.pdf    # visual network of which genomes and thereby species can be differentiated
    │   └── <genus>-<primer name>_heatmaps.pdf  # page 1 heatmap shows which genomes belong to which allele clusters and page 2 is the confusion matrix
    ├── full
    │   ├── <genus>.16S
    │   ├── <genus>.16sAln
    │   └── <genus>.16sTree
    ├── refseq
    │   └── bacteria
    │       ├── GCF_xxxxxxxxx.x
    │       │   ├── GCF_xxxxxxxxx.x_ASMxxxxv1_genomic.fna               # the original genome file, headers are renamed for downstream compatibility
    │       │   ├── GCF_xxxxxxxxx.x_ASMxxxxv1_genomic.fna.fai           # indexed fna file
    │       │   ├── GCF_xxxxxxxxx.x_ASMxxxxv1_genomic.fna.gz            # the gzipped original genome file
    │       │   ├── GCF_xxxxxxxxx.x_ASMxxxxv1_genomic.fna.rRNA          # rRNA file generated by barrnap
    │       │   ├── GCF_xxxxxxxxx.x_ASMxxxxv1_genomic.fna.rRNA.16S      # 16S genes taked from the .rRNA file
    │       │   ├── GCF_xxxxxxxxx.x_ASMxxxxv1_genomic.fna.rRNA.16sAln   # alligned .16S file (by mafft)
    │       │   ├── indiv_16S_dir                                       # folder for individual 16S filed if --ani is enabled
    │       │   ├──  ani                                                # files generated from pyani if --ani is enable
    │       │   └── MD5SUMS
    │       │   .
    │       │   .
    │       │   .
    │       └── GCF_xxxxxxxxx.x
    └── ribdif_logs
```
The most important files will be within the `Pseudoalteromonas` and then the `figures/` subdirectory.

If running as default on whole 16S genes extracted by barrnap you will see:
`<genus>_16S_summary.tsv` which contains a summery of each genomes and how many whole 16S genes were present
`<genus>_<primer name>-amp_summary.tsv` which contains a summery of each genomes and how many amplicons it produced

__OR__ if run with `--whole-genome` mode you will only see:
`<genus>_<primer name>-amp_summary.tsv` which contains a summery of each genomes and how many amplicons it produced

No matter your settings you will have `<genus>_<primer name>_overlap_report.txt` which show the text based summary of your run that is also printed to console

Within `figures/` are the heatmaps and network graphs for a more visual representation of the results.

`full/` contains the concatinated full 16S sequences if not using `--whole-genome`.

The individually downloaded genomes are found in `refseq/<domain>`







