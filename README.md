# RibDif

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

`conda install -c bioconda vsearch fasttree muscle barrnap -y`

Check your install:

`ribdif -h`

# Running

`ribdif/ribdif -g <some genus name>`

Use `-h` for full options list.






