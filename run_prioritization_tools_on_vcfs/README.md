# TNAMSE benchmark prioritization

## Installation

Create the conda environment tnamse_revision_benchmark and activate it

```
conda env create -f environment.yml
conda activate tnamse_revision_benchmark
```

Download and uncompress the datasets

VCFs.zip

VCFs_training.zip (Optional if you want to run on the training dataset)

all_cases.zip (Optional if you want to restrict the variant list)


Run the script or follow each step in prerun.sh. If you run the script, the last steps in the script
has a good change to fail due to limitation in anonymous download on google drive, in such case please download 
the file xrare-pub.2021.tar.gz from https://drive.google.com/uc?id=1wtwLiFNhbYhK30QkWgHA5l2HSrZiOuGP and run the
following commands

```
sudo docker load -i xrare-pub.2021.tar.gz
rm xrare-pub.2021.tar.gz
```

change in the file prio_tools/exomizer/exomiser-cli-13.2.1/application.properties
the lines from

```
exomiser.hg19.data-version=2109
exomiser.phenotype.data-version=2109
```

to

```
exomiser.hg19.data-version=2302
exomiser.phenotype.data-version=2302
```

## Run

In order to run the ranking prioritization benchmark you need:

- all_cases.csv input file. The columns are case, HPO term, and gene causing disease (Note the first number is unused).

- all_cases forder is optional, but if there is a case matching, the list of variant in the xlsx will be used

Make sure docker is running, and your username can use docker without sudo or root permission

To run:

```
python prioritize.py
```



