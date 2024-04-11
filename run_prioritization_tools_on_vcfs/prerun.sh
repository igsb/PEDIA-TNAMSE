#! /bin/bash

# We also doewnload the jar file for LIRICAL version 1.3.4 from https://github.com/TheJacksonLaboratory/LIRICAL
# if nothing has changed, it should be possible to get it with

wget https://github.com/TheJacksonLaboratory/LIRICAL/releases/download/v1.3.4/LIRICAL.jar -O prio_tools/lirical/LIRICAL.jar

# We also download the exomizer command line interface 13.2.1 from https://github.com/exomiser/Exomiser/releases
# if nothing has chenged

wget https://github.com/exomiser/Exomiser/releases/download/13.2.1/exomiser-cli-13.2.1-distribution.zip -O prio_tools/exomizer/exomiser-cli-13.2.1-distribution.zip

# and uncompress it

cd prio_tools/exomizer/
unzip -o exomiser-cli-13.2.1-distribution.zip
rm exomiser-cli-13.2.1-distribution.zip
cd -

# Download the database 2302 from https://data.monarchinitiative.org/exomiser/data/index.html
# if nothing has changed

wget https://data.monarchinitiative.org/exomiser/data/2302_hg19.zip -O 2302_hg19.zip
wget https://data.monarchinitiative.org/exomiser/data/2302_phenotype.zip -O 2302_phenotype.zip

# unzip the database in prio_tools/exomizer/exomiser-cli-13.2.1/data and link the database in prio_tools/lirical/data

mkdir prio_tools/exomizer/exomiser-cli-13.2.1/data
mv 2302_hg19.zip prio_tools/exomizer/exomiser-cli-13.2.1/data
mv 2302_phenotype.zip prio_tools/exomizer/exomiser-cli-13.2.1/data
cd prio_tools/exomizer/exomiser-cli-13.2.1/data
unzip 2302_hg19.zip
unzip 2302_phenotype.zip
rm 2302_hg19.zip
rm 2302_phenotype.zip
cd -
cd prio_tools/lirical/
java -jar LIRICAL.jar download --overwrite
cd -
ln -s $PWD/prio_tools/exomizer/exomiser-cli-13.2.1/data/2302_hg19 $PWD/prio_tools/lirical/2302_hg19

# download hg19 reference and dictionary

aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta references/human_g1k_v37_decoy.fasta
aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai references/human_g1k_v37_decoy.fasta.fai
aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict references/human_g1k_v37_decoy.dict
aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta references/Homo_sapiens_assembly38.fasta
aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai references/Homo_sapiens_assembly38.fasta.fai
aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict references/Homo_sapiens_assembly38.dict

# downloaaws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fai references/human_g1k_v37_decoy.faid also gene_info.gz from https://www.ncbi.nlm.nih.gov

wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz -O references/gene_info.gz

# We download the xrare docker from https://web.stanford.edu/group/wonglab/Xrare/xrare-pub.2021.html
# if nothing has changed we should be able to get it with
# Note Goggle limit the number of anonymous download, so in case of failire log-in with your browser in google and download the file

gdown https://drive.google.com/uc?id=1wtwLiFNhbYhK30QkWgHA5l2HSrZiOuGP

# we than load the image into docker (docker must be installed)

sudo docker load -i xrare-pub.2021.tar.gz
rm xrare-pub.2021.tar.gz

