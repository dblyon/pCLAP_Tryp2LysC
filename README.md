# Readme for pCLAP_Tryp2LysC
This is a proteomics tool to map tryptic peptides to Lys-C peptides

## Prerequisites
# Python 2.x
# Python packages: numpy, pandas

## Usage
# type "-h" or "--help" for help.
# e.g. type the following command in the terminal/command-line 
python tryptic2LysC.py -h
# example for OSX 
python tryptic2LysC.py -fasta /Directory/FastaFileName.fasta -mq /Directory/MaxQuant_peptides.txt -lysc /Directory/Results_LYSC.txt -ptm /Directory/Results_PTM.txt

# example for Windows
python tryptic2LysC.py -fasta C:\Directory\FastaFileName.fasta -mq C:\Directory\MaxQuant_peptides.txt -lysc C:\Directory\Results_LYSC.txt -ptm C:\Directory\Results_PTM.txt
