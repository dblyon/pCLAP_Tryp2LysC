# Readme for pCLAP_Tryp2LysC
This is a command-line script to map tryptic peptides to Lys-C peptides.
The program expects a FASTA-file and a "peptides.txt" file from <a href="http://www.biochem.mpg.de/5111795/maxquant">MaxQuant</a> as input,
and will generate two tab delimited text files as output. The first results file, 
definable using the "-lysc" option, generates the LysC peptide sequences as well as start and 
end positions within given proteins. The second results file, definable using the "-ptm" 
option, is used as input for <a href="http://www.biochem.mpg.de/5111795/maxquant">Perseus</a>.

# Citation
Mullari M, Lyon D, Jensen LJ, Nielsen ML, Specifying RNA-Binding Regions in Proteins by Peptide Cross-Linking and Affinity Purification. Journal of Proteome Research, 2017.
DOI: 10.1021/acs.jproteome.7b00042
https://www.ncbi.nlm.nih.gov/pubmed/28648085

# Prerequisites
This script can be run on any operating system (Windows, Mac OSX, Linux) running Python.
    - Python 2.x
    - Python packages: numpy, pandas

# Usage
type "-h" or "--help" for help.
e.g. type the following command in the terminal to see options.
<code>python tryptic2LysC.py -h</code>

## example for OSX 
<code>python tryptic2LysC.py -fasta /Directory/FastaFileName.fasta -mq /Directory/MaxQuant_peptides.txt -lysc /Directory/Results_LYSC.txt -ptm /Directory/Results_PTM.txt</code>

## example for Windows
<code>python tryptic2LysC.py -fasta C:\Directory\FastaFileName.fasta -mq C:\Directory\MaxQuant_peptides.txt -lysc C:\Directory\Results_LYSC.txt -ptm C:\Directory\Results_PTM.txt</code>

# Support
feel free to contact <david.lyon@cpr.ku.dk> for help!





