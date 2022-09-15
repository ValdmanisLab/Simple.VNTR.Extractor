Usage: You have a folder full of whole genomes, and know a consistent starting and ending motif (or what is directly before and after) a VNTR (or other genomic element) and wish you had a tool to extract exactly what is between those two sequences in each genome. Luckily you have Simple.VNTR.Extractor!

Dependencies are NOT yet included. You will need Python and the following Python modules:\
Bio os pandas argparse

Which can be installed with pip:\
    python3 -m pip install Bio pandas argparse

Additionally for MacOS you will need the module tables, and the homebrew package c-blosc before you can install the above required packages. Probably.

Windows Example: Execute this in terminal from the directory you want the output to be in:\
  C:/Users/username/AppData/Local/Programs/Python/Python310/python.exe C:/pythontutorial/Script1.0.py -f E:\wholegenomes -l 30 -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC

Mac Example: From the folder both containing the script and that you want the output to be in:\
  Python3 Script1.0.py -f /Volumes/wholegenomes -l 30 -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC

Outputs:\
out.csv                       – excel file with each allele named after the file, cut into simple motifs based on the length provided from the beginning of the allele \
alleles.txt                   – text file with each allele after the name of the file.\
errors_and_frequency.txt      – text file with any errors (ie. ‘Filename has only reverse end’ or ‘filename has no sign of tandem repeat in forward, reverse’) and ranked frequency of simple motifs.


Current help from stable build:\
Extract VNTR

options:\
  -h, --help          show this help message and exit\
  -f , --folder       Path to folder. ex: E:\HPRC\
  -l , --length       length of consensus motif. ex: 30\
  -s , --start        beginning of VNTR. ex: GCTA\
  -e , --end          end of VNTR. ex: GCTA\
  -b [], --before []  Is the start value of the VNTR before the VNTR? y OR n\
  -a [], --after []   Is the end value of the VNTR after the VNTR? y OR n\
  -t [], --type []    Filetype ex: .fasta (as default is .fa)


Current help from testing build:\
Extract VNTR

options:\
  -h, --help          show this help message and exit\
  -f , --folder       Path to folder. ex: E:\HPRC\
  -l , --length       length of consensus motif. ex: 30\
  -s , --start        beginning of VNTR. ex: GCTA\
  -e , --end          end of VNTR. ex: GCTA\
  -b [], --before []  Is the start value of the VNTR before the VNTR? y OR n\
  -a [], --after []   Is the end value of the VNTR after the VNTR? y OR n\
  -t [], --type []    Filetype ex: .fasta (as default is .fa)\
  -n [], --name []    output name)
