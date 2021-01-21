# DomainExtractor
Extracts A domain sequences from .fasta files containing whole polypeptide sequences

# Installation Instructions
git clone repository to desired location
```
git clone https://github.com/passm97/DomainExtractor
```
Download and install hmmer from http://hmmer.org
<br>
Download Pfam-A hmm file to same location using:
<br>
```
curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz -o Pfam-A.hmm.gz

gunzip Pfam-A.hmm

hmmpress Pfam-A.hmm
```

# Example Usage 

python DomainExtractor.py --hmmfile Pfam-A.hmm --file example_input.fasta --domainposition 2
