![genofunk_logo][logo]

[logo]: https://github.com/rmcolq/genofunk/blob/master/genofunk.png "Genofunk logo"
# 
Funky tools for viral consensus sequences

## Contents
* [Introduction](#introduction)
* [Quick Start](#quick-start)
* [Installation](#installation)
* [Usage](#usage)

## Introduction

## Quick Start

## Installation
Clone and enter the repository using
```
git clone https://github.com/rmcolq/genofunk.git
cd genofunk
```
then install using either
```
pip install .
```
or
```
python3 setup.py install
python3 setup.py test
```

## Usage
```
usage: genofunk <subcommand> <options>

Funky functions for virus consensus sequences

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

Available subcommands:
  
    annotate  Identify edits for a consensus sequence to remove frame shifts
              with respect to the reference. Assumes you know the closest
              reference for the provided sample and that this is the closest
              reference to all contigs in the consensus file
    merge     Collates the edits for several files and flags frame shifts
              which occur multiple times and could be real
    apply     Applies the edits to the consensus sequences and outputs an
              amino acid sequence for each
    translate
              Aligns the (corrected) nucleotide sequences based on an amino
              acid alignment

```

### Annotate
Annotate creates a list of edits to remove frameshifts between the open reading frames of the reference and the consensus sequence provided.
```
usage: genofunk annotate -r <reference_info_json> -c <consensus_file> -a ref_accession [-e <edit_file>]

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_FILE, --reference_file REFERENCE_FILE
                        Input file: a JSON formatted file containing
                        references. For each reference,the dictionary should
                        include nucleotide sequence and ORF coordinates
  -c CONSENSUS_FILE, --consensus_file CONSENSUS_FILE
                        Input file: a FASTA formatted file containing
                        nucleotide consensus sequences
  -a ACCESSION, --accession ACCESSION
                        Accession (as written in the reference_file) for the
                        closest reference sequence
  -e EDIT_FILE, --edit_file EDIT_FILE
                        CSV file containing edits already found for this
                        sample
  -v, --verbose         Run with high verbosity (debug level logging)
  ```
  
### Merge
Merge collects the lists of edits generated with annotator and combines them interactively, allowing the user to decide to keep or remove edits which occur multiple times and therefore might be real differences rather than sequencing artifacts.
```
usage: genofunk merge -d <directory> 

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Input directory: containing pairs of consensus files
                        and their edit lists (*.fasta, *.fasta.edits)
  -v, --verbose         Run with high verbosity (debug level logging)
  ```

### Apply
Apply adds the accepted edits to the consenusus sequences and outputs both nucleotide and amino acid fasta files ready for multiple sequence alignment.
```
usage: genofunk apply -d <directory> -e <edit_file> -o <output_file>

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Input directory: containing pairs of consensus files
                        and their edit lists (*.fasta, *.fasta.edits)
  -e EDIT_FILE, --edit_file EDIT_FILE
                        CSV file containing edits already found for this
                        sample
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output file name
  -v, --verbose         Run with high verbosity (debug level logging)
  ```
  
  ### Translate
  Translate uses an amino acid multiple sequence alignment to guide a nucleotide level multiple sequence alignment.
  ```
  usage: genofunk translate -c <consensus_file> -a <alignment> [ -o <output_file> ]

optional arguments:
  -h, --help            show this help message and exit
  -c CONSENSUS_FILE, --consensus_file CONSENSUS_FILE
                        Input file: a FASTA formatted file containing
                        (corrected) nucleotide consensus sequences as output
                        by `apply`
  -a ALIGNMENT_FILE, --alignment_file ALIGNMENT_FILE
                        Input file: a FASTA formatted amino acid alignment
                        file
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output file name
  -v, --verbose         Run with high verbosity (debug level logging)
  ```
