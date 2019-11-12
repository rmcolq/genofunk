![genofunk_logo][logo]

[logo]: https://github.com/rmcolq/genofunk/blob/master/genofunk.png "Genofunk logo"
# genofunk
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

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

Available subcommands:
  
    annotator
              Identify edits for a consensus sequence to remove frame shifts
              with respect to the reference. Assumes you know the closest
              reference for the provided sample and that this is the closest
              reference to all contigs in the consensus file
    gather    Collates the edits for several files and flags frame shifts
              which occur multiple times and could be real
```

### Annotator
Annotator creates a list of edits to remove frameshifts between the open reading frames of the reference and the consensus sequence provided.
```
usage: genofunk annotator -r <reference_info_json> -f <consensus_fasta> -a ref_accession [-e <edit_file>]

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_FILE, --reference_file REFERENCE_FILE
                        Input file: a JSON formatted file containing
                        references. For each reference,the dictionary should
                        include nucleotide sequence and ORF coordinates
  -f CONSENSUS_FILE, --consensus_file CONSENSUS_FILE
                        Input file: a FASTA formatted file containing the
                        consensus sequence(s) for a sample
  -a ACCESSION, --accession ACCESSION
                        Accession (as written in the reference_file) for the
                        closest reference sequence
  -e EDIT_FILE, --edit_file EDIT_FILE
                        CSV file containing edits already found for this
                        sample
  -v, --verbose         Run with high verbosity (debug level logging)
  ```
  
### Gather
Gather collects the lists of edits generated with annotator and combines them interactively, allowing the user to decide to keep or remove edits which occur multiple times and therefore might be real differences rather than sequencing artifacts.
```
usage: genofunk gather -d <directory> 

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Input directory: containing pairs of consensus files
                        and their edit lists (*.fasta, *.fasta.edits)
  -v, --verbose         Run with high verbosity (debug level logging)
  ```
