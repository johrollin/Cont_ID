<p align="center">
    <img src="https://github.com/johrollin/viral_contamination/blob/master/img/cont_id_logo_and_names.png" alt="Logo" width="600">
</p>
<p align="center">
        <a href="https://github.com/johrollin/Cont_ID/releases">
            <img src="https://img.shields.io/github/release/johrollin/Cont_ID.svg" />
        </a>
        <a href="https://github.com/johrollin/Cont_ID/blob/master/LICENSE" alt="License">
            <img src="https://img.shields.io/badge/License-GNUv3-purple.svg">
        </a>
        <a href="https://github.com/johrollin/Cont_ID/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/johrollin/Cont_ID">
        </a>     
        <br/>           
        <a href="https://www.biorxiv.org/content/XX">
          <img src="https://zenodo.org/badge/DOI/XX/XX.svg">
        </a>
    
 </p>

# About The Project
Cont-ID is a method designed to check for cross-contamination in viruses previously identified in metagenomic datasets. It relies on a simple principle, every sample in a sequencing batch should have been processed the same way with at least one alien control. It uses a voting system to classify every species prediction on every sample of the sequencing batch into (true) infection or (cross) contamination. This tool helps the virologist confirm viral detection in HTS data. For more information: https://www.biorxiv.org/content/XX 

# Getting Started

## Prerequisites

The code has been tested with Python V3.7.4 with the following dependancies:

```shell
pandas
numpy
```

## Installation

```shell
git clone https://github.com/johrollin/Cont_ID.git

cd Cont_ID/

python classify_contamination_vote.py -h
```


## Usage

#### Input file template

Cont-ID works with .csv input and accpet both ";" and "," as separator.

Template control & data file 

#### How to launch

```shell
python classify_contamination_vote.py -h

arguments:
  -h, --help            show this help message and exit
  -dr DATA_REPOSITORY, --data_repository DATA_REPOSITORY
                        repository where the data are
  -all ALL_FILES, --all_files ALL_FILES
                        If True, will go through every file of the
                        data_repository (default False) NOT COMPATIBLE with
                        file_data and file_control
  -s STANDARDISATION, --standardisation STANDARDISATION
                        expected number of reads for each sample to use for
                        standardisation default 5000000
  -tp THRESHOLD_PERSONALISATION [THRESHOLD_PERSONALISATION ...], --threshold_personalisation THRESHOLD_PERSONALISATION [THRESHOLD_PERSONALISATION ...]
                        personalisation threshold like "2:1000:1.5:5:1"
                        "0.002:500:1.5:5:1" with ":" as separator (the 5
                        control limit are mandatory)
  -fd FILE_DATA, --file_data FILE_DATA
                        file name containing data for processing NOT
                        COMPATIBLE with all_files
  -fc FILE_CONTROL, --file_control FILE_CONTROL
                        file name containing alien control data NOT COMPATIBLE
                        with all_files
```
There is two mode for Cont-ID:

- If you have few sequencing batch to check, you can launch Cont-ID individually for each of them

```shell
python classify_contamination_vote.py -dr /mnt/c/bioinfo/virus_test/ -fd Input_file_batch1 -fc control_file_batch1

```
The '-s' and '-tp' option can be used if you don't want to use the default parameter.

The '-all' option should never be used if '-fd' and '-fc' are given as '-all' asks the process of all the files while '-fd' and '-fc' give one specific to process.

- If you have several batch to check, you can do them all in one run:

The file containing the data should be named like "Input_file_XXX" with XXX as unique identifier (between data file). The control file should be name like "Control_file_XXX"  with XXX identical to the related datafile.

 ```shell
python classify_contamination_vote.py -dr /mnt/c/bioinfo/virus_test/ -all True 
```
The tool will process all the files present in the '-dr folder' when they follow the naming procedure describe above.

The '-s' and '-tp' option can be used if you don't want to use the default parameter.

The '-all' option should never be used if '-fd' and '-fc' are given as '-all' asks the process of all the files while '-fd' and '-fc' give one specific to process.

The '-tp' THRESHOLD_PERSONALISATION option needs two sets of threshold like "0.002:500:1.5:5:1" with ":" as separator each set need 5 numbers that will be used as threshold for each rule

<p align="center">
    <img src="https://github.com/johrollin/viral_contamination/blob/master/img/Cont-ID_formula_casesV2.png" alt="Logo">
</p>
For more information, you can read the artcile (link at the begining of this README)

## How to read the results

#### Result example

#### Advance statistic 

small explanation to the further analysis scripts

## Contributing



## Contact

Mail + Orc_ID 
</a>
<a href="https://twitter.com/intent/follow?screen_name=johan_rollin" alt="Author Twitter">
    <img src="https://img.shields.io/twitter/follow/johan_rollin?style=social&logo=twitter"
        alt="follow on Twitter">
</a>
<a href="https://twitter.com/intent/follow?screen_name=Be_Phytopath" alt="Author Twitter">
    <img src="https://img.shields.io/twitter/follow/Be_Phytopath?style=social&logo=twitter"
        alt="follow on Twitter">
</a>       
