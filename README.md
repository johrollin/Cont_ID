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
Cont-ID is a method designed to check for cross-contamination in viruses previously identified in metagenomic datasets. It relies on a simple principle, every sample in a sequencing batch should have been processed the same way with at least one alien control. It uses a voting system to classify every species prediction on every sample of the sequencing batch into (true) infection or (cross) contamination. This tool helps the virologist confirm viral detection in HTS data.

For more information: https://www.biorxiv.org/content/XX 

# Getting Started

## Prerequisites

The code has been tested with Python V3.7.4 with the following dependencies:

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

Cont-ID works with .csv input and accepts both ";" and "," as separator.

The input template file for control & data are csv files available in this repository, and should be completed as follows:
<img src="https://github.com/johrollin/viral_contamination/blob/master/img/input_file_column.PNG" alt="Logo" >

- **Virus detected** : Name (or abbreviation) of virus
- **Sample name** : Name of the sample (free text but no special character allowed)
- **Sample ID** : ID of the sample (the couple virus detected + sample ID need to be unique)
- **Mapped reads Nr.** : Reads number that map the specific virus on the specific sample.
- **Mapped reads deduplication...** : percentage of reads removed by deduplication (see article) can be fill with 'ND' (No Data) if not calculated or below the deduplication limit.
- **Total_Reads_Nr**: total reads number of the specific sample.

The control file should contain the same columns with an additional indexing one that can be filled with 'contamination' or 'infection'.

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
  -fd FILE_DATA, --file data FILE_DATA
                        file name containing data for processing NOT
                        COMPATIBLE with all files
  -fc FILE_CONTROL, --file control FILE_CONTROL
                        file name containing alien control data NOT COMPATIBLE
                        with all files
```
There are two modes for Cont-ID:

- If you have few sequencing batch to check, you can launch Cont-ID individually for each of them:

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

The '-tp' THRESHOLD_PERSONALISATION option needs two sets of thresholds like "0.002:500:1.5:5:1" with ":" as separator each set needs 5 numbers that will be used as threshold for each rule.

<p align="center">
    <img src="https://github.com/johrollin/viral_contamination/blob/master/img/Cont-ID_formula_casesV2.png" alt="Logo">
</p>
For more information, you can read the article (link at the beginning of this README).

## How to read the results

Results are displayed in csv format. In the results file you will find all the data present in the input file with the addition of several columns, the first three (standardize_reads_nb_mapped, mapping_ratio, standardize_mapping_ratio) contain metrics data after standardization. The next columns (Standard  classification (3 votes) (case 1/2), Total classification (2 votes) (case 1/2), Classification (case 1/2), Comment (case 1/2)) contains the tool prediction for each case, the comments provide the detail of the rules vote. Finally, the last column (Comparison both case) gives the prediction according to the comparison between the 2 cases with 3 possible outputs (contamination, infection or unconfirmed).

This result file provides a lot of details on the decision making of the tool as the user can decide to rely more or less on some rules/case for manual confirmation.

#### Advanced statistic 

In the further analysis repository, there is a script that can help you automatically obtain statistic? they are made available even if they are not formally part of Cont-ID. To use them you have to change manually the path in 'if __name__ == "__main__":'

The script compare_predictionv2_vote.py should be used first as it calculates on each result file additional statistic if you provide the indexing status. The indexing status can be given by using another file that contain the same information as the input file with the addition of the indexing status (the last column called 'Indexing' that contains only 'contamination' or 'infection'). The statistic obtained give you the amount of correct/unknown/wrong prediction for each case/vote level, the confusion matrix (TP, TN, FP, FN as explain in the article) and number like: Diagnostic sensitivity (DSE), Diagnostic specificity (DSP), accuracy, False Discovery rate (FDR), or False Omission Rate (FOR). 

The script Global_summary_table.py aim to give the previous statistic but between several results that may be connected. In our use case we use several banana sequencing, that script allow the comparison between the different sequencing batch to try to decide if the threshold used to make the prediction seems the most accurate for banana related data.

These scripts are here to help the user to test if the tool's option used work on their specific training datasets.

## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are greatly appreciated.

If you have a suggestion that would make this better, please fork the repo and create a pull request. Don't forget to give the project a star! Thanks again!

Special thanks to <a href="https://github.com/ValentinRollin">Valentin Rollin </a> for the help provided for the further_analysis scripts.


This work was part of the <a href="https://inextvir.eu/Home/WhatIsInextvir">INEXTVIR</a> project.

This project has received funding from the European Union’s Horizon2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement n° 813542.



## Contact

<a href="https://twitter.com/intent/follow?screen_name=johan_rollin" alt="Author Twitter">
    <img src="https://img.shields.io/twitter/follow/johan_rollin?style=social&logo=twitter"
        alt="follow on Twitter">
</a>
<a href="https://twitter.com/intent/follow?screen_name=Be_Phytopath" alt="Author Twitter">
    <img src="https://img.shields.io/twitter/follow/Be_Phytopath?style=social&logo=twitter"
        alt="follow on Twitter">
</a>       
  
