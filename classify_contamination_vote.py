#!/usr/bin/env python
"""Script to run cross-contamination analysis."""

import argparse
import os
import pandas as pd
import numpy as np
import statistics

def open_file(out_dir, file_name_data, file_name_control, col_name, col_name_control):
    """
    """
    try:
        virus_data = pd.read_csv(os.path.join(out_dir, file_name_data),
                                sep=";", header=0, names=col_name,
                                index_col=False)
    except pd.errors.ParserError:
        virus_data = pd.read_csv(os.path.join(out_dir, file_name_data),
                                sep=",", header=0, names=col_name,
                                index_col=False)        
    try:
        control_data = pd.read_csv(os.path.join(out_dir, file_name_control),
                        sep=";", header=0, names=col_name_control,
                        index_col=False)
    except pd.errors.ParserError:
        control_data = pd.read_csv(os.path.join(out_dir, file_name_control),
                        sep=",", header=0, names=col_name_control,
                        index_col=False)

    virus_data = virus_data.drop(virus_data.index[pd.isna(virus_data["Virus_detected"])])
    control_data = control_data.drop(control_data.index[pd.isna(control_data["Virus_detected"])])
    return virus_data, control_data

def calculate_mapping_ratio(virus_data,control_data, col_name, standardisation):
    """
    Calculate Mapped reads Nr./ Hihgest Reads Nr. Per sample in the same batch (for same virus)

    """
    list_mapping_ratio = []
    list_mapping_ratio_control = []
    standardize_list_mapping_ratio_control = []
    standardize_list_mapping_ratio = []
    standardize_list_reads_nb_mapped = []
    max_read_nb = max(control_data["Reads_nb_mapped"])
    row = control_data.iloc[control_data['Reads_nb_mapped'].idxmax()]
    total_read_nb_ref = row["Total_Reads_Nr"]
    standardisation = int(standardisation)

    standardize_max_read_nb = (max_read_nb*standardisation)/total_read_nb_ref

    #calculate mapping ratio for control (and standardize) 
    for element in control_data.itertuples():
        native_reads_nb_mapped = element.Reads_nb_mapped
        ratio = native_reads_nb_mapped/max_read_nb
        list_mapping_ratio_control.append(ratio)

        standardize_reads_nb_mapped = (native_reads_nb_mapped*standardisation)/element.Total_Reads_Nr
        standardize_ratio= standardize_reads_nb_mapped/standardize_max_read_nb
        standardize_list_mapping_ratio_control.append(standardize_ratio)

    control_data["mapping_ratio"] = list_mapping_ratio_control
    control_data["standardize_mapping_ratio"] = standardize_list_mapping_ratio_control


    #calculate mapping ratio for data (and standardize)
    for el in virus_data.itertuples():
        standardize_reads_nb_mapped = (el.Reads_nb_mapped*standardisation)/el.Total_Reads_Nr
        standardize_list_reads_nb_mapped.append(standardize_reads_nb_mapped)
    virus_data["standardize_reads_nb_mapped"] = standardize_list_reads_nb_mapped
    standardize_serie_virus_read = virus_data.groupby("Virus_detected")['standardize_reads_nb_mapped'].max()
    
    serie_virus_read = virus_data.groupby("Virus_detected")['Reads_nb_mapped'].max()

    for element in virus_data.itertuples():
        max_read_nb = serie_virus_read.loc[element.Virus_detected]
        native_reads_nb_mapped = element.Reads_nb_mapped
        ratio = native_reads_nb_mapped/max_read_nb
        list_mapping_ratio.append(ratio)

        standardize_max_read_nb = standardize_serie_virus_read.loc[element.Virus_detected]
        standardize_native_reads_nb_mapped = element.standardize_reads_nb_mapped
        standardize_ratio = standardize_native_reads_nb_mapped/standardize_max_read_nb
        standardize_list_mapping_ratio.append(standardize_ratio)        

    virus_data["mapping_ratio"] = list_mapping_ratio
    virus_data["standardize_mapping_ratio"] = standardize_list_mapping_ratio

    return virus_data, control_data

def calculate_standardize_threshold(control_data, case, standardisation):
    """
    calculate standardized threshold value from control for step 1
    """
    
    current_case = threshold_case[case]
    current_divider = current_case.split(":")
    ###### T1
    try:
        mean_conta = control_data.groupby("Indexing")['standardize_mapping_ratio'].mean().loc["ABSENT"]
    except KeyError:
        print("No contamination found in external control, you can be confident that you don't have cross-contamination")
        exit(0)
    # std standart deviation of contaminated sample from control virus
    try:
        std_conta = control_data.groupby("Indexing")['standardize_mapping_ratio'].std().loc["ABSENT"]
    except KeyError:
        print("No contamination found in external control, you can be confident that you don't have cross-contamination")
        exit(0)

    standardize_t1_threshold = (mean_conta+3*std_conta)/float(current_divider[0])
    if standardize_t1_threshold > 1 :
        standardize_t1_threshold = 1
    ######

    ###### T2
    max_read_nb = max(control_data["Reads_nb_mapped"])
    row = control_data.iloc[control_data['Reads_nb_mapped'].idxmax()]
    total_read_nb_ref = row["Total_Reads_Nr"]
    standardize_max_read_nb = (max_read_nb*int(standardisation))/total_read_nb_ref
    standardize_t2_threshold = int(standardize_max_read_nb/float(current_divider[1]))
    ######

    return standardize_t1_threshold, standardize_t2_threshold

def calculate_threshold(control_data, case, standardisation):
    """
    calculate threshold value from control
    Get all samplename from control to set target virus to contamination for those sample
    """
    

    control_name = []
    deduplication_ratio = []
    for el in control_data.itertuples():
        if el.Indexing == 1 or el.Indexing == "PRESENT":
            control_name.append(el.Sample_ID)
        try:
            if not pd.isna(el.deduplication):
                deduplication_ratio.append(float(el.deduplication))
        except ValueError:
            pass

    current_case = threshold_case[case]
    current_divider = current_case.split(":")

    standardize_t1_threshold, standardize_t2_threshold = calculate_standardize_threshold(control_data, case, standardisation)

    ###### T3
    # [de-duplication rate > X]
    # threshold is X
    # t3_threshold = 75.00
    # deduce X from control data
    # average control/threshold_case
    t3_threshold = (statistics.mean(deduplication_ratio))/float(current_divider[2])
    ######

    ## refinement threshold
    # conta 5
    nb_read_limit_conta = int(current_divider[3])
    # infection 1
    mapping_highest_ratio = int(current_divider[4])
    
    return standardize_t1_threshold, standardize_t2_threshold, t3_threshold, nb_read_limit_conta, mapping_highest_ratio, control_name

def check_comment(list_classification_comment, element, control_name, is_infection, is_contamination, t1_threshold):
    """
    """
    #col_name = ["Virus_detected","Sample_name","Sample ID","Reads_nb_mapped", "deduplication"]
    msg = "" 
    if element.Sample_ID in control_name and is_infection:
        msg += "The sample use for control has another virus labeled as infectious (consider using external control or changing threshold) "
    if element.Reads_nb_mapped <= 5 and is_infection:
        msg += "This infectious label is with less than 5 reads, the result validity is doubtful "
    if element.mapping_ratio == 1 and is_contamination:
        msg += "Contamination with a mapping ratio of 1, the result validity is doubtful "

    list_classification_comment.append(msg)
    return list_classification_comment

def classify_virus(virus_data, control_data, t1_threshold, t2_threshold, t3_threshold, 
        nb_read_limit_conta, mapping_highest_ratio, control_name):
    """
    """
    def reset_bool():
        """
        return boolean to characterize if infection (default: False)
        """
        return False, False, False, False, False, False, False

    list_classification_3vote = []
    list_classification_2vote = []
    list_classification_overall = []
    list_classification_comment = []

    for element in virus_data.itertuples():
        count_infect_vote = 0
        count_conta_vote = 0
        comment = ""
        r1_bool, r2_bool, r3_bool_infection, r3_bool_contamination, r3_nd,  r4_bool_contamination, r4_bool_infection = reset_bool()

        #R1
        #True => infection
        #False => contamination
        if element.standardize_mapping_ratio >= t1_threshold:
            r1_bool=True
            count_infect_vote+=1
            comment+="R1:infection, "
        else: 
            count_conta_vote+=1
            comment+="R1:contamination, "
        #R2
        if element.standardize_reads_nb_mapped >= t2_threshold:
            r2_bool=True
            count_infect_vote+=1
            comment+="R2:infection, "
        else: 
            count_conta_vote+=1
            comment+="R2:contamination, "
        #R3
        try:
            deduplication_ratio = float(element.deduplication)
            if deduplication_ratio <= t3_threshold:
                r3_bool_infection=True
                count_infect_vote+=1
                comment+="R3:infection"
            else:
                r3_bool_contamination=True
                count_conta_vote+=1
                comment+="R3:contamination"

        except ValueError: # RF (reference) or ND (Not enough Data)
            if element.mapping_ratio == mapping_highest_ratio:
                r4_bool_infection = True
                count_infect_vote+=1
                comment+="R3:infection"
            elif element.Reads_nb_mapped <= nb_read_limit_conta:
                r4_bool_contamination = True
                count_conta_vote+=1
                comment+="R3:contamination"
            else:
                print("unexpected result: Uncertain for R3-4")
            pass

        if count_infect_vote == 3 :
            list_classification_3vote.append("infection")
            list_classification_overall.append("infection")
        elif count_conta_vote == 3 :
            list_classification_3vote.append("contamination")
            list_classification_overall.append("contamination")
        else:
            list_classification_3vote.append("-")        


        if count_infect_vote == 2 :
            list_classification_2vote.append("infection")
            list_classification_overall.append("infection")
        elif count_conta_vote == 2 :
            list_classification_2vote.append("contamination")
            list_classification_overall.append("contamination")
        else:
            list_classification_2vote.append("-")

        if count_infect_vote == 1 and count_infect_vote >= count_conta_vote: # only ND for deduplication
            list_classification_overall.append("infection")
        elif count_conta_vote == 1 and count_conta_vote >= count_infect_vote: # only ND for deduplication:
            list_classification_overall.append("contamination")            
        
        list_classification_comment.append(comment)
   
    virus_data["Standart classification (3 votes)"] = list_classification_3vote
    virus_data["Total classification (2 votes)"] = list_classification_2vote
    virus_data["Comment"] = list_classification_comment
    virus_data["Overall classification"] = list_classification_overall

    

    return virus_data

def compare_virus_data(virus_data):
    """ compare case 1 and 2 results at step 3 
    """
    # for el in virus_data.itertuples():
    #     if el["Classification step3 (case 2)"] 
    virus_data["Comparison both case"] = np.where(virus_data["Classification (case 1)"] == virus_data["Classification (case 2)"] \
        , virus_data["Classification (case 1)"], "unconfirmed")

    return virus_data


def write_result(out_dir, file_name_data, case1_virus_data, case1_standardize_t1_threshold, \
        case1_standardize_t2_threshold, case1_t3_threshold, case1_nb_read_limit_conta, \
        case1_mapping_highest_ratio, virus_data, standardize_t1_threshold, standardize_t2_threshold, \
        t3_threshold, nb_read_limit_conta, mapping_highest_ratio, standardisation):
    """
    """
    path = file_name_data.split("/")
    result_folder_name = os.path.join(out_dir, "Result")

    if len(path)>1: # control/Input_file_control_batch6.csv
        core_name = str(path[1])
    else: #Input_file_control_batch6.csv
        core_name = str(path[0])
    main_name = core_name.split(".")
    new_core_name = main_name[0] +  "_vote." + main_name[1]

    file_name = "Result_" + new_core_name

    try:
        os.mkdir(result_folder_name)
    except FileExistsError:
        pass
    # "Confident classification (4 votes)": "Confident classification (4 votes) (case 1)", \
    case1_virus_data.rename(columns={ "Standart classification (3 votes)": "Standart classification (3 votes) (case 1)", \
        "Total classification (2 votes)": "Total classification (2 votes) (case 1)", \
        "Comment": "Comment (case 1)", \
        "Overall classification": "Classification (case 1)"}, inplace=True, copy=False)

    # case1_virus_data["Confident classification (4 votes) (case 2)"] = virus_data["Confident classification (4 votes)"].to_numpy()
    case1_virus_data["Standart classification (3 votes) (case 2)"] = virus_data["Standart classification (3 votes)"].to_numpy()
    case1_virus_data["Total classification (2 votes) (case 2)"] = virus_data["Total classification (2 votes)"].to_numpy()
    case1_virus_data["Comment (case 2)"] = virus_data["Comment"].to_numpy()
    case1_virus_data["Classification (case 2)"] = virus_data["Overall classification"].to_numpy()
    # compare case 1 and 2 results at step 3 
    final_virus_data = compare_virus_data(case1_virus_data)
    case1_virus_data.to_csv(os.path.join(result_folder_name, file_name), index=False, sep=';', encoding='utf-8')
    

    # write metric file
    current_divider1 = threshold_case[0].split(":")
    current_divider2 = threshold_case[1].split(":")
    file_name_metric =  "Metric_" + core_name
    with open(os.path.join(result_folder_name, file_name_metric), 'w') as out_file:
        
        out_file.write("Case 1 \n")
        if case1_standardize_t1_threshold == 1 :
            out_file.write("# mapping ratio threshold (r1) WARNING not enough contamination in your control for this metric; ((avg+3*std)/" \
                + current_divider1[0] +  ";" + str(case1_standardize_t1_threshold) + "\n")
        else:
            out_file.write("# mapping ratio threshold (r1); ((avg+3*std)/" \
                + current_divider1[0] + ";" + str(case1_standardize_t1_threshold) + "\n")
        out_file.write("# reads_nb_mapped threshold (r2); Hihgest Reads Nr. \
             Per sample in the same run/" + current_divider1[1] + ";" + str(case1_standardize_t2_threshold) + "\n")
        out_file.write("# deduplication threshold (r3); ((avg(deduplication_ratio))/"\
             + current_divider1[2] + ";" + str(case1_t3_threshold) + "\n")
        out_file.write("# nb_read_limit_conta (r4_1); read_number ;" + str(case1_nb_read_limit_conta) + "\n")
        out_file.write("# mapping_highest_ratio (r4_2); highest mapping ratio ;" + str(case1_mapping_highest_ratio) + "\n")
        out_file.write("\n")

        out_file.write("Case 2 \n")
        if standardize_t1_threshold == 1 :
            out_file.write("# mapping ratio threshold (r1) WARNING not enough contamination in your control for this metric; ((avg+3*std)/" \
                + current_divider2[0] +  ";" + str(standardize_t1_threshold) + "\n")
        else:
            out_file.write("# mapping ratio threshold (r1); ((avg+3*std)/" \
                + current_divider2[0] + ";" + str(standardize_t1_threshold) + "\n")
        out_file.write("# reads_nb_mapped threshold (r2); Hihgest Reads Nr. \
             Per sample in the same run/" + current_divider2[1] + ";" + str(standardize_t2_threshold) + "\n")
        out_file.write("# deduplication threshold (r3); ((avg(deduplication_ratio))/"\
             + current_divider1[2] + ";" + str(t3_threshold) + "\n")
        out_file.write("# nb_read_limit_conta (r4_1); read_number ;" + str(nb_read_limit_conta) + "\n")
        out_file.write("# mapping_highest_ratio (r4_2); highest mapping ratio ;" + str(mapping_highest_ratio) + "\n")
        out_file.write("\n")     
        out_file.write("# standardisation number ;" + str(standardisation) + "\n")
           



def run_analysis(out_dir, file_name_data, file_name_control, col_name, col_name_control, standardisation):
    """
    """
    # open input file
    virus_data, control_data = open_file(out_dir, file_name_data, file_name_control, col_name, col_name_control)
    # add mapping ratio data
    virus_data, control_data = calculate_mapping_ratio(virus_data,control_data, col_name, standardisation)

    # calculate threshold value from control CASE1
    standardize_t1_threshold, standardize_t2_threshold, t3_threshold, \
        nb_read_limit_conta, mapping_highest_ratio, control_name = calculate_threshold(control_data, 0, standardisation)
    # make virus classification CASE1
    virus_data = classify_virus(virus_data, control_data, standardize_t1_threshold, standardize_t2_threshold, t3_threshold, 
        nb_read_limit_conta, mapping_highest_ratio, control_name)

    # store values
    case1_virus_data = virus_data.copy()
    case1_standardize_t1_threshold = standardize_t1_threshold
    case1_standardize_t2_threshold = standardize_t2_threshold
    case1_t3_threshold = t3_threshold
    case1_nb_read_limit_conta = nb_read_limit_conta
    case1_mapping_highest_ratio = mapping_highest_ratio
    virus_data
    # calculate threshold value from control CASE2
    standardize_t1_threshold, standardize_t2_threshold, t3_threshold, \
        nb_read_limit_conta, mapping_highest_ratio, control_name = calculate_threshold(control_data, 1, standardisation)
    # make virus classification CASE2
    virus_data = classify_virus(virus_data, control_data, standardize_t1_threshold, standardize_t2_threshold, t3_threshold, 
        nb_read_limit_conta, mapping_highest_ratio, control_name)

    write_result(out_dir, file_name_data, case1_virus_data, case1_standardize_t1_threshold, \
        case1_standardize_t2_threshold, case1_t3_threshold, case1_nb_read_limit_conta, \
        case1_mapping_highest_ratio, virus_data, standardize_t1_threshold, standardize_t2_threshold, \
        t3_threshold, nb_read_limit_conta, mapping_highest_ratio, standardisation)

def manual_test(col_name, col_name_control, out_dir, standardisation):
    """"""

    filename_list = os.listdir(out_dir)

    count = 0
    for filename in filename_list:
        list_name_part = filename.split("_")
        bool_switch = False
        common_name_file = ""
        if filename.startswith("input") or filename.startswith("Input"):
            for i in range(1,len(list_name_part)):
                if bool_switch:
                    common_name_file = common_name_file + list_name_part[i]
                if list_name_part[i] == "file":
                    bool_switch = True
        #common_name_file = common_name_file + ".csv"  
        for filename_c in filename_list:
            list_name_part_control = filename_c.split("_")
            bool_switch_control = False
            common_name_file_control = ""
            if filename_c.startswith("control") or filename_c.startswith("Control"):
                for i in range(1,len(list_name_part_control)):
                    if bool_switch_control:
                        common_name_file_control = common_name_file_control + list_name_part_control[i] 
                    if list_name_part_control[i] == "file":
                        bool_switch_control = True
                    if common_name_file_control == common_name_file and common_name_file != "":
                        filename_control = filename_c
                        file_name_data = filename
                        count +=1
                        print(count)
                        print(file_name_data)
                        print(filename_control)
                        print("------")
                        run_analysis(out_dir, file_name_data, filename_control, col_name, col_name_control, standardisation)

def read_arg():
    """ check user arguments"""

    help_text = """
    Cont-ID is a tool to check for cross-contamination 
    using metric from previous analysis like reads number, and deduplication.
    """

    help_epilog = """
    Write description here

    See also https://github.com/
    """

    parser = argparse.ArgumentParser(description=help_text,
                                    epilog=help_epilog,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-dr', '--data_repository', type=str, required=True,
                    help='repository where the data are')
    parser.add_argument('-all', '--all_files', type=bool, default=False,
                    help='If True, will go through every file \
                    of data_repository (default False')   
    parser.add_argument('-s', '--standardisation', type=int, default=5000000,
                    help='expected number of reads for each sample to use for standardisation')                 
    parser.add_argument('-tp', '--threshold_personalisation', type=list, default=["2:1000:1.5:5:1", "0.002:500:1.5:5:1"],
                    help='personalisation threshold like ["2:1000:1.5:5:1", "0.002:500:1.5:5:1"] \
                        with ":" as separator')               
    parser.add_argument('-fd', '--file_data', type=str, default="",
                    help='file name containing data for processing')
    parser.add_argument('-fc', '--file_control', type=str, default="",
                    help='file name containing alien control data')     
    args = parser.parse_args()

    out_dir = args.data_repository
    # out_dir = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/ALL_final_version_vote/human/"
    
    standardisation = args.standardisation
    # standardisation = 5000000
    
    bool_all_file = args.all_files
    # bool_all_file = False
    
    file_name_data = args.file_data
    filename_control = args.file_control

    # threshold_case = ["2:1000:1.5"]
    # T1 will be divide by 2 or 0.002 
    # T2 divide by 1000 or 500
    # T3 divide by 1.5
    # ["2:1000:1.5", "0.002:500:1.5"]s  
    global threshold_case 
    threshold_case = args.threshold_personalisation
    #threshold_case = ["2:1000:1.5:5:1", "0.002:500:1.5:5:1"]
    #TODO make integrity check

    return out_dir, standardisation, file_name_data, filename_control, bool_all_file


if __name__ == "__main__":

    out_dir, standardisation, file_name_data, filename_control, bool_all_file = read_arg()

    col_name = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication", "Total_Reads_Nr"]
    col_name_control = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication", "Total_Reads_Nr", "Indexing"]

    if bool_all_file: 
        manual_test(col_name, col_name_control, out_dir, standardisation)
    else:
        run_analysis(out_dir, file_name_data, filename_control, col_name, col_name_control, standardisation)




  