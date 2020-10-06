#!/usr/bin/env python
"""Script for running Kodoja pipeline modules."""

import argparse
import os
import pandas as pd
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

    return virus_data, control_data

def calculate_mapping_ratio(virus_data,control_data):
    """
    Calculate Mapped reads Nr./ Hihgest Reads Nr. Per sample in the same batch (for same virus)

    """
    list_mapping_ratio = []
    list_mapping_ratio_control = []
    max_read_nb = max(control_data["Reads_nb_mapped"])
    

    for element in control_data.itertuples():
        ratio = element.Reads_nb_mapped/max_read_nb
        list_mapping_ratio_control.append(ratio)
    
    control_data["mapping_ratio"] = list_mapping_ratio_control

    serie_virus_read = virus_data.groupby("Virus_detected")['Reads_nb_mapped'].max()
    for element in virus_data.itertuples():
        max_read_nb = serie_virus_read.loc[element.Virus_detected]
        ratio = element.Reads_nb_mapped/max_read_nb
        list_mapping_ratio.append(ratio)
    virus_data["mapping_ratio"] = list_mapping_ratio

    return virus_data, control_data

def calculate_threshold(control_data):
    """
    calculate threshold value from control
    Get all samplename from control to set target virus to contamination for those sample
    """
    

    control_name = []
    deduplication_ratio = []
    for el in control_data.itertuples():
        if el.Indexing == 1 or el.Indexing == "PRESENT":
            control_name.append(el.Sample_name)
        try:
            deduplication_ratio.append(float(el.deduplication))
        except ValueError:
            pass
    ###### T1
    # T1= [nb read >5] + [(Mapped reads Nr./ Hihgest Reads Nr. Per sample in the same run) > ((avg+3*std)/2)]
    # avg is average mapping_ratio of contaminated sample from control virus 
    try:
        mean_conta = control_data.groupby("Indexing")['mapping_ratio'].mean().loc[0]
    except TypeError:
        try:
            mean_conta = control_data.groupby("Indexing")['mapping_ratio'].mean().loc["ABSENT"]
        except KeyError:
            print("No contamination found in external control, you can be confident that you don't have cross-contamination")
            exit(0)
    # std standart deviation of contaminated sample from control virus
    try:
        std_conta = control_data.groupby("Indexing")['mapping_ratio'].std().loc[0]
    except TypeError:
        try:
            std_conta = control_data.groupby("Indexing")['mapping_ratio'].std().loc["ABSENT"]
        except KeyError:
            print("No contamination found in external control, you can be confident that you don't have cross-contamination")
            exit(0)
            # threshold is (avg+3*std)/2
    #t1_threshold = (mean_conta+3*std_conta)/2
    # NOTE test other threshold (only when strain ?)
    t1_threshold = ((mean_conta+3*std_conta))*500
    ######

    ###### T2
    #[nb read > (Hihgest Reads Nr. Per sample in the same run/1000)]
    # threshold is Hihgest Reads Nr. Per sample in the same run/1000
    #t2_threshold = int(max(control_data["Reads_nb_mapped"])/1000)
    # NOTE test other threshold (only when strain ?)
    t2_threshold = int(max(control_data["Reads_nb_mapped"])/500)
    ######

    ###### T3
    # [de-duplication rate > X]
    # threshold is X
    #t3_threshold = 75.00
    # deduce X from control data
    # average control /2
    t3_threshold = (statistics.mean(deduplication_ratio))/1.5
    
    ######

    ## refinemant threshold
    # conta
    nb_read_limit_conta = 5
    # infection
    mapping_highest_ratio = 1
    
    return t1_threshold, t2_threshold, t3_threshold, nb_read_limit_conta, mapping_highest_ratio, control_name

def classify_virus(virus_data, control_data, t1_threshold, t2_threshold, t3_threshold, 
        nb_read_limit_conta, mapping_highest_ratio, control_name):
    """
    """
    def reset_bool():
        """
        return boolean to characterize if infection (default: False)
        """
        return False, False, False, False, False, False, False

    # TODO 
    # case BSV or not (treshold formula)
    # allow the user to change the formula 

    list_classification_step1 = []
    list_classification_step2 = []
    list_classification_step3 = []
    t1_bool, t2_bool, t3_bool_infection, t3_bool_contamination, is_uncertain, is_infection, is_contamination = reset_bool()
    for element in virus_data.itertuples():
        #T1
        if element.mapping_ratio >= t1_threshold:
            t1_bool=True
        #T2
        if element.Reads_nb_mapped >= t2_threshold:
            t2_bool=True
        #T3
        try:
            deduplication_ratio = float(element.deduplication)
            if deduplication_ratio <= t3_threshold:
                t3_bool_infection=True
            else:
                t3_bool_contamination=True

        except ValueError: # RF (reference) or ND (Not enough Data)
            pass
        # store status of current virus for step 1 
        if element.Sample_name in control_name or (t1_bool==False and t2_bool==False):
            is_contamination = True
            list_classification_step1.append("contamination")
        elif t1_bool and t2_bool:
            is_infection = True
            list_classification_step1.append("infection")
        else: 
            is_uncertain = True
            list_classification_step1.append("uncertain")
        # store status of current virus for step 2
        if is_infection:
            list_classification_step2.append("infection")
        if is_contamination:
            list_classification_step2.append("contamination")
        if is_uncertain:
            if t3_bool_infection:
                is_infection = True
                is_uncertain = False
                list_classification_step2.append("infection")
            if t3_bool_contamination:
                is_contamination = True
                is_uncertain = False
                list_classification_step2.append("contamination")
            if t3_bool_infection==False and t3_bool_contamination==False:
                is_uncertain = True
                list_classification_step2.append("uncertain")
        # store status for step 3
        if is_infection:
            list_classification_step3.append("infection")
        if is_contamination:
            list_classification_step3.append("contamination")
        if is_uncertain:
            #T4 (refinement)
            if element.mapping_ratio == mapping_highest_ratio:
                list_classification_step3.append("infection")
                is_uncertain = False
            elif element.Reads_nb_mapped <= nb_read_limit_conta:
                list_classification_step3.append("contamination")
                is_uncertain = False
            if is_uncertain:
                list_classification_step3.append("uncertain")
        t1_bool, t2_bool, t3_bool_infection, t3_bool_contamination, is_uncertain, is_infection, is_contamination = reset_bool()
                  
    virus_data["Confident classification (step1)"] = list_classification_step1
    virus_data["Standart classification (step2)"] = list_classification_step2
    virus_data["Total classification (step3)"] = list_classification_step3

    return virus_data

def write_result(out_dir, file_name_data, virus_data, t1_threshold, t2_threshold, t3_threshold, nb_read_limit_conta, mapping_highest_ratio):
    """
    """
    path = file_name_data.split("/")
    result_folder_name = os.path.join(out_dir, "Result_" + str(path[0]))
    file_name = "Result_" + str(path[1])
    try:
        os.mkdir(result_folder_name)
    except FileExistsError:
        pass

    virus_data.to_csv(os.path.join(result_folder_name, file_name), index=False, sep=';', encoding='utf-8')
    
    with open(os.path.join(result_folder_name, file_name), 'a') as out_file:
        out_file.write("\n")
        out_file.write("mapping ratio threshold (step1 t1);" + str(t1_threshold) + "\n")
        out_file.write("reads_nb_mapped threshold (step1 t2);" + str(t2_threshold) + "\n")
        out_file.write("deduplication threshold (step2 t3);" + str(t3_threshold) + "\n")
        out_file.write("nb_read_limit_conta (step3 t4_1);" + str(nb_read_limit_conta) + "\n")
        out_file.write("mapping_highest_ratio (step3 t4_2);" + str(mapping_highest_ratio) + "\n")

def run_analysis(out_dir, file_name_data, file_name_control, col_name, col_name_control):

    print(file_name_data)

    # open input file
    virus_data, control_data = open_file(out_dir, file_name_data, file_name_control, col_name, col_name_control)
    # add mapping ratio data
    virus_data, control_data = calculate_mapping_ratio(virus_data,control_data)
    # calculate thresshold value from control
    t1_threshold, t2_threshold, t3_threshold, nb_read_limit_conta, mapping_highest_ratio, \
        control_name = calculate_threshold(control_data)
    # make virus classification
    virus_data = classify_virus(virus_data, control_data, t1_threshold, t2_threshold, t3_threshold, 
        nb_read_limit_conta, mapping_highest_ratio, control_name)

    write_result(out_dir, file_name_data, virus_data, t1_threshold, t2_threshold, 
        t3_threshold, nb_read_limit_conta, mapping_highest_ratio)
        
    print(t1_threshold)
    print(t2_threshold)
    print(t3_threshold)
    print(nb_read_limit_conta)
    print(mapping_highest_ratio)

if __name__ == "__main__":

    #TODO argparse use argument
    out_dir = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/"
    # BSV strain by strain
    # file_name_data = "input_file_bsv_strain/Input_file_batch$_bsv_strains.csv"
    # BSV all strain
    # file_name_data = "input_file5_strain/Input_file_batch$.csv"
    # other virus cmv bbmmv bbtv ...

    # file_name_data = "other_virus/other_virus_batch$.csv"
    
    # file_name_control = "control/Input_file_control_batch$.csv"
    col_name = ["Virus_detected","Sample_name","Reads_nb_mapped", "deduplication"]
    col_name_control = ["Virus_detected","Sample_name","Reads_nb_mapped", "deduplication", "Indexing"]

    # for i in range(1,5):
    #     print( "round: " + str(i))
    #     file_name_data2 = file_name_data.replace("$", str(i))
    #     file_name_control2 = file_name_control.replace("$", str(i))
    #     run_analysis(out_dir, file_name_data2, file_name_control2, col_name, col_name_control)


    file_name_control = "control/Input_file_control_batch1.csv"
    file_name_data = "input_file5_strain/Input_file_batch1.csv"

    run_analysis(out_dir, file_name_data, file_name_control, col_name, col_name_control)