#!/usr/bin/env python
#author : Valentin Rollin 6 march 2021 
#author : Johan Rollin 6 march 2021 

import os
import pandas as pd


def write_output1(out_file, list1_name, list1_c1, list2_c1, list1_name2, list1_c2, list2_c2, list1_bc_name, list1_bc, list2_bc):
    """
    """
    size1 = len(list1_name)
    size2 = len(list1_bc_name)
    for i in range(0,size1):
        if i < size2:
            out_file.write(str(list1_name[i]) + ";" + str(list1_c1[i]) + ";" + str(list2_c1[i]) + ";" + str(list1_name2[i]) \
                + ";" + str(list1_c2[i]) + ";" + str(list2_c2[i]) + ";" + str(list1_bc_name[i]) + ";" + str(list1_bc[i])  
                + ";" + str(list2_bc[i]) + "\n")
        else: 
            out_file.write(str(list1_name[i]) + ";" + str(list1_c1[i]) + ";" + str(list2_c1[i]) + ";" + str(list1_name2[i]) \
                + ";" + str(list1_c2[i]) + ";" + str(list2_c2[i]) + "\n")

def write_output2(out_file, list_name, list_c1, list_name2, list_c2, list_bc_name, list_bc):
    """
    """
    size1 = len(list_name)
    size2 = len(list_bc_name)
    for i in range(0,size1):
        if i < size2:
            out_file.write(str(list_name[i]) + ";" + str(list_c1[i]) + ";" + str(list_name2[i]) + ";" + str(list_c2[i]) 
                + ";" + str(list_bc_name[i]) + ";" + str(list_bc[i]) + "\n")    
        else: 
            out_file.write(str(list_name[i]) + ";" + str(list_c1[i]) + ";" + str(list_name2[i]) + ";" + str(list_c2[i]) + "\n")


def make_addcount(count_current_correct, count_current_uncertain, count_current_wrong, count_previous_correct, count_previous_uncertain, count_previous_wrong):
    """
    """
    new_count_correct = "+" + str(count_current_correct - count_previous_correct)
    new_count_uncertain = "-" + str(count_previous_uncertain - count_current_uncertain)
    new_count_wrong = "+" + str(count_current_wrong - count_previous_wrong)

    return new_count_correct, new_count_uncertain, new_count_wrong

def make_calc(count_FN, count_FP, count_TN, count_TP):
    """
    """
    if count_TP + count_FN == 0:
        dse = "unknown"

    else:
        dse = count_TP / (count_TP + count_FN)
    if count_TN + count_FP == 0:
        dsp ="unknown"
    else:
        dsp = count_TN / (count_TN + count_FP)

    accuracy = (count_TP + count_TN) / (count_TP + count_TN + count_FP + count_FN)

    ############### TODO change to make 
    # FDR FP/(FP+TP)
    # FOR Fn/(FN+TN)

    if count_TP + count_FP == 0:
        fdr = "unknown"
    else:
        fdr = count_FP / (count_TP + count_FP)
    if count_TN + count_FN == 0:
        f_omission_r = "unknown"
    else:
        f_omission_r = count_FN / (count_TN + count_FN)




    return dse, dsp, accuracy, fdr, f_omission_r


def run_simple_comparison(index_data, batch_data, col_name):
    """
    """
    list_index_status = []
    list_class_prediction = []
    count_bc_TP = 0
    count_bc_TN = 0
    count_bc_FP = 0
    count_bc_FN = 0
    count_bc_unknown = 0
    count_bc_correct = 0
    count_bc_wrong = 0
    for index,element in index_data.iterrows():
        for i,el in batch_data.iterrows():
            if element.Virus_detected == el.Virus_detected:
                if element.Sample_ID == el.Sample_ID:
                #if element.Sample_name == el.Sample_name:
                    check_col = el[col_name]
                    index_status = element.Indexing
                    list_index_status.append(index_status)
                    if check_col == element.Indexing:
                        count_bc_correct += 1 
                        list_class_prediction.append("TRUE")
                        if check_col == "infection":
                            count_bc_TP += 1
                        elif check_col == "contamination": 
                            count_bc_TN += 1
                    else:
                        list_class_prediction.append("FALSE")
                        if check_col == "unconfirmed":
                            count_bc_unknown += 1
                        elif check_col == "uncertain":
                            count_bc_unknown += 1
                        elif check_col == "infection":
                            count_bc_wrong += 1
                            count_bc_FP += 1
                        elif check_col == "contamination":
                            count_bc_wrong += 1
                            count_bc_FN += 1

    return count_bc_FN, count_bc_FP, count_bc_TN,count_bc_TP, count_bc_correct, count_bc_unknown, count_bc_wrong, list_class_prediction, list_index_status

def run_analysis(out_dir, out_dir_res, out_dir_index, out_dir_res_towrite, col_name_index, col_name, filename_result, filename_index):
    """
    """
    try:
        index_data = pd.read_csv(os.path.join(out_dir_index, filename_index), sep = ",", header = 0, names = col_name_index, index_col = False)
    except pd.errors.ParserError:
        index_data = pd.read_csv(os.path.join(out_dir_index, filename_index), sep = ";", header = 0, names = col_name_index, index_col = False)
        
    index_data.sort_values(by = ["Virus_detected", "Sample_name"], ignore_index = True, inplace = True)

    batch_data = pd.read_csv(os.path.join(out_dir_res, filename_result), sep = ";", header = 0, \
                names = col_name, index_col = False)

    batch_data.sort_values(by = ["Virus_detected", "Sample_name"], ignore_index = True, inplace = True)

    #Comparison both case
    count_bc_FN, count_bc_FP, count_bc_TN, count_bc_TP, count_bc_correct, count_bc_unknown, count_bc_wrong, \
    list_class_prediction, list_index_status = run_simple_comparison(index_data, batch_data, "comparison")

    #step 1 case 1
    count_s1_c1_FN, count_s1_c1_FP, count_s1_c1_TN, count_s1_c1_TP, count_s1_c1_correct, count_s1_c1_unknown, count_s1_c1_wrong, \
    list_class_prediction_tmp, list_index_status_tmp = run_simple_comparison(index_data, batch_data, "classification_step1_case1")

    #step 2 case 1
    count_s2_c1_FN, count_s2_c1_FP, count_s2_c1_TN, count_s2_c1_TP, count_s2_c1_correct, count_s2_c1_unknown, count_s2_c1_wrong, \
    list_class_prediction_tmp, list_index_status_tmp = run_simple_comparison(index_data, batch_data, "classification_step2_case1")
    
    #step 3 case 1
    count_s3_c1_FN, count_s3_c1_FP, count_s3_c1_TN,count_s3_c1_TP, count_s3_c1_correct, count_s3_c1_unknown, count_s3_c1_wrong, \
    list_class_prediction_tmp, list_index_status_tmp = run_simple_comparison(index_data, batch_data, "classification_step3_case1")

    #step 1 case 2
    count_s1_c2_FN, count_s1_c2_FP, count_s1_c2_TN, count_s1_c2_TP, count_s1_c2_correct, count_s1_c2_unknown, count_s1_c2_wrong, \
    list_class_prediction_tmp, list_index_status_tmp = run_simple_comparison(index_data, batch_data, "classification_step1_case2")

    #step 2 case 2
    count_s2_c2_FN, count_s2_c2_FP, count_s2_c2_TN, count_s2_c2_TP, count_s2_c2_correct, count_s2_c2_unknown, count_s2_c2_wrong, \
    list_class_prediction_tmp, list_index_status_tmp = run_simple_comparison(index_data, batch_data, "classification_step2_case2")

    #step 3 case 2
    count_s3_c2_FN, count_s3_c2_FP, count_s3_c2_TN, count_s3_c2_TP, count_s3_c2_correct, count_s3_c2_unknown, count_s3_c2_wrong, \
    list_class_prediction_tmp, list_index_status_tmp = run_simple_comparison(index_data, batch_data, "classification_step3_case2")

    add_s2_c1_correct, add_s2_c1_unknown, add_s2_c1_wrong = make_addcount(count_s2_c1_correct, count_s2_c1_unknown, count_s2_c1_wrong, count_s1_c1_correct, count_s1_c1_unknown, count_s1_c1_wrong)
    add_s3_c1_correct, add_s3_c1_unknown, add_s3_c1_wrong = make_addcount(count_s3_c1_correct, count_s3_c1_unknown, count_s3_c1_wrong, count_s2_c1_correct, count_s2_c1_unknown, count_s2_c1_wrong)
    add_s2_c2_correct, add_s2_c2_unknown, add_s2_c2_wrong = make_addcount(count_s2_c2_correct, count_s2_c2_unknown, count_s2_c2_wrong, count_s1_c2_correct, count_s1_c2_unknown, count_s1_c2_wrong)
    add_s3_c2_correct, add_s3_c2_unknown, add_s3_c2_wrong = make_addcount(count_s3_c2_correct, count_s3_c2_unknown, count_s3_c2_wrong, count_s2_c2_correct, count_s2_c2_unknown, count_s2_c2_wrong)
    
    list1_name = ["Step 1 correct", "Step 1 unknown", "Step 1 wrong", "Step 2 correct", "Step 2 unknown", "Step 2 wrong", "Step 3 correct", "Step 3 unknown", "Step 3 wrong"]
    list1_bc_name = ["Correct", "Unconfirmed", "Wrong"]
    list1_c1 = [count_s1_c1_correct, count_s1_c1_unknown, count_s1_c1_wrong, count_s2_c1_correct, count_s2_c1_unknown, count_s2_c1_wrong, count_s3_c1_correct, count_s3_c1_unknown, count_s3_c1_wrong]
    list1_c2 = [count_s1_c2_correct, count_s1_c2_unknown, count_s1_c2_wrong, count_s2_c2_correct, count_s2_c2_unknown, count_s2_c2_wrong, count_s3_c2_correct, count_s3_c2_unknown, count_s3_c2_wrong]
    list1_bc = [count_bc_correct, count_bc_unknown, count_bc_wrong]   
    list2_c1 = [count_s1_c1_correct, count_s1_c1_unknown, count_s1_c1_wrong, add_s2_c1_correct, add_s2_c1_unknown, add_s2_c1_wrong, add_s3_c1_correct, add_s3_c1_unknown, add_s3_c1_wrong]
    list2_c2 = [count_s1_c2_correct, count_s1_c2_unknown, count_s1_c2_wrong, add_s2_c2_correct, add_s2_c2_unknown, add_s2_c2_wrong, add_s3_c2_correct, add_s3_c2_unknown, add_s3_c2_wrong]
    list2_bc = [count_bc_correct, count_bc_unknown, count_bc_wrong]
    
    list3_name = ["Step 1 TP", "Step 1 TN", "Step 1 FP", "Step 1 FN", "Step 2 TP", "Step 2 TN", "Step 2 FP", "Step 2 FN", "Step 3 TP", "Step 3 TN", "Step 3 FP", "Step 3 FN"]
    list3_bc_name = ["TP", "TN", "FP", "FN"]
    list3_c1 = [count_s1_c1_TP, count_s1_c1_TN, count_s1_c1_FP, count_s1_c1_FN, count_s2_c1_TP, count_s2_c1_TN, count_s2_c1_FP, count_s2_c1_FN, count_s3_c1_TP, count_s3_c1_TN, count_s3_c1_FP, count_s3_c1_FN]
    list3_c2 = [count_s1_c2_TP, count_s1_c2_TN, count_s1_c2_FP, count_s1_c2_FN, count_s2_c2_TP, count_s2_c2_TN, count_s2_c2_FP, count_s2_c2_FN, count_s3_c2_TP, count_s3_c2_TN, count_s3_c2_FP, count_s3_c2_FN]
    list3_bc = [count_bc_TP, count_bc_TN, count_bc_FP, count_bc_FN]

    #calc comparison
    dse_bc, dsp_bc, accuracy_bc, fdr_bc, for_bc = make_calc(count_bc_FN, count_bc_FP, count_bc_TN,count_bc_TP)

    #calc step 1 case 1
    dse_s1_c1, dsp_s1_c1, accuracy_s1_c1, fdr_s1_c1, for_s1_c1 = make_calc(count_s1_c1_FN, count_s1_c1_FP, count_s1_c1_TN,count_s1_c1_TP)

    #calc step 2 case 1
    dse_s2_c1, dsp_s2_c1, accuracy_s2_c1, fdr_s2_c1, for_s2_c1 = make_calc(count_s2_c1_FN, count_s2_c1_FP, count_s2_c1_TN,count_s2_c1_TP)

    #calc step 3 case 1
    dse_s3_c1, dsp_s3_c1, accuracy_s3_c1, fdr_s3_c1, for_s3_c1 = make_calc(count_s3_c1_FN, count_s3_c1_FP, count_s3_c1_TN,count_s3_c1_TP)

    #calc step 1 case 2
    dse_s1_c2, dsp_s1_c2, accuracy_s1_c2, fdr_s1_c2, for_s1_c2 = make_calc(count_s1_c2_FN, count_s1_c2_FP, count_s1_c2_TN,count_s1_c2_TP)

    #calc step 2 case 2 
    dse_s2_c2, dsp_s2_c2, accuracy_s2_c2, fdr_s2_c2, for_s2_c2 = make_calc(count_s2_c2_FN, count_s2_c2_FP, count_s2_c2_TN,count_s2_c2_TP)

    #calc step 3 case 2
    dse_s3_c2, dsp_s3_c2, accuracy_s3_c2, fdr_s3_c2, for_s3_c2 = make_calc(count_s3_c2_FN, count_s3_c2_FP, count_s3_c2_TN,count_s3_c2_TP)
    list4_name = ["Step 1 DSE", "Step 1 DSP", "Step 1 accuracy", "Step 1 FDR", "Step 1 FOR", "Step 2 DSE", "Step 2 DSP", "Step 2 accuracy", "Step 2 FDR", "Step 2 FOR", "Step 3 DSE", "Step 3 DSP", "Step 3 accuracy", "Step 3 FDR", "Step 3 FOR"]
    list4_bc_name = ["DSE", "DSP", "accuracy", "FDR", "FOR"]
    list4_c1 = [dse_s1_c1, dsp_s1_c1, accuracy_s1_c1, fdr_s1_c1, for_s1_c1, dse_s2_c1, dsp_s2_c1, accuracy_s2_c1, fdr_s2_c1, for_s2_c1, dse_s3_c1, dsp_s3_c1, accuracy_s3_c1, fdr_s3_c1, for_s3_c1]
    list4_c2 = [dse_s1_c2, dsp_s1_c2, accuracy_s1_c2, fdr_s1_c2, for_s1_c2, dse_s2_c2, dsp_s2_c2, accuracy_s2_c2, fdr_s2_c2, for_s2_c2, dse_s3_c2, dsp_s3_c2, accuracy_s3_c2, fdr_s3_c2, for_s3_c2]
    list4_bc = [dse_bc, dsp_bc, accuracy_bc, fdr_bc, for_bc]

    try:
        os.mkdir(out_dir_res_towrite)
    except FileExistsError:
        pass

    file_name = "Analysis_"+ filename_result
    batch_data["Prediction analysis"] = list_class_prediction
    batch_data["Index status"] = list_index_status
    batch_data.to_csv(os.path.join(out_dir_res_towrite, file_name), index=False, sep=";", encoding="utf-8")
    try:
        os.mkdir(out_dir_res_towrite)
    except FileExistsError:
        pass
    with open(os.path.join(out_dir_res_towrite, file_name), "a") as out_file:
        out_file.write("\n") 
        out_file.write("Case 1 ;;;Case2 ;;;Comparison\n")
        write_output1(out_file, list1_name, list1_c1, list2_c1, list1_name, list1_c2, list2_c2, list1_bc_name, list1_bc, list2_bc)
        out_file.write("\n")
        out_file.write("Case 1 ;;Case2 ;;Comparison\n")
        write_output2(out_file, list3_name, list3_c1, list3_name, list3_c2, list3_bc_name, list3_bc)
        out_file.write("\n")
        out_file.write("Case 1 ;;Case2 ;;Comparison\n")
        write_output2(out_file, list4_name, list4_c1, list4_name, list4_c2, list4_bc_name, list4_bc)
        out_file.write("\n")
        
                                       
               
if __name__ == "__main__":
    
    out_dir = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/ALL_final_version/banana/"
    out_dir_res = out_dir + "Result/"
    out_dir_index = out_dir + "indexing/"
    out_dir_res_towrite = out_dir + "Analysis/"
    
    col_name_index = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication", "Total_Reads_Nr", "Indexing"]

    col_name = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication",
    "Total_Reads_Nr", "standardize_reads_nb_mapped", "mapping_ratio", "standardize_mapping_ratio", 
    "classification_step1_case1", "classification_step2_case1", "classification_step3_case1", "Comment_case1",
    "classification_step1_case2", "classification_step2_case2", "classification_step3_case2", "Comment_case2",
    "comparison"] 

    # filename_index = "Indexing_file_batch5_dsRNA_Musa.csv"
    # filename_result = "Result_Input_file_batch5_dsRNA_Musa.csv"
    filename_list = os.listdir(out_dir_res)
    filename_list_index = os.listdir(out_dir_index)
    count = 0
    for filename in filename_list:
        
        list_name_part = filename.split("_")
        bool_switch = False
        common_name_file = ""
        if filename.startswith("Result"):
            for i in range(0,len(list_name_part)):
                if bool_switch:
                    common_name_file = common_name_file + list_name_part[i]
                if list_name_part[i] == "file":
                    bool_switch = True
        for filename_i in filename_list_index:
            list_name_part_index = filename_i.split("_")
            bool_switch_index = False
            common_name_file_index = ""
            for i in range(0,len(list_name_part_index)):
                if bool_switch_index:
                    common_name_file_index = common_name_file_index + list_name_part_index[i]
                if list_name_part_index[i] == "file":
                    bool_switch_index = True
                if common_name_file_index == common_name_file and common_name_file != "":
                    filename_index = filename_i
                    filename_result = filename
                    count +=1
                    print(count)
                    print(filename_result)
                    print(filename_index)
                    print("------")
                    run_analysis(out_dir, out_dir_res, out_dir_index, out_dir_res_towrite, col_name_index, col_name, filename_result, filename_index)                    
    



