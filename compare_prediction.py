#!/usr/bin/env python

import os
import pandas as pd

def write_result(out_dir_res_towrite, filename, batch_data, dic_result_stat, threshold_used_data):
    """
    """

    try:
        os.mkdir(out_dir_res_towrite)
    except FileExistsError:
        pass

    file_name = "Analysis_"+ filename

    batch_data.to_csv(os.path.join(out_dir_res_towrite, file_name), index=False, sep=';', encoding='utf-8')

    with open(os.path.join(out_dir_res_towrite, file_name), 'a') as out_file:
        out_file.write("\n") 

        for el in threshold_used_data.itertuples():
            out_file.write(str(el.Virus_detected) + ";" + str(el.Sample_name) + ";" + str(el.Sample_ID) + "\n")
        out_file.write("\n") 
        list_metric = ["step1_correct", "step1_uncertain", "step1_wrong", "step2_correct", "step2_uncertain", \
            "step2_wrong", "step3_correct", "step3_wrong", "all_correct", "all_wrong"]
        for el in list_metric:
            out_file.write(str(el) + ";" + str(dic_result_stat[el][0]) + ";" + str(dic_result_stat[el][1]) + "\n")

        # for key, value in dic_result_stat.items():
        #     out_file.write(str(key) + ";" + str(value[0]) + ";" + str(value[1]) + "\n")



def make_comparison(out_dir_res, filename, col_name, index_data, out_dir_res_towrite):
    """
    """

    batch_data = pd.read_csv(os.path.join(out_dir_res, filename), sep=";", header=0, \
            names=col_name, index_col=False, engine='python', skipfooter=6)

    batch_data.sort_values(by=['Virus_detected', 'Sample_name'], ignore_index=True, inplace=True)

    list_class_prediction = []
    list_index_status = []
    count_correct=0
    count_wrong=0
    count_step1_uncertain = 0
    count_step1_correct=0
    count_step2_uncertain = 0
    count_step2_correct=0
    total_class_nb = len(batch_data)
    for element in index_data.itertuples():
        for el in batch_data.itertuples():
            if element.Virus_detected == el.Virus_detected:
                if element.Sample_name == el.Sample_name:
                    step1_status = el.classification_step1
                    step2_status = el.classification_step2
                    step3_status = el.classification_step3
                    index_status = element.Indexing
                    list_index_status.append(index_status)
                    if index_status == step3_status:
                        list_class_prediction.append("True")
                        count_correct+=1
                    else:
                        list_class_prediction.append("False")
                        count_wrong+=1
                    if step1_status == "uncertain":
                        count_step1_uncertain+=1
                    else:
                        if step1_status == index_status:
                            count_step1_correct+=1
                    if step2_status == "uncertain":
                        count_step2_uncertain+=1
                    else:
                        if step2_status == index_status:
                            count_step2_correct+=1
    try:
        batch_data["result_analysis"] = list_class_prediction
    except ValueError:
        print(filename)

    batch_data["index_status"] = list_index_status

    count_step2_wrong = total_class_nb- (count_step2_correct + count_step2_uncertain)
    count_step1_wrong = total_class_nb - (count_step1_correct + count_step1_uncertain)
    count_step3_correct = count_correct
    count_step3_wrong = count_wrong

    new_count_step1_correct = "+" + str(count_step1_correct)
    new_count_step1_uncertain = "+" + str(count_step1_uncertain)
    new_count_step1_wrong = "+" + str(count_step1_wrong)
    new_count_step2_correct = "+" + str(count_step2_correct - count_step1_correct)
    new_count_step2_uncertain = "-" + str(count_step1_uncertain - count_step2_uncertain)
    new_count_step2_wrong = "+" + str(count_step2_wrong - count_step1_wrong)
    new_count_step3_wrong = "+" + str(count_wrong - count_step2_wrong)
    new_count_step3_correct = "+" + str(count_correct - count_step2_correct)

    dic_result_stat = {"step1_correct":[count_step1_correct,new_count_step1_correct], 
    "step2_correct":[count_step2_correct,new_count_step2_correct],
    "step1_uncertain":[count_step1_uncertain,new_count_step1_uncertain],
    "step2_uncertain":[count_step2_uncertain,new_count_step2_uncertain],
    "step1_wrong":[count_step1_wrong,new_count_step1_wrong],
    "step2_wrong":[count_step2_wrong,new_count_step2_wrong],
    "step3_correct":[count_step3_correct,new_count_step3_correct],
    "step3_wrong":[count_step3_wrong,new_count_step3_wrong],
    "all_correct":[count_correct,""],
    "all_wrong":[count_wrong,""]}
    
    threshold_used_data = pd.read_csv(os.path.join(out_dir_res, filename), sep=";", engine='python', names=col_name, skiprows=total_class_nb+1)
    
    write_result(out_dir_res_towrite, filename, batch_data, dic_result_stat, threshold_used_data)

def run_comparison(batch_file_list, index_bsv, index_other, dataset_name, out_dir_res, out_dir_res_towrite):
    """
    """
    col_name_index = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication", "Total_Reads_Nr", "Indexing"]
    col_name = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication", \
        "Total_Reads_Nr", "standardize_reads_nb_mapped", "mapping_ratio", "standardize_mapping_ratio", \
        "classification_step1", "classification_step2", "classification_step3", "Comment"]

    filename_list = os.listdir(index_bsv)
    for filename in filename_list: 
        if dataset_name in filename:
            try:
                index_bsv_data = pd.read_csv(os.path.join(index_bsv, filename), sep=",", header=0, names=col_name_index, index_col=False)
            except pd.errors.ParserError:
                index_bsv_data = pd.read_csv(os.path.join(index_bsv, filename), sep=";", header=0, names=col_name_index, index_col=False)
    filename_list = os.listdir(index_other)
    for filename in filename_list: 
        if dataset_name in filename:
            try:
                index_other_data = pd.read_csv(os.path.join(index_other, filename), sep=",", header=0, names=col_name_index, index_col=False)
            except pd.errors.ParserError:
                index_other_data = pd.read_csv(os.path.join(index_other, filename), sep=";", header=0, names=col_name_index, index_col=False)
    index_bsv_data.sort_values(by=['Virus_detected', 'Sample_name'], ignore_index=True, inplace=True)
    index_other_data.sort_values(by=['Virus_detected', 'Sample_name'], ignore_index=True, inplace=True)
    for filename in batch_file_list:
        # case1, case2, case1_standard, case2_standard
        if 'case_1' in filename:
            if "standart" in filename:
                if "bsv" in filename:
                    make_comparison(out_dir_res, filename, col_name, index_bsv_data, out_dir_res_towrite)
                else:
                    make_comparison(out_dir_res, filename, col_name, index_other_data, out_dir_res_towrite)
            else:
                if "bsv" in filename:
                    make_comparison(out_dir_res, filename, col_name, index_bsv_data, out_dir_res_towrite)
                else:
                    make_comparison(out_dir_res, filename, col_name, index_other_data, out_dir_res_towrite)
        else:
            if "standart" in filename:
                if "bsv" in filename:
                    make_comparison(out_dir_res, filename, col_name, index_bsv_data, out_dir_res_towrite)
                else:
                    make_comparison(out_dir_res, filename, col_name, index_other_data, out_dir_res_towrite)                
            else:
                if "bsv" in filename:
                    make_comparison(out_dir_res, filename, col_name, index_bsv_data, out_dir_res_towrite)
                else:
                    make_comparison(out_dir_res, filename, col_name, index_other_data, out_dir_res_towrite)

def run_simple_comparison(batch_file_list, out_dir_index, out_dir_res, out_dir_res_towrite):
    """
    """
    col_name_index = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication", "Total_Reads_Nr", "Indexing"]
    col_name = ["Virus_detected","Sample_name","Sample_ID","Reads_nb_mapped", "deduplication", \
        "Total_Reads_Nr", "standardize_reads_nb_mapped", "mapping_ratio", "standardize_mapping_ratio", \
        "classification_step1", "classification_step2", "classification_step3", "Comment"]

    filename_list = os.listdir(out_dir_index)
    for filename in filename_list: 
        try:
            index_data = pd.read_csv(os.path.join(out_dir_index, filename), sep=",", header=0, names=col_name_index, index_col=False)
        except pd.errors.ParserError:
            index_data = pd.read_csv(os.path.join(out_dir_index, filename), sep=";", header=0, names=col_name_index, index_col=False)
    
    index_data.sort_values(by=['Virus_detected', 'Sample_name'], ignore_index=True, inplace=True)

    for filename in batch_file_list:
        # case1, case2, case1_standard, case2_standard
        if 'case_1' in filename:
            if "standart" in filename:
                make_comparison(out_dir_res, filename, col_name, index_data, out_dir_res_towrite)
            else:
                make_comparison(out_dir_res, filename, col_name, index_data, out_dir_res_towrite)
        else:
            if "standart" in filename:
                make_comparison(out_dir_res, filename, col_name, index_data, out_dir_res_towrite)
            else:
                make_comparison(out_dir_res, filename, col_name, index_data, out_dir_res_towrite)
               
if __name__ == "__main__":


    out_dir_res = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/Human_data_V2/Result2/"
    out_dir_index = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/Human_data_V2/indexing2/"
    out_dir_res_towrite = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/Human_data_V2/Result_analysis2/"
    batch_file_list = []
    filename_list = os.listdir(out_dir_res)
    for filename in filename_list:
        batch_file_list.append(filename)
    run_simple_comparison(batch_file_list, out_dir_index, out_dir_res, out_dir_res_towrite)

    #############################################################################################################

    # out_dir_res = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/all_banana_data_V3/Result/"
    # out_dir_index = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/all_banana_data_V3/indexing/"
    # out_dir_res_towrite = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/all_banana_data_V3/Result_analysis/"
    # batch1_file_list = []
    # batch2_file_list = []
    # batch3_file_list = []
    # batch4_file_list = []
    # batch5_file_list = []    
    # smallrna_file_list = []  

    # filename_list = os.listdir(out_dir_res)
    # for filename in filename_list:
    #     if "batch1" in filename:
    #         batch1_file_list.append(filename)
    #     elif "batch2" in filename:
    #         batch2_file_list.append(filename)
    #     elif "batch3" in filename:
    #         batch3_file_list.append(filename)
    #     elif "batch4" in filename:
    #         batch4_file_list.append(filename)
    #     elif "batch5" in filename:
    #         batch5_file_list.append(filename)
    #     elif "smallrna" in filename:
    #         smallrna_file_list.append(filename)
    #     else:
    #         continue 

    # index_bsv = os.path.join(out_dir_index, "BSV")
    # index_other = os.path.join(out_dir_index, "Other")
    
    # #batch 1 case1, case2, case1_standard, case2_standard
    # run_comparison(batch1_file_list, index_bsv, index_other, "batch1", out_dir_res, out_dir_res_towrite)

    # run_comparison(batch2_file_list, index_bsv, index_other, "batch2", out_dir_res, out_dir_res_towrite)

    # run_comparison(batch3_file_list, index_bsv, index_other, "batch3", out_dir_res, out_dir_res_towrite)
    
    # run_comparison(batch4_file_list, index_bsv, index_other, "batch4", out_dir_res, out_dir_res_towrite)

    # run_comparison(batch5_file_list, index_bsv, index_other, "batch5", out_dir_res, out_dir_res_towrite)

    # run_comparison(smallrna_file_list, index_bsv, index_other, "smallrna", out_dir_res, out_dir_res_towrite)

