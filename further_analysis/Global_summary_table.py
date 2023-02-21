#!/usr/bin/env python
#author : Valentin Rollin 6 march 2021 
#author : Johan Rollin 27 june 2021 

import os
import pandas as pd

def make_calc_comp(data_base, data_conf_matrix):
    """
    """
    data = []
    for index,element in data_base.iterrows():
        for i,el in data_conf_matrix.iterrows():
            if element.filename == el.filename:
                unknown = element.Unconfirmed
                count_TP = el.TP
                count_TN = el.TN
                count_FP = el.FP
                count_FN = el.FN

                if count_TP + count_FN == 0:
                    dse = "unknown"
                else:
                    dse = count_TP / (count_TP + count_FN)
                
                if count_TN + count_FP == 0:
                    dsp ="unknown"
                else:                                       
                    dsp = count_TN / (count_TN + count_FP)

                denominator = count_TP + count_TN + count_FP + count_FN 
                if denominator == 0:
                    accuracy ="unknown"
                else:
                    accuracy = (count_TP + count_TN) / (count_TP + count_TN + count_FP + count_FN)

                denominator_u = count_TP + count_TN + count_FP + count_FN 
                if denominator_u == 0:
                    accuracy_u ="unknown"
                else:
                    accuracy_u = (count_TP + count_TN) / (count_TP + count_TN + count_FP + count_FN + unknown)

            ############### 
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

                data.append({"filename": el.filename, "Case": "-", "confidence": "comparison", "dse": dse, "dsp": dsp, "accuracy": accuracy,\
                    "accuracy_u": accuracy_u, "fdr": fdr, "f_omission_r": f_omission_r})                        
    df = pd.DataFrame(data)

    return df

def make_calc(data_base, data_conf_matrix):
    """
    """
    data = []
    for index,element in data_base.iterrows():
        for i,el in data_conf_matrix.iterrows():
            if element.filename == el.filename:
                if element.Case == el.Case:
                    unknown3 = element.votes3_unknown
                    unknown2 = element.votes2_unknown
                    unknown = element.overall_unknown
                    count_TP3 = el.votes3_TP
                    count_TP2 = el.votes2_TP
                    count_TP = el.overall_TP
                    count_TN3 = el.votes3_TN
                    count_TN2 = el.votes2_TN
                    count_TN = el.overall_TN
                    count_FP3 = el.votes3_FP
                    count_FP2 = el.votes2_FP
                    count_FP = el.overall_FP
                    count_FN3 = el.votes3_FN
                    count_FN2 = el.votes2_FN
                    count_FN = el.overall_FN

                    # unknown = int(unknown)
                    # count_TP = int(count_TP)
                    # count_TN = int(count_TN)
                    # count_FP = int(count_FP)
                    # count_FN = int(count_FN)

                    if count_TP + count_FN == 0:
                        dse = "unknown"
                    else:
                        dse = count_TP / (count_TP + count_FN)
                    if count_TP3 + count_FN3 == 0:
                        dse3 = "unknown"
                    else:
                        dse3 = count_TP3 / (count_TP3 + count_FN3)
                    if count_TP2 + count_FN2 == 0:
                        dse2 = "unknown"
                    else:
                        dse2 = count_TP2 / (count_TP2 + count_FN2)


                    if count_TN + count_FP == 0:
                        dsp ="unknown"
                    else:
                        dsp = count_TN / (count_TN + count_FP)
                    if count_TN3 + count_FP3 == 0:
                        dsp3 ="unknown"
                    else:
                        dsp3 = count_TN3 / (count_TN3 + count_FP3)
                    if count_TN2 + count_FP2 == 0:
                        dsp2 ="unknown"
                    else:
                        dsp2 = count_TN2 / (count_TN2 + count_FP2)                                                

                    denominator = count_TP + count_TN + count_FP + count_FN 
                    if denominator == 0:
                        accuracy ="unknown"
                    else:
                        accuracy = (count_TP + count_TN) / (count_TP + count_TN + count_FP + count_FN)
                    denominator3 = count_TP3 + count_TN3 + count_FP3 + count_FN3 
                    if denominator3 == 0:
                        accuracy3 ="unknown"
                    else:
                        accuracy3 = (count_TP3 + count_TN3) / (count_TP3 + count_TN3 + count_FP3 + count_FN3)
                    denominator2 = count_TP2 + count_TN2 + count_FP2 + count_FN2 
                    if denominator2 == 0:
                        accuracy2 ="unknown"
                    else:
                        accuracy2 = (count_TP2 + count_TN2) / (count_TP2 + count_TN2 + count_FP2 + count_FN2)                                                

                    denominator_u = count_TP + count_TN + count_FP + count_FN 
                    if denominator_u == 0:
                        accuracy_u ="unknown"
                    else:
                        accuracy_u = (count_TP + count_TN) / (count_TP + count_TN + count_FP + count_FN + unknown)
                    denominator_u3 = count_TP3 + count_TN3 + count_FP3 + count_FN3
                    if denominator_u3 == 0:
                        accuracy_u3 ="unknown"
                    else:
                        accuracy_u3 = (count_TP3 + count_TN3) / (count_TP3 + count_TN3 + count_FP3 + count_FN3 + unknown3)
                    denominator_u2 = count_TP2 + count_TN2 + count_FP2 + count_FN2 
                    if denominator_u2 == 0:
                        accuracy_u2 ="unknown"
                    else:
                        accuracy_u2 = (count_TP2 + count_TN2) / (count_TP2 + count_TN2 + count_FP2 + count_FN2 + unknown2)

                    ############### 
                    # FDR FP/(FP+TP)
                    # FOR Fn/(FN+TN)

                    if count_TP + count_FP == 0:
                        fdr = "unknown"
                    else:
                        fdr = count_FP / (count_TP + count_FP)
                    if count_TP3 + count_FP3 == 0:
                        fdr3 = "unknown"
                    else:
                        fdr3 = count_FP3 / (count_TP3 + count_FP3)
                    if count_TP2 + count_FP2 == 0:
                        fdr2 = "unknown"
                    else:
                        fdr2 = count_FP2 / (count_TP2 + count_FP2)                        

                    if count_TN + count_FN == 0:
                        f_omission_r = "unknown"
                    else:
                        f_omission_r = count_FN / (count_TN + count_FN)
                    if count_TN3 + count_FN3 == 0:
                        f_omission_r3 = "unknown"
                    else:
                        f_omission_r3 = count_FN3 / (count_TN3 + count_FN3)
                    if count_TN2 + count_FN2 == 0:
                        f_omission_r2 = "unknown"
                    else:
                        f_omission_r2 = count_FN2 / (count_TN2 + count_FN2)                        

                    data.append({"filename": el.filename, "Case": el.Case, "confidence": "3_votes", "dse": dse3, "dsp": dsp3, "accuracy": accuracy3,\
                        "accuracy_u": accuracy_u3, "fdr": fdr3, "f_omission_r": f_omission_r3})
                    data.append({"filename": el.filename, "Case": el.Case, "confidence": "2_votes", "dse": dse2, "dsp": dsp2, "accuracy": accuracy2,\
                        "accuracy_u": accuracy_u2, "fdr": fdr2, "f_omission_r": f_omission_r2})
                    data.append({"filename": el.filename, "Case": el.Case, "confidence": "overall", "dse": dse, "dsp": dsp, "accuracy": accuracy,\
                        "accuracy_u": accuracy_u, "fdr": fdr, "f_omission_r": f_omission_r})                        
    df = pd.DataFrame(data)

    return df


def run_analysis(data_base, data_base_comp, data_conf_matrix, data_conf_matrix_comp, out_dir_res_towrite):
    """
    """
 
    #calc comparison
    diag_df_comp = make_calc_comp(data_base_comp, data_conf_matrix_comp)

    diag_df= make_calc(data_base, data_conf_matrix)

    try:
        os.mkdir(out_dir_res_towrite)
    except FileExistsError:
        pass    
    diagnostic_file = out_dir_res_towrite + "diagnostic.csv"
    diag_df.to_csv(diagnostic_file, sep=";", index=False)
    diag_df_comp.to_csv(diagnostic_file, sep=";", index=False, mode='a', header=False)                                   
               
if __name__ == "__main__":

    out_dir = "/mnt/c/Users/johan/OneDrive/Bureau/bioinfo/Wei_virus_test/Key_sample/ALL_final_version_vote/ContID_human_relax/Global/"
    out_dir_res_towrite = out_dir + "Analysis/"
    file1 = out_dir + "project_analysis_base.csv"
    file2 = out_dir + "project_analysis_base_comp.csv"
    file3 = out_dir + "project_analysis_conf_matrix.csv"
    file4 = out_dir + "project_analysis_conf_matrix_comp.csv"

    # header_table1 = "filename;Case;3_votes_correct;3_votes_unknown;3_votes_wrong;2_votes_correct;2_votes_unknown;2_votes_wrong;overall_correct;overall_unknown;overall_wrong \n"
    # header_table2 = "filename;Correct;Unconfirmed;Wrong \n"
    # header_table3 = "filename;Case;3_votes_TP;3_votes_TN;3_votes_FP;3_votes_FN;2_votes_TP;2_votes_TN;2_votes_FP;2_votes_FN;overall_TP;overall_TN;overall_FP;overall_FN \n"
    # header_table4 = "filename;TP;TN;FP;FN \n"
    col_name_index1 = ["filename","Case","votes3_correct","votes3_unknown", "votes3_wrong", "votes2_correct", "votes2_unknown", "votes2_wrong", "overall_correct", "overall_unknown", "overall_wrong"]
    col_name_index3 = ["filename","Case","votes3_TP","votes3_TN", "votes3_FP", "votes3_FN", "votes2_TP", "votes2_TN", "votes2_FP", "votes2_FN", "overall_TP", "overall_TN", "overall_FP", "overall_FN"]
    col_name_index4 = ["filename","TP","TN","FP", "FN"]

    data_base = pd.read_csv(file1, sep = ";", header = 0, names = col_name_index1, index_col = False)
    data_base_comp = pd.read_csv(file2, sep = ";", header = 0, index_col = False)
    data_conf_matrix = pd.read_csv(file3, sep = ";", header = 0, names = col_name_index3, index_col = False)
    data_conf_matrix_comp = pd.read_csv(file4, sep = ";", header = 0, names = col_name_index4, index_col = False)
    
    run_analysis(data_base, data_base_comp, data_conf_matrix, data_conf_matrix_comp, out_dir_res_towrite)
    print("job done for: " + out_dir)
    



