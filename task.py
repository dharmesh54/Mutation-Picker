import sys
import os
import pandas as pd
import numpy
import numpy as np
from pandas.api.types import is_string_dtype
import warnings
import csv
import time

def run_script(file_name, path, Wild_id, total_gen):
        
        #File_name = self.file_input.text()
        #if "/" in File_name:
             #file_name = File_name.split("/")[-1]
        #else:
            #file_name =  File_name
        #path = self.path_input.text()
        #Wild_id = self.wild_id_input.text()
        #total_gen = int(self.gen_input.text())
        #time = self.time_input.text() # new input field for time
        #load = self.load_input.text() # new input field for load
        #print(f'File Name: {File_name}')
        #print(f'Path: {path}')
        #print(f'Wild ID: {wild_id}')
        #file_name = sys.argv[1]
        #path = sys.argv[2]
        #Wild_id = sys.argv[3]
        #total_gen = int(sys.argv[4])
        files=os.listdir(path)
        def dictionery(file, path, filename):
           "reading text file to dictionery"
           protein = {}
           for i in files:   ## files in the working directory
             if i  == filename:
                ip = os.path.join(path,i)
                with open(ip) as f:
                  for line in f:
                        key = line.strip() #
                        value = f.readline().strip()
                        protein[key] = list(value) ##apending key:value pair to the directory
                  return protein    #converting protein text file into the dictionery 
        
        dict_protein = dictionery(files, path, file_name)
        length = len(dict_protein)
        df=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in dict_protein.items()]))
        df3 =  df.copy(deep=True)
        
        
        def length_protein(df_s, identifier):
                list_output = []
                p_name = df_s.filter(regex= identifier,axis=1)
                if len(p_name.columns) == 1:
                  colname = p_name.columns[0]
                  for col in df_s.columns:
                    if col == colname:  
                        pass
                    else:
                        parent_column = [x for x in df_s[colname] if str(x) != 'nan']
                        query_column = [x for x in df_s[col] if str(x) != 'nan']
                    if len(parent_column) != len(query_column):
                        if colname not in list_output: 
                           list_output.append(colname)
                        list_output.append(col)
                        list_output.append(len(query_column))
                return list_output
        
        def truncation(df_dummy, identifier):
                 list_truncated_mutant = []
                 p_name = df_dummy.filter(regex= identifier,axis=1)
                 if len(p_name.columns) == 1: 
                  colname = p_name.columns[0]
                  for col in df_dummy.columns:
                    if col == colname:  
                        pass
                    else:
                        #list_output = []
                        parent_column = [x for x in df_dummy[colname] if str(x) != 'nan']
                        query_column = [x for x in df_dummy[col] if str(x) != 'nan']
                        M_position = [i for i,val in enumerate(parent_column) if val =="M"]     
                        for j in range(1, len(M_position)): 
                            list_of_match_amino_acids = [] 
                            for i in range(len(query_column)):
                                l = M_position[j]+i
                                if len(parent_column) < l:
                                    list_of_match_amino_acids.append("False")
                                if len(parent_column) > l:
                                         if(parent_column[l]==query_column[i]): 
                                                list_of_match_amino_acids.append("True")
                                         else:
                                             list_of_match_amino_acids.append("False")
                                           
                            fal_eval = list_of_match_amino_acids.count("False")
                            percent_identity = fal_eval/len(list_of_match_amino_acids)
                            if percent_identity <= 0.05:
                               if colname not in list_truncated_mutant:
                                   list_truncated_mutant.append(colname)
                               if col != colname:
                                   list_truncated_mutant.append(col)
                               list_truncated_mutant.append(M_position[j])
                               break
                 return list_truncated_mutant
                 
                 
        def mutation(df_dummy, identifier):
                    df_dummy.loc[:,'position'] = pd.Series(list(range(1,length +1)))
                    p_name = df_dummy.filter(regex= identifier,axis=1)
                    if len(p_name.columns) == 1:
                      colname = p_name.columns[0] 
                      for col in df_dummy.columns:
                        if col == colname or col == "position":  
                            pass
                        else:
                            new_col = "new" + col
                            list_output = []
                            parent_column = [x for x in df_dummy[colname] if str(x) != 'nan']
                            query_column = [x for x in df_dummy[col] if str(x) != 'nan']
                            if len(parent_column) >= len(query_column):
                                process_length = len(query_column)
                            if len(parent_column) <= len(query_column):
                                process_length = len(parent_column)
                            for i in range(0,process_length):
                                     if(parent_column[i]==query_column[i]):
                                          list_output.append("wild")
                                     else:
                                          list_output.append("mutate")
                            
                            df_dummy[new_col] = pd.Series(list_output)
                            warnings.simplefilter(action='ignore', category=FutureWarning)
                            del list_output
                    df_mismatched = df_dummy.loc[np.array([df_dummy[col].str.contains("mutate", na=False) for col in df_dummy.columns if is_string_dtype(df_dummy[col])]).any(axis=0)]
                    df_final = df_mismatched.drop(df_mismatched.filter(regex= 'new').columns, axis=1)
                    return df_final
		
        
        
        
        d = {}
        final_list  = []
        list_truncation = [] 
        n = 0
        while n < length :
              df_v = df[df.columns[n:n+total_gen]]  #subsetting dataframe passing protein from each generation 
              if len(df_v.columns) == total_gen:
                 similarity_index = []
                 p_name = df_v.filter(regex= 'SP11',axis=1)
                 if len(p_name.columns) == 1:
                    colname = p_name.columns[0]
                    for col in df_v.columns:
                            if col == colname:  
                                     pass
                            else:
                                     list_output = []
                                     parent_column = [x for x in df_v[colname] if str(x) != 'nan']
                                     query_column = [x for x in df_v[col] if str(x) != 'nan']
                                     if len(parent_column) >= len(query_column):
                                                    process_length = len(query_column)
                                     if len(parent_column) <= len(query_column):
                                                    process_length = len(parent_column)
                                     for i in range(process_length):
                                           if(parent_column[i]==query_column[i]):
                                                        list_output.append("True")
                                           else:
                                                        list_output.append("False")
                                     percent_identity = list_output.count("False")/len(list_output)                
                                     similarity_index.append(percent_identity)
                 if len(similarity_index) >= 1:
                  if similarity_index[0] <= 0.04 and similarity_index[1] <= 0.04 and similarity_index[2] <= 0.04 and similarity_index[3] <= 0.04 and similarity_index[4] <= 0.04 and similarity_index[5] <=0.04:
                   list_truncation.append(length_protein(df_v, Wild_id))
                   l = list(df_v.columns)
                   df3.drop(l, axis=1, inplace=True)
                   df_mutation = mutation(df_v, Wild_id)
                   df_transpose = df_mutation.T
                   row1 = df_mutation.columns.values.tolist()
                   row2 = df_mutation.values.tolist()
                   if len(row2) >= 1:
                    d[str(row1)] = str(row2)
                    final_list.append(row1)
                    final_list.append(row2) 
                   n += (total_gen-1)      
              n += 1 
        

        n = 0
        list_methionine = []
        while n <= length:
                df_subset = df3[df3.columns[n:n+total_gen]] #subsetting dataframe passing protein from each generation 
                if len(df_subset.columns) == total_gen:
                   list_methionine.append(truncation(df_subset, Wild_id))
                   #n += 4 
                n += 1
                
        z = []
        for value in list_methionine:
            if len(value) >=1:
               z.append(value)
        df_methionine = pd.DataFrame(z)
        if len(df_methionine.columns) > 1:
          #print(df_methionine)
          df_methionine=df_methionine.groupby(0, as_index=False).agg({col: list for col in df_methionine.columns[1:]})
          for col in df_methionine.columns:
            if col != 0:
             df_methionine[col] = [','.join(map(str, x)) for x in df_methionine[col]]

        df_point = pd.Series(d).to_frame()
        length_trunc = []
        for value in list_truncation:
            if len(value) != 0:
                length_trunc.append(value)
        df_trunc2 = pd.DataFrame(length_trunc)

        df_point.to_csv("Point_mutation.csv")
        df_methionine.to_csv("missed_meth.csv")
        df_trunc2.to_csv("Truncation.csv")     
