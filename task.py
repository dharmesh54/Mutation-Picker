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
	# List all files in the specified directory
        files=os.listdir(path)

	# Function to read a text file and convert it into a dictionary
        def dictionery(file, path, filename):
           "reading text file to dictionery"
           protein = {}
           for i in files:   ## files in the working directory
             if i  == filename:
                ip = os.path.join(path,i)
                with open(ip) as f:
                  for line in f:
                        key = line.strip() # Remove leading/trailing whitespace
                        value = f.readline().strip() # Read the next line as value
                        protein[key] = list(value) ##apending key:value pair to the directory
                  return protein    #converting protein text file into the dictionery 

	# Convert the protein dictionary into a DataFrame
        dict_protein = dictionery(files, path, file_name)
        length = len(dict_protein)
        df=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in dict_protein.items()]))
        df3 =  df.copy(deep=True)  # Convert the protein dictionary into a DataFrame
        
        # Function to identify sequences with different lengths compared to the wild type
        def length_protein(df_s, identifier):
                list_output = []
                p_name = df_s.filter(regex= identifier,axis=1)   # Filter columns matching the identifier
                if len(p_name.columns) == 1:     # Check if exactly one column matches the identifier (e.g., the wild type)
                  colname = p_name.columns[0]   # Store the name of that single column (wild type)
		  # Iterate over all columns in the original dataframe
                  for col in df_s.columns:
                    if col == colname:      # Skip comparison if it's the wild type column itself
                        pass
                    else:
                        parent_column = [x for x in df_s[colname] if str(x) != 'nan']   # Get the wild type sequence, removing any 'nan' entries
                        query_column = [x for x in df_s[col] if str(x) != 'nan']     # Get the mutant sequence, also removing 'nan' entries  
                    if len(parent_column) != len(query_column):     # Compare lengths of wild type and mutant sequences
                        if colname not in list_output:     # Add wild type name to output if not already included
                           list_output.append(colname)
                        list_output.append(col)      # Add mutant column name and its sequence length
                        list_output.append(len(query_column))  
                return list_output

	# Function to detect truncation mutations based on methionine (M) positions
        def truncation(df_dummy, identifier):
                 list_truncated_mutant = []   # Initialize an empty list to store results for truncated mutants
                 p_name = df_dummy.filter(regex= identifier,axis=1)   # Filter the dataframe to get the column(s) matching the identifier (e.g., wild type)
                 if len(p_name.columns) == 1: # Ensure only one parent/wild-type column is selected
                  colname = p_name.columns[0]  # Get the name of the wild-type column
                  for col in df_dummy.columns: 
                    if col == colname:      # Skip comparison with the wild type itself
                        pass
                    else:
                        #list_output = []
                        parent_column = [x for x in df_dummy[colname] if str(x) != 'nan']  # Extract sequences from the wild type, removing any 'nan' values
                        query_column = [x for x in df_dummy[col] if str(x) != 'nan']    # Extract sequences from the mutant column, also removing 'nan' values
                        
			# Find all indices where Methionine ("M") occurs in the wild-type sequence
			M_position = [i for i,val in enumerate(parent_column) if val =="M"]   
			# Iterate through Methionine positions (skipping the first one)
                        for j in range(1, len(M_position)): 
                            list_of_match_amino_acids = [] 
			    # Align mutant sequence with wild type starting from the j-th Methionine position
                            for i in range(len(query_column)):
                                l = M_position[j]+i    # Shifted index in the wild-type sequence
                                if len(parent_column) < l:  # If shifted index exceeds wild-type sequence length, mark mismatch
                                    list_of_match_amino_acids.append("False")
                                if len(parent_column) > l:  # If still within bounds, compare amino acids
                                         if(parent_column[l]==query_column[i]): 
                                                list_of_match_amino_acids.append("True")
                                         else:
                                             list_of_match_amino_acids.append("False")

			    # Count number of mismatches
                            fal_eval = list_of_match_amino_acids.count("False") 
			    # Calculate mismatch ratio (False rate)
                            percent_identity = fal_eval/len(list_of_match_amino_acids)
                            if percent_identity <= 0.05:  # If identity is high (≤ 5% mismatch), we suspect truncation or shift
                               if colname not in list_truncated_mutant: # Record the wild-type column name (only once)
                                   list_truncated_mutant.append(colname)
			       # Record the mutant column name (if not wild-type)
                               if col != colname:
                                   list_truncated_mutant.append(col)
                               list_truncated_mutant.append(M_position[j])  # Record the Methionine start position indicating likely truncation
                               break
                 return list_truncated_mutant   # Return the list of columns where a truncation-like pattern is detected
                 
         # Function to identify point mutations in mutant sequences compared to a wild-type sequence        
        def mutation(df_dummy, identifier):
                    df_dummy.loc[:,'position'] = pd.Series(list(range(1,length +1))) # Add a 'position' column (1-based indexing) to track amino acid positions
                    p_name = df_dummy.filter(regex= identifier,axis=1)  # Identify the wild-type column using the identifier pattern
                    if len(p_name.columns) == 1: # Ensure exactly one wild-type column is selected
                      colname = p_name.columns[0] 
		      # Iterate over all other columns to compare with the wild-type    
                      for col in df_dummy.columns:
                        if col == colname or col == "position":  
                            pass
                        else:
                            new_col = "new" + col  # Initialize new column for marking mutation status
                            list_output = []
			    # Remove NaNs and extract sequences
                            parent_column = [x for x in df_dummy[colname] if str(x) != 'nan']
                            query_column = [x for x in df_dummy[col] if str(x) != 'nan']

			    # Compare up to the shortest length between wild-type and mutant
                            if len(parent_column) >= len(query_column):
                                process_length = len(query_column)
                            if len(parent_column) <= len(query_column):
                                process_length = len(parent_column)
			     # Label each position as "wild" (same) or "mutate" (different)
                            for i in range(0,process_length):
                                     if(parent_column[i]==query_column[i]):
                                          list_output.append("wild")
                                     else:
                                          list_output.append("mutate")
					     
                            # Add result column to dataframe
                            df_dummy[new_col] = pd.Series(list_output)
                            warnings.simplefilter(action='ignore', category=FutureWarning)  # Avoid future warnings on chained assignment
                            del list_output   # Clear list to free memory
                    df_mismatched = df_dummy.loc[np.array([df_dummy[col].str.contains("mutate", na=False) for col in df_dummy.columns if is_string_dtype(df_dummy[col])]).any(axis=0)]   # Identify rows containing any mutation ("mutate" label)
                    # Drop intermediate columns that start with 'new' (mutation markers)
		    df_final = df_mismatched.drop(df_mismatched.filter(regex= 'new').columns, axis=1)
                    return df_final
		
        
        
        
        d = {}   # Initialize containers to hold mutation results
        final_list  = []
        list_truncation = [] # Truncation results (based on full alignment)
        n = 0
	# Loop through the DataFrame `df` in segments (each segment represents a generation group)
        while n < length :
              df_v = df[df.columns[n:n+total_gen]]  #subsetting dataframe passing protein from each generation 
              if len(df_v.columns) == total_gen:  # Ensure the selected group has the correct number of sequences
                 similarity_index = []
		 # Identify the wild-type column (typically named using "SP11" or similar)
                 p_name = df_v.filter(regex= 'SP11',axis=1)
                 if len(p_name.columns) == 1:
                    colname = p_name.columns[0] 
		    # Compare each mutant column against the wild-type
                    for col in df_v.columns:
                            if col == colname:  
                                     pass
                            else:
                                     list_output = [] 
				    # Remove NaNs and extract sequences
                                     parent_column = [x for x in df_v[colname] if str(x) != 'nan']
                                     query_column = [x for x in df_v[col] if str(x) != 'nan']
                                     if len(parent_column) >= len(query_column):
                                                    process_length = len(query_column)
                                     if len(parent_column) <= len(query_column):
                                                    process_length = len(parent_column)
				     # Identify mutations by comparing residues
                                     for i in range(process_length):
                                           if(parent_column[i]==query_column[i]):
                                                        list_output.append("True")
                                           else:
                                                        list_output.append("False")
				     # Calculate mutation rate (as a percent identity)
                                     percent_identity = list_output.count("False")/len(list_output)                
                                     similarity_index.append(percent_identity)
		 # If the first 6 mutants are at least 96% similar (≤ 4% mismatch), continue analysis
                 if len(similarity_index) >= 1:
                  if similarity_index[0] <= 0.04 and similarity_index[1] <= 0.04 and similarity_index[2] <= 0.04 and similarity_index[3] <= 0.04 and similarity_index[4] <= 0.04 and similarity_index[5] <=0.04:
                   # Analyze truncation (length comparison)
		   list_truncation.append(length_protein(df_v, Wild_id))
                   l = list(df_v.columns)
                   df3.drop(l, axis=1, inplace=True)   # Drop analyzed columns from df3 to prevent re-analysis
                   df_mutation = mutation(df_v, Wild_id)  # Run mutation analysis for point mutations
                   df_transpose = df_mutation.T  # Transpose result for output formatting
                   row1 = df_mutation.columns.values.tolist()  # Column headers
                   row2 = df_mutation.values.tolist()  # Mutation values
                   if len(row2) >= 1:  # Save non-empty mutation result
                    d[str(row1)] = str(row2)
                    final_list.append(row1)
                    final_list.append(row2) 
		   # Skip to next non-overlapping group
                   n += (total_gen-1)      
              n += 1 
        

        n = 0
        list_methionine = []
        while n <= length:    # Loop through remaining columns in df3 in groups of total_gen
                df_subset = df3[df3.columns[n:n+total_gen]] #subsetting dataframe passing protein from each generation 
                # If subset has full group, check for missed Methionine-induced truncation
		if len(df_subset.columns) == total_gen:
                   list_methionine.append(truncation(df_subset, Wild_id))
                   #n += 4 
                n += 1
        # Filter out empty truncation results from methionine-based analysis        
        z = []
        for value in list_methionine:
            if len(value) >=1:
               z.append(value)
        df_methionine = pd.DataFrame(z)
	# Group by wild-type and concatenate corresponding mutants as comma-separated values
        if len(df_methionine.columns) > 1:
          #print(df_methionine)
          df_methionine=df_methionine.groupby(0, as_index=False).agg({col: list for col in df_methionine.columns[1:]})
          for col in df_methionine.columns:
            if col != 0:
             df_methionine[col] = [','.join(map(str, x)) for x in df_methionine[col]]
        # Convert dictionary of point mutations into DataFrame
        df_point = pd.Series(d).to_frame()

	# Filter out empty results from full-length truncation analysis
        length_trunc = []
        for value in list_truncation:
            if len(value) != 0:
                length_trunc.append(value)
        df_trunc2 = pd.DataFrame(length_trunc)

        df_point.to_csv("Point_mutation.csv")   # Contains detailed point mutation data
        df_methionine.to_csv("missed_meth.csv")  # Captures methionine-initiated truncations
        df_trunc2.to_csv("Truncation.csv")     # Captures length-based truncations
