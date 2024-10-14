import pandas as pd
import numpy as np
from scipy.stats import ks_2samp

norm_count = pd.read_csv('/group/pmc020/binsawang/gene_project/csv/v2/unified_data.csv', low_memory=False)
norm_count.set_index(norm_count.columns[0], inplace=True)
norm_count.index.name = None

#### Dynamic Input ####
min_subtype = 1
subtype = 'All'
#######################

if subtype == 'All':
    # Use all data
    df_filtered = norm_count
else:
    # Only keep columns where the subtype is 'subtype'
    df_filtered = norm_count.loc[:, norm_count.iloc[-1] == subtype]

print('Subtype:', subtype)
print('Dim df:', df_filtered.shape)

# Get survival time, status, and subtype from the dataframe
survival_time = pd.to_numeric(df_filtered.iloc[-3].values, errors='coerce')
# status = pd.to_numeric(norm_count.iloc[-2].values, errors='coerce')
# subtypes = df_filtered.iloc[-1].values
# subtype_row = df_filtered.loc['PAM50']
# subtype_count = subtype_row.value_counts().to_dict()

gene_list = df_filtered.index[:-3]
results = []

for gene in gene_list:

    best_threshold = None
    best_p_value = np.inf
    best_ks_stat = -np.inf
    size_le = None
    size_he = None

    # Get gene count
    row = pd.to_numeric(df_filtered.loc[gene].values, errors='coerce')
    
    # Iterate over count values to find the best threshold
    for threshold in row:
        
        # Split time based on current threshold
        time_g1 = survival_time[row <= threshold]
        time_g2 = survival_time[row > threshold]

        # Perform KS test
        if len(time_g1) >= min_subtype and len(time_g2) >= min_subtype:
            ks_test = ks_2samp(time_g1, time_g2)
            ks_stat = round(ks_test.statistic, 4)
            ks_pval = round(ks_test.pvalue, 4)
        
            # Update best threshold if current ks_pval is smaller
            if ks_pval < best_p_value:
                best_p_value = ks_pval
                best_ks_stat = ks_stat
                best_threshold = round(threshold, 6)
                size_le = len(time_g1)
                size_he = len(time_g2)
    
    # Store the results in the list of dictionaries
    results.append({
        'gene': gene,
        f'threshold_{subtype}': best_threshold,
        f'pvalue_{subtype}': best_p_value,
        f'ks_statistic_{subtype}': best_ks_stat,
        'size_le': size_le,
        'size_he': size_he
    })

# Checking
print(min_subtype)
print(subtype)

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)
results_df.to_csv(f'/group/pmc020/binsawang/gene_project/result3/ks_{subtype}_{min_subtype}.csv', index=False)

#################################################################################################
# import pandas as pd
# import numpy as np
# from scipy.stats import ks_2samp

# # Read .csv files
# norm_data = pd.read_csv('/group/pmc020/binsawang/gene_project/csv/v2/norm_data.csv')
# meta_data = pd.read_csv('/group/pmc020/binsawang/gene_project/csv/v2/meta_data.csv')

# norm_data.set_index(norm_data.columns[0], inplace=True)
# norm_data.index.name = None
# norm_data.columns = norm_data.columns.astype(int)

# # Creat mapping dictionary
# key_to_time = meta_data.set_index('sample_id_key')['survival_time'].to_dict()
# key_to_time = {int(key): int(value) for key, value in key_to_time.items()}
# # key_to_status = meta_data.set_index('sample_id_key')['status'].to_dict()
# key_to_subtype = meta_data.set_index('sample_id_key')['subtype'].to_dict()

# # Create a dictionary to for gene id
# gene_id = norm_data.index
# gene_to_index = {gene: idx for idx, gene in enumerate(gene_id)}

# #### Dynamic Input ####
# min_group_size = 1
# subtype = 'LumA'
# results = []

# # Extract column names where the value is subtype
# col_to_keep = [key for key, value in key_to_subtype.items() if value == subtype]

# df_filtered = norm_data[col_to_keep]
# print('The dimensions of the filtered DataFrame are', df_filtered.shape)

# # Convert to array
# data = df_filtered.values
# columns_array = df_filtered.columns.to_numpy()

# for gene in gene_id:
    
#     # Get the index for the gene
#     gene_idx = gene_to_index[gene]

#     # Get the row where index = gene
#     row = data[gene_idx]
    
#     # Initialize variables
#     optimal_threshold = None
#     lowest_pval = float('inf')
#     highest_ks = float('-inf')
#     size_le = None
#     size_he = None
    
#     # Loop through each value in the row to find the optimal threshold
#     for threshold in row:
        
#         # Create masking
#         le_mask = row <= threshold
#         gt_mask = row > threshold

#         # Extract column indices for each group
#         group1_col = columns_array[le_mask]
#         group2_col = columns_array[gt_mask]
        
#         if len(group1_col) >= min_group_size and len(group2_col) >= min_group_size:
        
#             # Extract the survival time where index matches the column names
#             time_group1 = [key_to_time.get(key) for key in group1_col]
#             time_group2 = [key_to_time.get(key) for key in group2_col]

#             # Perform KS test
#             ks_test = ks_2samp(np.array(time_group1), np.array(time_group2))
#             ks_stat = round(ks_test.statistic, 4)
#             ks_pval = round(ks_test.pvalue, 4)

#             # Check if this p-value is the lowest
#             if ks_pval < lowest_pval:
#                 highest_ks = ks_stat
#                 lowest_pval = ks_pval
#                 optimal_threshold = threshold
#                 size_le = len(time_group1)
#                 size_he = len(time_group2)

#     # Store the results in the dictionary
#     results.append({
#         'gene': gene,
#         'threshold': optimal_threshold,
#         'pvalue': lowest_pval,
#         'ks_statistic': highest_ks,
#         'size_le': size_le,
#         'size_he': size_he
#     })

# # Checking
# print(min_group_size)
# print(subtype)

# # Convert the results to a DataFrame
# results_df = pd.DataFrame(results)
# results_df.to_csv('./result2/ot_all_ks_nomin.csv', index=False)