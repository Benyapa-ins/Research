import pandas as pd
import numpy as np

norm_count = pd.read_csv('/group/pmc020/binsawang/gene_project/csv/v2/unified_data.csv', low_memory=False)
norm_count.set_index(norm_count.columns[0], inplace=True)
norm_count.index.name = None

#### Dynamic Input ####
min_subtype = 1
subtype = 'Normal'
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
epsilon = 1e-10

for gene in gene_list:

    # Get gene count
    row = pd.to_numeric(df_filtered.loc[gene].values, errors='coerce')

    best_threshold = None
    highest_jsd = -np.inf
    size_le = None
    size_he = None

    # Iterate over count values to find the best threshold
    for threshold in row:
        # Split time based on current threshold
        time_g1 = survival_time[row <= threshold]
        time_g2 = survival_time[row > threshold]

        # Perform KS test
        if len(time_g1) >= min_subtype and len(time_g2) >= min_subtype:
            # Combine both time groups to get all unique time values
            # all_times = list(set(time_g1 + time_g2))
            all_times = set(time_g1).union(set(time_g2))
            # Calculate the probability distribution for each group
            p = pd.Series(time_g1).value_counts(normalize=True).reindex(all_times, fill_value=0)
            q = pd.Series(time_g2).value_counts(normalize=True).reindex(all_times, fill_value=0)
            p = np.clip(p, epsilon, None)
            q = np.clip(q, epsilon, None)
            
            # Compute the average distribution
            m = 0.5 * (p + q)
            # Calculate KL Divergences
            KL_P_M = np.sum(p * np.log(p / m))
            KL_Q_M = np.sum(q * np.log(q / m))
            # Calculate JS divergence
            JSD = 0.5 * (KL_P_M + KL_Q_M)

            # Check if this JSD is the highest
            if JSD > highest_jsd:
                highest_jsd = JSD
                best_threshold = threshold
                size_le = len(time_g1)
                size_he = len(time_g2)

    # Store the results in the dictionary
    results.append({
        'gene': gene,
        'threshold': best_threshold,
        'JSD': highest_jsd,
        'size_le': size_le,
        'size_he': size_he
    })

# Checking
print(min_subtype)
print(subtype)

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)
results_df.to_csv(f'/group/pmc020/binsawang/gene_project/result2/jsd_{subtype}_{min_subtype}.csv', index=False)