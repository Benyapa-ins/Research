import pandas as pd
from lifelines import CoxPHFitter

# Read .csv files
norm_data = pd.read_csv('/Users/tk/Desktop/Research-CITS5014/csv/norm_data.csv')
meta_data = pd.read_csv('/Users/tk/Desktop/Research-CITS5014/csv/meta_data.csv')

norm_data.set_index(norm_data.columns[0], inplace=True)
norm_data.index.name = None
norm_data.columns = norm_data.columns.astype(int)

# Filter norm data to get indices where subtype is not NaN
valid_id = meta_data[meta_data['subtype'].notna()]['sample_id_key'].to_list()
norm_data = norm_data.loc[:, valid_id]
meta_data = meta_data[meta_data['sample_id_key'].isin(valid_id)]

# Creat mapping dictionary
key_to_time = meta_data.set_index('sample_id_key')['survival_time'].to_dict()
key_to_time = {int(key): int(value) for key, value in key_to_time.items()}
key_to_status = meta_data.set_index('sample_id_key')['status'].to_dict()
key_to_subtype = meta_data.set_index('sample_id_key')['subtype'].to_dict()

time_row = {col: key_to_time.get(col, None) for col in norm_data.columns}
df = pd.concat([pd.DataFrame(time_row, index=['time']), norm_data], ignore_index=False)

event_row = {col: key_to_status.get(col, None) for col in norm_data.columns}
df = pd.concat([pd.DataFrame(event_row, index=['event']), df], ignore_index=False)

df_cox = df.transpose()
# Remove genes where all values are 0
df_cox = df_cox.loc[:, (df_cox != 0).any(axis=0)]

cph = CoxPHFitter()
cph.fit(df_cox, duration_col='time', event_col='event')
cph.print_sumsmsary()