from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import numpy as np
import os
import math
# os.chdir('/gpfs52/data/g_gamazon_lab/zhoud2/data/cmap/raw/LINCS/LINCS2_beta/')

gctx_file = "level5_beta_trt_cp_n720216x12328.gctx"
print('start parsing')
gctx_data = parse(gctx_file)
print('parsing is ok')
# Access the expression matrix
exp_mat = gctx_data.data_df
exp_mat = pd.DataFrame(exp_mat)

# Access the sample information
pd.DataFrame(exp_mat.columns).to_csv('parse_gctx/col.txt', index = 0, sep = '\t')
# Access the row (gene) information
pd.DataFrame(exp_mat.index).to_csv('parse_gctx/row.csv', index = 0, sep = '\t')

ncol = len(exp_mat.columns)
group_k = math.ceil(ncol / 100)
for i in range(0, 100):
  print(i)
  if i == 99:
    group_col = exp_mat.columns[(i*group_k):]
  else:
    group_col = exp_mat.columns[(i*group_k):((i + 1)*group_k)]
  group_mat = exp_mat[group_col]
  group_mat.to_csv('parse_gctx/exp_mat_{}.csv'.format(i), index = 0, sep = ',')
#np.save('exp_mat.npy', expression_matrix)
#exp_mat.to_csv('parse_gctx/exp_mat.csv', index = 0, sep = ',')

print('ok')
