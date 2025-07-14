t1 = Sys.time()
setwd('/gpfs52/data/g_gamazon_lab/zhoud2/data/cmap/raw/LINCS/LINCS2_beta/')

.libPaths("~/R/rlib-4.0.5")
library(cmapR)
args = as.numeric(commandArgs(TRUE))
print('start parsing!')
file = parse_gctx('level5_beta_trt_cp_n720216x12328.gctx')
print('parsing is ok!')

rows = file@rdesc
cols = file@cdesc
exp_mat = file@mat
write.table(rows, file = 'parse_gctx/row.txt', sep = '\t', row.names = F)
print('row is saved')
write.table(cols, file = 'parse_gctx/col.txt', sep = '\t', row.names = F)
print('col is saved')
write.table(exp_mat, file = 'parse_gctx/exp_mat.csv', sep = ',', row.names = F)
print('exp mat is saved')
t2 = Sys.time()
print(t2-t1)
