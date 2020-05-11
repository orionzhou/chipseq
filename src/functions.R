require(devtools)
load_all("~/git/rmaize")
require(ggforce)
lib = 'chipseq'
dirp = sprintf("~/projects/%s", lib)
dird = file.path(dirp, 'data')
t_cfg = read_gspread_master(lib = c("chipseq",'dapseq'))

