##############################################
########### reshape tensorqtl output to z matrix ###########
##############################################
# load packages -----
library(data.table)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

inp = parameter[1]

## output -----
outp = parameter[2]


# reshape tensorqtl output to z matrix -----
z = fread(
  inp,
  header = TRUE, sep = "\t"
)

z.mat = dcast(z, snp~gene, value.var = "zscore", fun.aggregate = mean, drop = FALSE)

fwrite(
  z.mat, outp,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

