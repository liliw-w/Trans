rm(list = ls())
require(data.table)

## Download SU's info for each job
system("rcchelp usage --byjob > SU.txt")

## read SU file
Dat = fread("SU.txt")

## Check input
str(Dat)
head(Dat)
Dat[1:10, ]

## Data wrangling, add date, year, month, day
Dat$Start_Date = as.Date(Dat$Start)
Dat$Start_Y = format(Dat$Start_Date, "%y")
Dat$Start_M = format(Dat$Start_Date, "%m")
Dat$Start_D = format(Dat$Start_Date, "%d")

## Row slicing
Dat_Xuanyao_21 = Dat[Dat$Account == "pi-xuanyao" & Dat$Start_Y == "21", ]
min(Dat_Xuanyao_21[Dat_Xuanyao_21$Start_M == "08", Start_Date])
max(Dat_Xuanyao_21[Dat_Xuanyao_21$Start_M == "08", Start_Date])

s = as.Date("2021-07-18"); e = as.Date("2021-09-01")
ind_job = Dat_Xuanyao_21$Start_Date > s & Dat_Xuanyao_21$Start_Date < e
sum(ind_job)

## write out
res = Dat_Xuanyao_21[ind, ]
sum(res$"Charge (SU)")
fwrite(res, "SU.XY.21.M73.txt", quote=FALSE, sep="\t")
