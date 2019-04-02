## probs

library(readr)

# read original example data (SNPTEST software) data with all decimals
SNPTEST <- readr::read_table2("inst/exdata/SNPTESToriginal.probs", col_names = FALSE) # white space separated
SNPTEST <- as.data.frame(SNPTEST)
SNPTEST[,5] <- as.character(SNPTEST[,5])
SNPTEST[,5] <- ifelse(SNPTEST[,5]=="TRUE","T",SNPTEST[,5])

# select 50 first SNPs
SNPTEST <- SNPTEST[1:50,]

# put 4 decimals
for (i in 6:ncol(SNPTEST)) SNPTEST[,i] <- round(SNPTEST[,i], 4)

cases <- SNPTEST[,c(1:5, 6:1505)]
controls <- SNPTEST[,c(1:5, 1506:3005)]

write.table(SNPTEST, file="inst/exdata/SNPTEST.probs", sep=" ", row.names = FALSE, col.names=FALSE)


save(cases, controls, file="./data/SNPTEST.rda",compress = TRUE)



# 
# 
# 
# 
# # write to probs format (all decimals)
# write.table(SNPTEST, file="inst/exdata/SNPTEST.probs", sep=" ", col.names=FALSE, row.names=FALSE)
# 
# # write to fst format (all decimals)
# fst::write_fst(SNPTEST, "./inst/exdata/SNPTEST.fst",uniform_encoding = TRUE,compress = 90)
# 
# # round to 4 decimals
# for (j in 6:ncol(SNPTEST)) SNPTEST[,j] <- round(as.double(SNPTEST[,j]),4)
# 
# # write to fst format (4 decimals)
# fst::write_fst(SNPTEST, "./inst/exdata/SNPTEST2.fst",uniform_encoding = TRUE,compress = 90)
# 
# 
# # write to probs format (4 decimals)
# write.table(SNPTEST, file="inst/exdata/SNPTEST2.probs", sep=" ", col.names=FALSE, row.names=FALSE)
# 
# 
# ## file info
# file.info("./inst/exdata/SNPTESToriginal.probs")
# file.info("./inst/exdata/SNPTEST.probs")
# file.info("./inst/exdata/SNPTEST.fst")
# file.info("./inst/exdata/SNPTEST2.probs")
# file.info("./inst/exdata/SNPTEST2.fst")
# 
# save(SNPTEST, file="./inst/exdata/SNPTESTcompressed.rda", compress = TRUE)
# save(SNPTEST, file="./inst/exdata/SNPTEST.rda")
# 
# 



# 
# 
# 
# # check
# 
# xxx <- fst::read.fst("./inst/exdata/SNPTEST.fst")
# yyy <- read_table2("inst/exdata/SNPTEST.probs", col_names = FALSE)
# as.data.frame(xxx)[1:10,1:10]
# as.data.frame(yyy)[1:10,1:10]
# 
# xxx <- fst::read.fst("./inst/exdata/SNPTEST2.fst")
# yyy <- read_table2("inst/exdata/SNPTEST2.probs", col_names = FALSE)
# as.data.frame(xxx)[1:10,1:10]
# as.data.frame(yyy)[1:10,1:10]
# 

