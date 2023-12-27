total.looms <- grep(".loom$",list.files("../../CellRanger-Outputs/",recursive = TRUE),value = TRUE)
grep("E_SI",total.looms,value = TRUE)

total.rds <- grep("rds$",list.files("../SingleSample/",recursive = TRUE),value = TRUE)
grep("^E_SI",total.rds,value = TRUE)

# remove bad samples
bad.samples <- c("E_LI_D135","E_SI_D135","E_SI_D135_B4")

cat(paste("\"../../CellRanger-Outputs/",grep("E_SI",total.looms,value = TRUE),"\",\n",sep = ""))
