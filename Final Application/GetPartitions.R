library(data.table)

idx <- seq(1, 750, by = 3)

outfile <- "used_partitions.csv"

if (file.exists(outfile)) file.remove(outfile)

for (i in 1:200){
  
  file3 <- paste0("Final Application/Final Draws_cont/Final_Draws_2cont", i, ".csv")
  
  draws3 <- fread(
    file3,
    nrows = 750
  )
  
  part3 <- draws3[idx, 3:ncol(draws3)]
  fwrite(part3, outfile, append = TRUE)
  
  rm(draws3, part3)
  gc()
  
  file4 <- paste0("Final Application/Final Draws_cont/Final_Draws_4cont", i, ".csv")
  
  draws4 <- fread(
    file4,
    nrows = 750
  )
  
  part4 <- draws4[idx, 3:ncol(draws4)]
  fwrite(part4, outfile, append = TRUE)
  
  rm(draws4, part4)
  gc()
  
  cat("Finished i =", i, "\n")
}
