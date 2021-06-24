#!/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript --vanilla

library(data.table)

results <- data.table()

for (i in seq(1,1596)) {

	file <- paste0("outfiles/result.",i,".rdat")
#	if (file.exists(file)) {
		curr.table <- readRDS(file)
		results <- rbind(results, curr.table)
#	}

}

saveRDS(results, "models.rdat")
