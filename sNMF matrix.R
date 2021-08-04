
#Read in the ClumppIndFile files produced by Clumpak for each value of k
#these are the assignment values for the major modes only

setwd("C:\\Users\\zshaik\\Desktop\\clumpak\\all")

list <- as.list(list.files(pattern = "\\.output"))
files <- list()
n <- c(1:length(list))
for (i in n) {
  files[[i]] <- read.table(list[[i]])
}
names(files) <- list

#drop unnecessary columns
for (i in n) {
  files[[i]] <- files[[i]][, -c(1:5)]
}

##find the column with the highest ancestry coefficient for each accession (row)
for (i in n) {
  files[[i]]$max <- colnames(files[[i]])[max.col(files[[i]], ties.method = "first")]
}

maxs <- list()
for (i in n) {
  maxs[[i]] <- data.frame(files[[i]]$max)  
}

#all accessions are assigned to one ancestral pool for k = 1, so it should not be included 
df <- cbind(maxs[[2]], maxs[[3]], maxs[[4]], maxs[[5]],
            maxs[[6]], maxs[[7]], maxs[[8]], maxs[[9]], maxs[[10]],
            maxs[[11]], maxs[[12]], maxs[[13]], maxs[[14]], maxs[[15]],
            maxs[[16]], maxs[[17]], maxs[[18]], maxs[[19]], maxs[[20]])

library(stringr)
names <- str_remove(unlist(list), ".output")
names <- names[2:20]
colnames(df) <- names

##calculate distance matrix
#most columns have >2 levels of a nominal categorical variable,
#dummify to make compatible with daisy

library(fastDummies)
dum_df <- dummy_cols(df)
dum_df <- dum_df[,20:222]
colnames(dum_df)

library(cluster)
gower.dist <- as.matrix(daisy(dum_df, metric = c("gower"), 
                              type = list(asymm = c(1:203))))
gower.dist <- round(gower.dist, 2)

library(plot.matrix)
windows()
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(gower.dist, border = NA, 
     col = gray.colors(100, start=0, end = 1), 
    axis.col=NULL, axis.row=NULL, 
     xlab='', main = "",
     key=list(side=3, cex.axis=1, font = 2),
     fmt.key="%.2f", 
     polygon.key=NULL, axis.key=NULL, spacing.key=c(3,2,2),
     ylab = "Population")

##make custom population labels
axis(side = 1, at = c(seq(0, 288, 6)), cex.axis = 0.7,
     labels = c("","2", "1", "3", "9", "10", "11", "13", "28", "4", "5", "6", "7", "8", "14", "15", "16", "24", "25", "46", "47", "17", "18", "20", "21", "26", "27", "22", "30", "44", "48", "19", "23", "29", "31", "42", "43", "33", "35", "36", "37", "38", "39", "40", "41", "34", "45", "32", ""), 
     lwd = 0.5, las = 2)




