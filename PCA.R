library(rcdk)
library(ggplot2)
library(plotly)
library(dplyr)


# --------
# Load data: EDIT to path with your data!
# --------

inFile_datapath = '/Users/marianagonzmed/Downloads/PBOX.csv'
compounds2 <- read.csv(inFile_datapath, header = TRUE)

# --------
# Save data and plots: EDIT path to save in your directory
# --------

# Descriptors data here we just need the path, don't add the file name
descriptors_path = '/Users/marianagonzmed/Downloads/stat-'

# Descriptors plot here we just need the path, don't add the file name
descriptors_plot_path = '/Users/marianagonzmed/Downloads/plot-'

# PCA data
PCAresults_path = '/Users/marianagonzmed/Downloads/PCA_results.csv'
PCA_PROPS_WEIGHTS_path = '/Users/marianagonzmed/Downloads/weights_pca_results.csv'
PCA_summary_path = '/Users/marianagonzmed/Downloads/pca_sumary_results.csv'

# PCA plot
PCA_2D_plot = '/Users/marianagonzmed/Downloads/pca_2D_plot.tiff'

# --------
# Descriptors I want: EDIT to include the descriptors you want
# --------

# List of all descriptors available 
descriptor_names <- c(get.desc.names("all"))
print(descriptor_names)

# clasic ones: MW, TopoPSA, nRotB, nHBDon, nHBAcc,ALogP
# EDIT these numbers ONLY look at output of print(descriptor_names)
descriptor_numbers <- c(9,12,14,27,28,50,23)

# -------
# Data-input processing
# -------

# third column with ID should not be numeric
ALPHAID <- as.data.frame(paste("ID_",1:nrow(compounds2),sep=""))
if (is.numeric(compounds2[,3])){
  compounds2[,3] <- ALPHAID 
  }

# we only need the first 3 columns
t_data <- compounds2[,1:3]
# we will call the columns by name so we rename them
colnames(t_data) <- c("SMILES","DB","ID")
# we reorder so compunds from same database are grouped
t_data <- t_data[order(t_data$DB),]

# --------
# Descriptors calculation, this part takes time if you have a lot of molecules 
# --------

# We choose clasic ones
descNames <- c()
for (p in descriptor_numbers) {
  descNames <- c(descNames, descriptor_names[p])
  }

print('selected descriptors:')
print(descNames)

mols2 <- parse.smiles(as.vector(t_data[,1]))
# evaluate molecules
descriptors <- eval.desc(mols2, descNames)
descriptors <- na.omit(descriptors)
# remove NAS
descriptors <- na.omit(descriptors[,-(7:8)])
# save column names
names_cols <- c(colnames(descriptors))
#remove row names
row.names(descriptors) <-seq.int(nrow(descriptors))
descriptors2 <- cbind(t_data, descriptors)

# ------
# Descriptors distribution
# ------

group_by_database <- subset(descriptors2, select=-c(SMILES,ID)) %>% group_by(DB)

stats_descriptors <- group_by_database %>%
  summarise_all(list(min = min, max = max, mean = mean, median = median,
                     Q1 = ~ quantile(x = ., probs = 0.25),
                     Q3 = ~ quantile(x = ., probs = 0.75),
                     sd=sd
                     ))

# remove DB name column
names_rows <- stats_descriptors$DB
stats_descriptors[1] <- NULL
row.names(stats_descriptors) <- names_rows

# divide long database
divide_num = length(stats_descriptors)/length(descriptor_numbers)
split_statistics <- split.default(stats_descriptors, rep(1:divide_num, divide_num))

# save data frames with database statistics, 
# calculated descriptors are exported with PCA results later on
for (i in 1:divide_num) {
  path_to_save = paste(descriptors_path, names_cols[i], '.csv', sep = "")
  df <- as.data.frame(split_statistics[i])
  row.names(df) <- names_rows
  write.csv(df, path_to_save)
  }

# ------
# Descriptors plots. The size of the text on the plots can be changed if needed
# ------

for (i in 1:divide_num) {
  index_property <- grep(names_cols[i], colnames(descriptors2))
  
  path_to_save_density = paste(descriptors_plot_path, names_cols[i],'density' ,'.tiff', sep = "")
  
  # Save plots
  
  tiff(path_to_save_density, units="in", width=6, height=5, res=300)
  print(ggplot(descriptors2, aes(descriptors2[,index_property])) +
          labs(x = names_cols[i], y = "Density") +
          geom_density(aes(fill=DB, colour=DB), alpha = 0.2) +
          theme(text = element_text(size=25),
                axis.text.y = element_text(vjust=1,size=18),
                axis.text.x = element_text(vjust=1,size=18),
                legend.title=element_blank(), 
                legend.key.size = unit(1.0,"cm"),
                legend.key = element_rect(fill = "white"),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black")))
  dev.off()
  
  # Histogram plot
  path_to_save_hist = paste(descriptors_plot_path, names_cols[i],'hist' ,'.tiff', sep = "")
  
  tiff(path_to_save_hist, units="in", width=6, height=5, res=300)
  print(ggplot(descriptors2, aes(descriptors2[,index_property])) +
          labs(x = names_cols[i], y = "Frequency") +
          geom_histogram(aes(fill=DB, colour=DB), alpha = 0.2, position = position_dodge()) +
          theme(text = element_text(size=25),
                axis.text.y = element_text(vjust=1,size=18),
                axis.text.x = element_text(vjust=1,size=18),
                legend.key.size = unit(1.0,"cm"),
                legend.title=element_blank(), 
                legend.key = element_rect(fill = "white"),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black")))
  dev.off()
  
  # Box plot
  path_to_save_bp = paste(descriptors_plot_path, names_cols[i],'pm' ,'.tiff', sep = "")
  
  tiff(path_to_save_bp, units="in", width=6, height=5, res=300)
  print(ggplot(descriptors2, aes(DB, descriptors2[,index_property])) +
          labs(x = "Data Sets", y = names_cols[i]) +
          geom_boxplot(aes(fill=DB)) +
          theme(text = element_text(size=25),
                axis.text.y = element_text(size=14),
                axis.text.x = element_text(vjust=1,size=14),
                legend.title=element_blank(), 
                legend.key.size = unit(1.0,"cm"),
                legend.key = element_rect(fill = "white"),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black")))
  dev.off()
  
  # Interactive plot

  m <- list(l = 110, r = 60, b = 100, t = 50)
  l <- list(font= list(size=11))
  
  # Density plot
  print(
    ggplotly(ggplot(descriptors2, aes(descriptors2[,index_property])) +
                   labs(x = names_cols[i], y = "Density") +
                   geom_density(aes(fill=DB, colour=DB), alpha = 0.2) +
                   theme(text = element_text(size=18),
                         axis.text.y = element_text(vjust=1,size=12),
                         axis.text.x = element_text(vjust=1,size=12),
                         legend.title=element_blank(), 
                         legend.key = element_rect(fill = "white"),
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"))
             ) %>% layout(autosize = F, width = 700, height = 700, margin = m, legend = l)
    )
  
  # Histogram
  print(
    ggplotly(ggplot(descriptors2, aes(descriptors2[,index_property])) +
                   labs(x = names_cols[i], y = "Frequency") +
                   geom_histogram(aes(fill=DB, colour=DB), alpha = 0.2, position = position_dodge()) +
                   theme(text = element_text(size=18),
                         axis.text.y = element_text(vjust=1,size=12),
                         axis.text.x = element_text(vjust=1,size=12),
                         legend.title=element_blank(), 
                         legend.key = element_rect(fill = "white"),
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"))
             ) %>% layout(autosize = F, width = 700, height = 700, margin = m, legend = l)
    )
  
  # Box plot
  print(
    ggplotly(ggplot(descriptors2, aes(DB, descriptors2[,index_property])) +
                   labs(x = "Data sets", y = names_cols[i]) +
                   geom_boxplot(aes(fill=DB)) +
                   theme(text = element_text(size=18),
                         axis.text.y = element_text(size=12),
                         axis.text.x = element_text(vjust=1,size=12),
                         legend.title=element_blank(), 
                         legend.key = element_rect(fill = "white"),
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"))
             ) %>% layout(autosize = F, width = 700, height = 700, margin = m, legend = l)
  )
  }

# --------
# PCA 
# --------

# Function to remove descriptors with stdev of 0
SDcutoff <- function(block, SDval) {
  standev <- as.data.frame(sapply(block, sd))
  names(standev) <- "column"
  standev <- subset(standev, column > SDval)
  return(block[,rownames(standev)])
}

# apply sdev func
descriptors3 <- SDcutoff(descriptors,0)
# apply PCA
pca <- prcomp(descriptors3, center = TRUE, scale. = TRUE) 
# get data frame with SMILES, principal components and descriptors
results <- cbind(descriptors2, pca$x[,1:6])
# contribution of each property to each component
loadns <- as.data.frame(pca$rotation)
#summary of PCA
sumpc <-summary(pca) 
sumpc <- do.call(rbind.data.frame,sumpc)
sumpc2 <- data.frame(c(sumpc[1,]))
sumpc2[nrow(sumpc2) + 1,] = c(sumpc[8,])
sumpc2[nrow(sumpc2) + 1,] = c(sumpc[9,])
sumpc3 <- tail(sumpc, n = 2L)
summary <- rbind(sumpc2,sumpc3)
rownames(summary)<- c("sdev","center","scale","Proportion.of.variance","Cumulative.Proportion")

# --------
# Export PCA results 
# --------

# export SMILES, DB, IDS, descriptors and PCA results
write.csv(results, PCAresults_path)
# export PCA summary
write.csv(loadns, PCA_PROPS_WEIGHTS_path)
# export properties contribution to each component
write.csv(summary, PCA_summary_path)

# ------
# Plot PCA
# ------

#   2D Plot

tiff(PCA_2D_plot, units="in", width=6, height=5, res=300)
print(ggplot(data = results, aes(x = PC1, y = PC2, colour = DB)) +
        geom_point(size=3) +
        theme(legend.title=element_blank(), 
              legend.key = element_rect(fill = "white"),
              legend.key.size = unit(1.5,"cm"),
              legend.text = element_text(size = 20),
              text = element_text(size=30),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black")))
dev.off()

# interactive 2D plot
m <- list(l = 50,r = 50,b = 50,t = 50,pad = 4)

print(
  ggplotly(
    ggplot(data = results, aes(x = PC1, y = PC2, colour = DB)) +
      geom_point(aes(text=paste("ID:",ID)), size = .9) +
      theme(legend.title=element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
  ) %>% layout(autosize = F, width = 700, height = 700, margin = m)
)

# interactive 3D plot
print(plot_ly(results, x = ~PC1, y = ~PC2, z = ~PC3, color = ~DB,
        text = results$ID, width = 800, height = 800) %>%
  add_markers(marker=list(size=4)) %>%
  layout(scene = list(xaxis = list(title = 'PCa'),
                      yaxis = list(title = 'PCb'),
                      zaxis = list(title = 'PCc'))))

# ------
# Properties distance calculation
# ------

# all against all
row.names(descriptors) <- compounds2$Compound.ID
# matrix with similarity values
matrix_distances_all = as.matrix(dist(descriptors, method = "euclidean", diag = FALSE, upper = FALSE, p = 2))
# remove upper matrix to avoid duplicates
upper.tri(matrix_distances_all, diag = FALSE)
matrix_distances_all[upper.tri(matrix_distances_all)] <- NA
# matrix to data frame
matrix_distances_all2 <- as.data.frame(as.table(matrix_distances_all))
# remove NAs
matrix_distances_all2 <- na.omit(matrix_distances_all2)

# add SMILES and DB columns
colnames(compounds2) <- c("SMILES","DB","Var1")
CHECK <- merge(matrix_distances_all2, compounds2, by = "Var1", all.x=TRUE, all.y=FALSE)
colnames(compounds2) <- c("SMILES","DB","Var2")
CHECK <- merge(CHECK, compounds2, by = "Var2", all.x=TRUE, all.y=FALSE)

# INTRA database distance


# INTER database distance

# ------
# Properties distance plot
# ------


# ------
# Smiles distance calculation
# ------


# ------
# Smiles distance plot
# ------


# ------
# Scaffolds curves
# ------

# ------
# Scaffolds Shannon Entropy
# ------

# ------
# Fingerprints curves
# ------