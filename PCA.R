library(rcdk)
library(ggplot2)
library(plotly)

# --------
# Load data: EDIT to path with your data!
# --------

inFile_datapath = '/Users/marianagonzmed/Downloads/PBOX.csv'
compounds2 <- read.csv(inFile_datapath, header = TRUE)

# --------
# Save data and plots: EDIT path to save in your directory
# --------

#data
PCAresults_path = '/Users/marianagonzmed/Downloads/PCA_results.csv'
PCA_PROPS_WEIGHTS_path = '/Users/marianagonzmed/Downloads/weights_pca_results.csv'
PCA_summary_path = '/Users/marianagonzmed/Downloads/pca_sumary_results.csv'

#plot
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
# remove NAS
descriptors <- na.omit(descriptors[,-(7:8)])
#remove row names
row.names(descriptors) <-seq.int(nrow(descriptors))

# --------
# PCA 
# --------

# Function to remove properties with stdev of 0
SDcutoff <- function(block, SDval) {
  standev <- as.data.frame(sapply(block, sd))
  names(standev) <- "column"
  standev <- subset(standev, column > SDval)
  return(block[,rownames(standev)])
}

# apply sdev func
descriptors2 <- SDcutoff(descriptors,0)
# apply PCA
pca <- prcomp(descriptors2, center = TRUE, scale. = TRUE) 
# get dataframe with SMILES, principal components and descriptors
results <- cbind(t_data, pca$x[,1:6],descriptors)
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

# export results
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