library(limma)
library(dplyr)
library(stringr)

FILE <- "RAMOS_ProteinGroup.tsv"
df <- read.delim(FILE, row.names = 1, sep = "\t")
cat("Before decontamination:", dim(df), "\n")
contaminants <- c("P02768", "P0DUB6", "P0DTE7", "P0DTE8", "P01008", "P08758", 
                  "P61769", "P55957", "P00915", "P00918", "P04040", "P07339", 
                  "P08311", "P01031", "P02741", "P00167", "P99999", "P01133", 
                  "P05413", "P06396", "P08263", "P09211", "P69905", "P68871", 
                  "Q99075", "P01344", "P10145", "P08476", "P09529", "P06732", 
                  "P00709", "P41159", "P61626", "P02144", "Q15843", "P15559", 
                  "P16083", "P01127", "P62937", "Q06830", "P01112", "P02753", 
                  "P62979", "P00441", "P63165", "P12081", "P10636", "P10599", 
                  "P01375", "P02787", "P02788", "P51965", "O00762", "Q7Z3Y9", 
                  "P12035", "P19012", "Q5XKE5", "P04259", "Q04695", "P19013", 
                  "Q86Y46", "Q3SY84", "P08729", "P05787", "Q6A163", "P05783", 
                  "Q7Z3Y7", "P08779", "P04264", "P13647", "Q6KB66", "Q7Z794", 
                  "O95678", "P02538", "Q14CN4", "Q2M2I5", "P08727", "Q7RTS7", 
                  "Q7Z3Y8", "P35900", "P35527", "Q9C075", "Q99456", "P02533", 
                  "P48668", "P13645", "P13646", "Q01546", "Q7Z3Z0", "Q8N1N4", 
                  "P35908", "Q6A162", "P78386", "O76015", "O76011", "O43790", 
                  "Q92764", "Q14525", "Q14533", "O76014", "O76009", "P78385", 
                  "Q15323", "Q9NSB4", "Q14532", "O76013", "Q9NSB2")

# Filter out contaminants
df <- df[!rownames(df) %in% contaminants, ]
cat("After decontamination:", dim(df), "\n")

df <- df[df[, 1] != "", ] 
rownames(df) <- df[, 1]
df <- df[, -1]
head(df)

# Define group columns
group_columns <- list(
  DMSO = grep("1500ngRD", colnames(df), value = TRUE) %>% .[!str_detect(., "TR1")],
  PS1_200nM = grep("1500ngR1_", colnames(df), value = TRUE) %>% .[!str_detect(., "TR1")],
  PS2_200nM = grep("1500ngR2", colnames(df), value = TRUE) %>% .[!str_detect(., "TR1")],
  PS3_200nM = grep("1500ngR3", colnames(df), value = TRUE) %>% .[!str_detect(., "TR1")],
  PS5_200nM = grep("1500ngR5", colnames(df), value = TRUE) %>% .[!str_detect(., "TR1")],
  PS10_200nM = grep("1500ngR10", colnames(df), value = TRUE) %>% .[!str_detect(., "TR1")]
)
# group_columns <- list(
#   DMSO = grep("1500ng231_D", colnames(df), value = TRUE),
#   PS1_200nM = grep("1500ng231_1_", colnames(df), value = TRUE),
#   PS2_200nM = grep("1500ng231_2", colnames(df), value = TRUE),
#   PS3_200nM = grep("1500ng231_3", colnames(df), value = TRUE),
#   PS5_200nM = grep("1500ng231_5", colnames(df), value = TRUE),
#   PS10_200nM = grep("1500ng231_10", colnames(df), value = TRUE)
# )

# Print group columns and their lengths
print(group_columns)
for (key in names(group_columns)) {
  cat(key, ": ", length(group_columns[[key]]), "\n")
}

# Extract data for each group into a matrix for limma
expression_data <- df[, unlist(group_columns)]
# Log-transform the data (log2 transformation is common in proteomics)
#expression_data_log <- log2(expression_data + 1)  # Adding 1 to avoid log(0) issues
expression_data_log <- log2(as.matrix(expression_data) + 1)
# Scale the data (z-score normalization)
expression_data_scaled <- scale(expression_data_log)
# Visualizing the log-transformed data
boxplot(expression_data_log, main = "Boxplot of Log-Transformed Data")
hist(expression_data_log, breaks = 50, main = "Histogram of Log-Transformed Data")
# Visualizing the scaled data
boxplot(expression_data_scaled, main = "Boxplot of Scaled Data")
hist(expression_data_scaled, breaks = 50, main = "Histogram of Scaled Data")

expression_data <- expression_data_log  # Don't use _scaled. It leads to skewed FC ratio
expression_data <- cbind(PG.Genes = df$PG.Genes, expression_data)


# Create a factor for each sample's condition
group_factor <- factor(rep(names(group_columns), sapply(group_columns, length)))

# Set up the design matrix, with DMSO as the reference
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- sub("group_factor", "", colnames(design))
print(design)
# Use limma to fit the model
fit <- lmFit(expression_data, design)

# Define contrasts for each treatment group vs. DMSO
contrast_matrix <- makeContrasts(
  PS1_200nM_vs_DMSO = PS1_200nM - DMSO,
  PS2_200nM_vs_DMSO = PS2_200nM - DMSO,
  PS3_200nM_vs_DMSO = PS3_200nM - DMSO,
  PS5_200nM_vs_DMSO = PS5_200nM - DMSO,
  PS10_200nM_vs_DMSO = PS10_200nM - DMSO,
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Example: Get results for each contrast
results_list <- list()
for (contrast in colnames(contrast_matrix)) {
  results_list[[contrast]] <- topTable(fit2, coef = contrast, adjust.method = "BH", number = Inf)
}

# Print results for each contrast
lapply(names(results_list), function(contrast) {
  cat("\nResults for", contrast, ":\n")
  print(head(results_list[[contrast]]))
})




plot_data <- data.frame()

# Extract the results for each contrast
for (contrast in colnames(contrast_matrix)) {
  # Get the topTable results for the contrast
  results <- topTable(fit2, coef = contrast, adjust.method = "BH", number = Inf)
  
  # Select logFC and adj.P.Val, and add them to the plot_data
  plot_data <- as.data.frame(cbind(
    Gene = rownames(expression_data), 
    PG.Genes = df$PG.Genes, 
    logFC = fit2$coefficients[, 1], 
    logAdjP = -log10(fit2$p.value[, 1])
  ))
}

# Check the first few rows of the plot data
head(plot_data)

# Load the necessary libraries
library(ggplot2)

# Create an empty data frame to store logFC and log(adj.P.Val)
plot_data <- data.frame()

# Extract the results for each contrast
for (contrast in colnames(contrast_matrix)) {
  # Get the topTable results for the contrast
  results <- topTable(fit2, coef = contrast, adjust.method = "BH", number = Inf)
  
  # Select logFC and adj.P.Val, and add them to the plot_data
  plot_data <- rbind(plot_data, 
                     data.frame(
                       gene = row.names(results),
                       logFC = results$logFC, 
                       logAdjP = -log10(results$adj.P.Val),  # log10 to avoid log(0)
                       Contrast = contrast
                     ))
}

# Add color column based on conditions
plot_data$Color <- ifelse(plot_data$logAdjP >5 & plot_data$logFC < -1, 'blue', 
                          ifelse(plot_data$logAdjP > 5 & plot_data$logFC > 1, 'red', 'black'))
label_data <- plot_data[plot_data$Color %in% c('blue', 'red'), ]

# Plot with color based on the condition
library(ggplot2)
library(ggrepel)

ggplot(plot_data, aes(x = logFC, y = logAdjP)) +
  geom_point(aes(color = Color)) +
  facet_wrap(~Contrast) +
  scale_color_identity() +  # Use the colors as they are (blue, red, black)
  geom_text_repel(data = label_data, aes(label = gene), 
                  size = 3, box.padding = 0.5, point.padding = 1, 
                  max.overlaps = 50, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "",
       x = "LogFC", y = "-log10(Adj.P.Val)") +
  scale_x_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) +
  theme(legend.position = "none", plot.subtitle = element_text(face = "bold"))+
  facet_wrap(~ Contrast, scales = "free", ncol = 3)  # Create separate plots for each contrast, with free axes and 2 columns
  