# Import libraries
library(hsdar)
library(car)
library(caret)
library(psych)
library(ggridges)
library(ggpubr)
library(reshape2)
library(stringr)

#-------------------------------------------------------------------------------
# Setup working directory
in_fol <- "E:/OneDrive - University of Toronto/EDUCATION/UW-Madison/Projects/Aspen genotype/DrySpec/Leaf_level_dry/inputs"
out_fol <- "E:/OneDrive - University of Toronto/EDUCATION/UW-Madison/Projects/Aspen genotype/DrySpec/Leaf_level_dry/outputs/result_SWIR"


# Load spectra, trait, and wavelength data
ref_fn <- read.csv(paste(in_fol, "/", "WisAsp_2021_LeafSpec_LindrothLab.csv", sep = ""), header = TRUE, sep = ",")
trait_fn <- read.csv(paste(in_fol, "/", "WisAsp_2021_Traits_LindrothLab.csv", sep = ""), header = TRUE, sep = ",")
wavelen_fn <- read.csv(paste(in_fol, "/", "wave_length_res.csv", sep = ""), header = TRUE, sep = ",")


# Setup wavelength and mask
wavelen <- as.numeric(na.omit(wavelen_fn$wl_asd)) # wavelengths of spectra
spec_mask <- c(330, 1000) # mask out unused bands
spec_scale <- 0.8


trait_names = c("Pg_sum", "Salicin", "Salicortin", "Tremuloidin", "Tremulacin", "Ct", "Nitrogen")
for (trait_name in trait_names){
  trait_sub = subset(trait_fn, trait_fn[trait_name] != "NA" & ref_fn[,20] != "NA")
  spec_trait <- subset(ref_fn, trait_fn[trait_name] != "NA" & ref_fn[,20] != "NA")
  
  #-----------------------------------------------------------------------------
  # Create spectral library for each trait
  spec_trait <- spec_trait[, 11:ncol(spec_trait)]
  spec_trait_matrix <- as.matrix(spec_trait)
  ids_trait <- seq(1, nrow(spec_trait), 1)
  speclib_trait <- speclib(spec_trait_matrix, wavelen)
  idSpeclib(speclib_trait) <- as.character(ids_trait)
  

  # Mask out band bands (e.g., noisy bands) from spectra
  mask(speclib_trait) <- spec_mask
  

  # Plot spectra and Save the plot to *.png file
  out_spec <- paste(out_fol, "/", "WisAsp_2021_Spec_", trait_name, ".png", sep = "")
  plot(speclib_trait, FUN=c(1:nrow(spec_trait)), xlim=c(min(wavelength(speclib_trait)), max(wavelength(speclib_trait))), ylim=c(0, spec_scale),
       col=1:nrow(spec_trait), cex.lab=1, cex.axis=1)
  
  png(out_spec, width = 9, height = 11, units = "cm", res = 300)
  plot(speclib_trait, FUN=c(1:nrow(spec_trait)), xlim=c(min(wavelength(speclib_trait)), max(wavelength(speclib_trait))), ylim=c(0, spec_scale),
       col=1:nrow(spec_trait), cex.lab=1, cex.axis=1)
  dev.off()
  
  #-----------------------------------------------------------------------------
  # Calculate available vegetation indices in the hsdar package
  vegIndex_all <- vegindex(speclib_trait, vegindex())
  
  # Correlate available indices (Pearson-r) with plant traits
  vi_trait <- cor(vegIndex_all, trait_sub[trait_name])
  
  # Create a dataframe of the correlation results
  exist_index <- na.omit(as.data.frame(round(vi_trait, 3)), na.rm =TRUE)
  
  # Output Pearson-r values to CSV file
  write.csv(exist_index, paste(out_fol, "/", "Cor_Existing_Indices_", trait_name, ".csv", sep = ""))
  
  
  #-----------------------------------------------------------------------------
  # Calculated all possible normalized vegetation indices
  nri_All <- nri(speclib_trait, recursive = TRUE)
  
  
  # Construct linear model with Pg & Ct content
  lm_f <- as.formula(paste("nri_All ~", trait_name))
  lm_trait <- lm.nri(lm_f, preddata = trait_sub)
  
  # Export the map of R2 to a png file
  png(paste(out_fol, "/", "Index_trait_R2_", trait_name, ".png", sep=""), width = 8, height = 10, units = 'cm', res = 300)
  
  plot(lm_trait, coefficient = "r.squared", range=c(0, 0.8), main=paste("Index-", trait_name, " R2", sep=""),
       cex.lab=1, cex.axis=1, cex=1)
  
  dev.off()
  
  
  # Output n (n=100 in this case) regression models with the highest R2
  best_index <- nri_best_performance(lm_trait, n = 100)
  
  nri_Best <- getNRI(nri_All, best_index) # Get NRI values for best models
  nri_Best_final <- append(trait_sub[trait_name], nri_Best)
  
  # Write these model results to a CSV file
  write.csv(nri_Best_final, paste(out_fol, "/", "R2_nri_best100_", trait_name, ".csv", sep=""), row.names = FALSE)
  
  
  # Print trait name after processing
  print(paste(trait_name, ": completed", sep=""))
  
}




#-------------------------------------------------------------------------------
# Train the regression model for the best index using K-fold cross validation
# Setup plot theme
my_theme = theme(
  strip.text = element_text(size=12),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text.x = element_text(size = 10, colour="black", angle = 0),
  axis.text.y = element_text(size = 10, colour="black"),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major = element_line(colour = "grey60", linetype="dashed", size=0.2),
  legend.title = element_text(size=12),
  legend.text = element_text(size=12, colour="black"),
  legend.position = "top")


# Names of traits to be analyzed
trait_names = c("Pg_sum", "Salicin", "Salicortin", "Tremuloidin", "Tremulacin", "Ct", "Nitrogen")
index_cols = c(2, 2, 2, 2, 2, 2, 2) # Column in the index file (>=2)
model_result_all = NULL

for (trait_id in seq_along(trait_names)){
  best_vegindex <- read.csv(paste(out_fol, "/", "R2_nri_best100_", trait_names[trait_id], ".csv", sep=""), header = TRUE, sep = ",")
  
  # Train the trait models
  set.seed(1)
  k.fold <- trainControl(method = "repeatedcv",
                         repeats = 100,
                         number = 10,
                         verboseIter = TRUE)
  
  lm_f <- as.formula(paste(trait_names[trait_id], "~", colnames(best_vegindex[index_cols[trait_id]])))
  model_trait <- train(lm_f,               # model to fit
                    data = best_vegindex,
                    trControl = k.fold,    # folds
                    method = "lm",         # specifying regression model
                    na.action = na.pass)   # pass missing data to model - some models will handle
  
  # Extract model's parameters
  model_trait
  model_trait$finalModel
  model_trait$bestTune
  model_trait$results
  model_trait$resample
  
  # calculate mean trait values
  mean_trait <- mean(best_vegindex[[trait_names[trait_id]]])
  
  model_result <- data.frame(trait = setNames(data.frame(replicate(nrow(model_trait$resample), trait_names[trait_id])), "Trait"),
                             Rsquared = setNames(model_trait$resample["Rsquared"], "R2"),
                             RMSE = model_trait$resample["RMSE"]/mean_trait
                             )
  
  model_result_all <- rbind(model_result_all, model_result)
  
  
  # Plot linear regression result of selected index
  spe_index <- best_vegindex[[index_cols[trait_id]]]
  trait_val <- best_vegindex[[trait_names[trait_id]]]
  
  png(paste(out_fol, "/", "WisAsp_2021_regression_", trait_names[trait_id], ".png", sep=""),
      width = 11, height = 14, units = "cm", res = 300)
  
  col_name <- colnames(best_vegindex[index_cols[trait_id]])
  
  plot(spe_index, trait_val, xlab=paste("VI (", sub("[.]", ", ", substring(col_name, 3)), " nm", ")", sep=""),
       ylab=paste(trait_names[trait_id], "(mg/g)"), xlim=c(min(spe_index), max(spe_index)),
       ylim = c(min(trait_val), max(trait_val)), pch=21, col="orange", bg="blue", cex=1.3)
  
  abline(lsfit(spe_index, trait_val), col="orange", lwd=2)
  
  trait_range <- max(trait_val) - min(trait_val)
  text(min(spe_index), max(trait_val) - 0.03*trait_range, paste("y = ", round(model_trait$finalModel$coefficients[1],3), " - ",
                                    abs(round(model_trait$finalModel$coefficients[2],3)),
                                    "*", "x",sep=""), pos=4, cex=1.1)
  text(min(spe_index), max(trait_val) - 0.10*trait_range, paste("R2 = ", round(model_trait$results$Rsquared, 3)), pos=4, cex=1.1)
  text(min(spe_index), max(trait_val) - 0.17*trait_range, paste("RMSE = ", round(model_trait$results$RMSE, 3)), pos=4, cex=1.1)
  
  dev.off()
}

model_result_melt <- data.frame(melt(model_result_all, id.vars = "Trait",
                                     variable.name = "metrics",
                                     value.name = c("value")))

# Plot model training results for all traits
model_result_plt <- ggplot(model_result_melt, aes(x=Trait, y=value, fill=Trait)) +
  geom_violin()+
  geom_boxplot(width=0.1, fill="white")+
  facet_wrap( ~ metrics, scales="free_x", labeller=as_labeller(c("R2"="R2", "RMSE"="RMSE (%)")))+
  scale_fill_brewer(palette="Dark2")+
  xlab("Plant traits")+
  ylab("")+
  coord_flip()+
  theme_bw() +
  my_theme
model_result_plt

# Save the plot to a png file
png(paste(out_fol, "/", "WisAsp_2021_model_training.png", sep=""), width = 20, height = 20, units = "cm", res = 300)
model_result_plt
dev.off()

# The End

