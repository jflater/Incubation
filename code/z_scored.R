# For distance based redundancy analysis, variables must be z-scored
data <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")
meta_data <- data.frame(sample_data(data))
meta_data[7:13] <- data.frame(lapply(meta_data[7:13], function(x) scale(x)))
sample_data(data) <- meta_data  

