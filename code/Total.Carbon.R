library(nlme)
library(ggplot2)
library(emmeans)
library(phyloseq)
library(ggplot2)
library(plyr)
library(tidyverse)
library(xtable)
library(ggpubr)
library(purrr)
library(broom)
library(agricolae)

physeq <- readRDS("Data/incubation_raw.RDS")
inc.data <- subset_samples(physeq, day %in% c("7", "14", "21", "35", "49", "97"))
mydata <- data.frame(sample_data(inc.data))

mydata$day <- as.factor(mydata$day)

mydata <- mydata %>%
  mutate(c_n = C_flash / N_flash)

mydata$treatment <- relevel(mydata$treatment, ref = "Control") #The levels of a factor are re-ordered so that the level specified by ref is first and the others are moved down. This is useful for contr.treatment contrasts which take the first level as the reference.
colnames(mydata)

# Day 7
data <- mydata %>%
  filter(day == 7)

aov <-  aov(C_flash ~ treatment, data = data)
HSD <- HSD.test(aov, "treatment", group = T) 

groups <- data.frame(HSD$groups) %>%
  rownames_to_column("Treatment")
std <- data.frame(HSD$means) %>%
  rownames_to_column("Treatment")

d.7 <- full_join(groups, std) %>%
  mutate(Day = 7) %>%
  select(Treatment, Day, Total_Carbon = C_flash, Standard_Deviation = std, Sig. = groups) 

# Day 14
data <- mydata %>%
  filter(day == 14)

aov <-  aov(C_flash ~ treatment, data = data)
HSD <- HSD.test(aov, "treatment", group = T) 

groups <- data.frame(HSD$groups) %>%
  rownames_to_column("Treatment")
std <- data.frame(HSD$means) %>%
  rownames_to_column("Treatment")

d.14 <- full_join(groups, std) %>%
  mutate(Day = 14) %>%
  select(Treatment, Day, Total_Carbon = C_flash, Standard_Deviation = std, Sig. = groups)

# Day 21
data <- mydata %>%
  filter(day == 21)

aov <-  aov(C_flash ~ treatment, data = data)
HSD <- HSD.test(aov, "treatment", group = T) 

groups <- data.frame(HSD$groups) %>%
  rownames_to_column("Treatment")
std <- data.frame(HSD$means) %>%
  rownames_to_column("Treatment")

d.21 <- full_join(groups, std) %>%
  mutate(Day = 21) %>%
  select(Treatment, Day, Total_Carbon = C_flash, Standard_Deviation = std, Sig. = groups)

# Day 35
data <- mydata %>%
  filter(day == 35)

aov <-  aov(C_flash ~ treatment, data = data)
HSD <- HSD.test(aov, "treatment", group = T) 

groups <- data.frame(HSD$groups) %>%
  rownames_to_column("Treatment")
std <- data.frame(HSD$means) %>%
  rownames_to_column("Treatment")

d.35 <- full_join(groups, std) %>%
  mutate(Day = 35) %>%
  select(Treatment, Day, Total_Carbon = C_flash, Standard_Deviation = std, Sig. = groups)

# Day 49
data <- mydata %>%
  filter(day == 49)

aov <-  aov(C_flash ~ treatment, data = data)
HSD <- HSD.test(aov, "treatment", group = T) 

groups <- data.frame(HSD$groups) %>%
  rownames_to_column("Treatment")
std <- data.frame(HSD$means) %>%
  rownames_to_column("Treatment")

d.49 <- full_join(groups, std) %>%
  mutate(Day = 49) %>%
  select(Treatment, Day, Total_Carbon = C_flash, Standard_Deviation = std, Sig. = groups)

# Day 97
data <- mydata %>%
  filter(day == 97)

aov <-  aov(C_flash ~ treatment, data = data)
HSD <- HSD.test(aov, "treatment", group = T) 

groups <- data.frame(HSD$groups) %>%
  rownames_to_column("Treatment")
std <- data.frame(HSD$means) %>%
  rownames_to_column("Treatment")

d.97 <- full_join(groups, std) %>%
  mutate(Day = 97) %>%
  select(Treatment, Day, Total_Carbon = C_flash, Standard_Deviation = std, Sig. = groups)

vert.table <- rbind(d.7, d.14, d.21, d.35, d.49, d.97) 

vert.table$Total_Carbon <- round(vert.table$Total_Carbon, digits = 2) 
vert.table$Standard_Deviation <- round(vert.table$Standard_Deviation, digits = 2)

test <- paste0("(", format(unlist(vert.table[,4])), ")")
vert.table$Standard_Deviation <- test 


c.table <- unite(data = vert.table, col = Total_Carbon, c("Total_Carbon", "Standard_Deviation", "Sig."), sep = " ")

c.table <- unite(data = vert.table, col = Total_Carbon, c("Total_Carbon", "Standard_Deviation", "Sig."), sep = " ")
