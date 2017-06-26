# Regressions with R

library(tidyverse)
library(ggplot2)
library(scales)
library(broom)
library(car)

setwd("C://Users/Ryo Nakagawara/Documents/R materials/RegressionsWithR")
load("Regressions with R datasets.Rdata")

dim(tango)  # 18 rows/observations of 11 variables (including interactions...)
names(tango) # all variables: Genotype, RU, Tango, etc. etc.
tango$RU   # Levels: 'A' or 'P'
class(tango$RU)  # FACTOR!
summary(tango$RU)  # A: 9, P: 9

table(tango$Genotype)
# # of each type of genotypes: Abeta = 6, RNAi = 6, UAS Arm = 6

par(pty = 's')    # shape of graphical window set to SQUARE!
plot(tango$RU, tango$Tango)

# OR in a more detailed graph...

plot(tango[ , 6], tango$Tango,       # column 6 is already in numeric! so use that too....
     ylab = "Tango", xlab = "Treated Group", 
     cex = 2, pch = "X", xaxt = "n", col = "green", cex.axis = 2, cex.lab = 1.5)
axis(1, 0:1, labels = c("Not treated", "Treated"), cex.axis = 2)

# add mean tango expression for each group
points(c(0, 1),
       c(mean(tango$Tango[tango$RU == 'A']),
         mean(tango$Tango[tango$RU == 'P'])), pch = "-", cex = 4)
# From graph: tratement group = avg higher values of Tango expression vs. NOT treatment

tango %>% t.test(Tango ~ RU, data = .)
# RU: treatment present 'P' or NOT present 'A'
# Difference in group means: 0.71-1.01 = ~0.3
# 95% CI = (-0.61, 0.008) with p-value of 0.05586, marginally significant!


m <- lm(Tango ~ RU, data = tango)
summary(m)$r.squared     # r-squared: 0.2106799
summary(m)$adj.r.squared # adjusted: 0.1613474   (taking into account additional parameters)

# fitted regression model: Tango = 0.7095 + 0.3001(Treated: 0 or 1)
# For UNTREATED: 0.71 + (0.3*0) = 0.71
# For TREATED  : 0.71 + (0.3*1) = 1.01


confint(m)   # (-0.008, 0.61)
deviance(m)  # 1.518641
BIC(m)       # 15.2469

fit.m <- fitted(m)   # mean tango expression for treated/non-treated groups
pred.m <- predict(m)
summary(fit.m)
summary(pred.m)

# add again into previous graph:
plot(tango[ , 6], tango$Tango,      
     ylab = "Tango", xlab = "Treated Group", 
     cex = 2, pch = "X", xaxt = "n", col = "green", cex.axis = 2, cex.lab = 1.5)

axis(1, 0:1, labels = c("Not treated", "Treated"), cex.axis = 2)

points(tango[ , 6], fit.m, pch = "--", cex = 4)

# Predictions for specific value of independent variable
pred.m1 <- predict(m, data.frame(RU = 'P'))
pred.m1 # 1.009606

# Residuals:




































