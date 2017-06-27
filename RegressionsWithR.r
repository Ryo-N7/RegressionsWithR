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

tango %>% ggplot(aes(RU, Tango)) + 
                 geom_point(aes(color = RU, shape = RU), size = 4) +
                 theme(axis.text.x = element_text(face = "bold", size = 15)) +
                 scale_shape_discrete(name = "Treatment \nCondition", labels = c("No Treatment", "Treatment")) +
                 scale_color_discrete(name = "Treatment \nCondition", labels = c("No Treatment", "Treatment")) +
                 scale_x_discrete(name = "Treatment Condition", labels = c("No Treatment", "Treatment"))
# LEGEND IS REDUNDANT SO CAN JUST DELETE AS WELL WITH: legend.position = "none" in theme()


tango %>% t.test(Tango ~ RU, data = .)
# RU: treatment present 'P' or NOT present 'A'
# Difference in group means: 0.71-1.01 = ~0.3
# 95% CI = (-0.61, 0.008) with p-value of 0.05586, marginally significant!


m <- lm(Tango ~ RU, data = tango)
names(summary(m))
summary(m)$r.squared     # r-squared: 0.2106799
summary(m)$adj.r.squared # adjusted: 0.1613474   (taking into account additional parameters)
summary(m)$fstatistic    # f = 4.27 df 1,16
summary(m)$residuals
summary(m)$coefficients

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

# OR using scatterplot

plot(as.numeric(tango$RU), tango$Tango)
points(as.numeric(tango$RU), fitted(m), col = 'red', 
                  pch = '-', cex = 3)

# Predictions for specific value of independent variable
pred.m1 <- predict(m, data.frame(RU = 'P'))
pred.m1 # 1.009606

# Residuals:
qqnorm(residuals(m))
qqline(residuals(m), col = 6)

plot(predict(m), residuals(m))
abline(h = 0, col = "blue")

# Influence measures:
inf.m <- influence.measures(m)
par(mfrow = c(3,3))
for (i in 1:(dim(inf.m$infmat) [2]-1))

# add ID variable for ease of identification
tangoID <- 1:dim(tango)[1]
tango2 <- data.frame(tangoID, tango)
tango2$tangoID[cooks.distance(m) > 0.2]   # 2!

# BIC more conservative and takes into account sample size + penalizes for ^# of parameters
# Can use for comparison of non-nested models (lower BIC = better fit!)
BIC(m)  # 15.2469..........     USELESS BY ITSELF
AIC(m)  # 12.57578.........     NEED OTHER MODELS FOR COMPARISON
deviance(m) # 1.518641
cooks <- cooks.distance(m)
which(cooks > 1)   # 0!

layout(matrix(c(1,2,3,4), 2,2))
plot(m)
# obsv 2, 14, 10 seem like outliers necessary for further inspections!!
tango[2 ,]
tango[14, ]
tango[10, ]

glm(Tango ~ RU + Genotype, data = tango, family = "gaussian")
tango %>% do(tidy(glm(Tango ~ RU + Genotype, data = .)))
# RUP and RNAi = significant, UAS Arm = NOT signif. 

plot(tango$Genotype, tango$Tango)

m2 <- lm(Tango ~ RU:Genotype, data = tango)
summary(m2)
# RU = A|RNAi significant
# RU = P|RNAi significant
tango %>% ggplot(aes(Genotype, Tango)) + geom_point()

# Exercises #1:
m.tango2 <- glm(Tango~RU, data = tango)
summary(m.tango2) # intercept: 0.70, slope = 0.3    SAME AS lm()
# null deviance: deviance of model with no parameters (the null model)
# residual deviance: 
# AIC: 

plot(as.numeric(tango$Genotype), tango$Tango)
plot(tango$Genotype, tango$Tango)

m.tango.geno <- lm(Tango ~ Genotype, data = tango)
summary(m.tango.geno)   # ABeta as reference.
# on avg. tango expression is lower for RNAi comparison to ABeta
# on avg. tango expression is higher for UAS Arm comparison to ABeta
# adj.R = 0.352 (improve from 0.16 with treatment variable)

plot(as.numeric(tango$Genotype), tango$Tango,
     xaxt = "n",
     col = c("red", "green")[tango$RU],
     pch = c("X", "O")[tango$RU])
axis(1, 1:3, labels = levels(tango$Genotype), cex.axis = 1)

tango %>% ggplot(aes(Genotype, Tango)) + 
  geom_point(aes(color = RU, shape = RU), size = 4) +
  theme(axis.text.x = element_text(face = "bold", size = 15)) +
  scale_x_discrete(name = "Genotype", labels = c("UAS Arm", "ABeta", "RNAi")) +
  scale_shape_discrete(name = "Treatment \nCondition", labels = c("No Treatment", "Treatment")) +
  scale_color_discrete(name = "Treatment \nCondition", labels = c("No Treatment", "Treatment"))


m.tango.gen.RU <- lm(Tango ~ Genotype + RU, data = tango) 
summary(m.tango.gen.RU)
# Genotype RNAi = significant, lower comparison to ABeta.
# Adj.R = 0.5615 (improve from previous models...!)
# p = 0.002079

m.tango.int <- lm(Tango ~ Genotype*RU, data = tango)
summary(m.tango.int)
# NOT significant.
# Interactions are not necessary for these predictions.
# Adj.R = 0.49......
# p = 0.01706

# Change reference with relevel(), redefine in Genotype variable in dataset:
tango$Genotype <- relevel(tango$Genotype, ref = "UAS Arm")
m.tango.int <- lm(Tango ~ Genotype*RU, data = tango)
summary(m.tango.int) # from re-running, ABeta appears with UAS Arm as reference instead.



# LOGISTIC REGRESSION #

titanic
View(titanic)
str(titanic)
table(titanic$survived)   # 0 = NOT Survive, 1 = Survived

plot(titanic$survived)   # not very useful.....
summary(titanic$age)     # avg age: 29.8, NA = 263..., Median = 28

plot(jitter(titanic$age), jitter(titanic$survived))  # add noise to see more clusters...

class(titanic$survived)   # as integer!
titanic$survived <- factor(titanic$survived,
                           levels = c(0, 1),
                           labels = c("Not_Survive", "Survive"))
class(titanic$survived)   # now as factor!

table(titanic$survived)
plot(titanic$age, titanic$survived)

prop.table(table(titanic$survived))*100      # NOT = 61.8%, Survive = 38.1%
tapply(titanic$age, titanic$survived, mean, na.rm = TRUE) 
# NOT = 30.54 AGE, Survive = 28.9 AGE

m.age <- glm(survived ~ age, data = titanic,
             family = 'binomial')
summary(m.age)
# Coefficients: 
# Age: p-value = 0.07, NOT significant
# Z-vaule: distance of estimate from 0 in terms of s.e.

exp(coef(m.age))        # intercept: 0.87, age = 0.99 
# age: any unit change in age +1, odds of survival are the same (1)
# odds ratio for survival for ex. age = 40 is the same for age = 41, etc.

exp(confint(m.age))     # intercept: (0.656, 1.15), age: (0.983, 1.0007)

# Practice #2:
# explore variables: sex, age, fare, passenger class
table(titanic$sex) # female = 466, male = 843
hist(as.numeric(titanic$sex))    # fix the x-axis labels but yeah...

class(titanic$sex)
str(titanic$sex)
titanic$sex <- factor(titanic$sex,
                      levels = c(0, 1),
                      labels = c("Female", "Male"))
titanic$sex
table(titanic$sex)
prop.table(table(titanic$sex))

plot(titanic$sex, titanic$survived) # Female: 35.6%, Male: 64.4% 


titanic %>% ggplot(aes(x = sex, y = survived)) + 
                  geom_histogram(stat = "identity")
# Female survivors > Male survivors!

prop.table(table(titanic$sex, titanic$survived), 1)*100 # proportions by row? or columns ',2'

with(titanic,
     table(sex, survived, pclass))

pairs(~ survived + age + fare, data = titanic)

m.AS <- glm(survived ~ age + sex, data = titanic, family = 'binomial')
summary(m.AS)
# Sex(Male): -2.46, p = 2e^-16, signif., Age = NOT significant.
exp(coef(m.AS)) # age = 0.99 odds, sex(male) = 0.08 odds
# odds of survival if male is VERY LOW 0.08 compared to female!
# (1-0.08)*100 = 92%, Odds of survival for male is 92% LOWER than Female! :O
exp(confint(m.AS))

class(titanic$pclass)   # integer...
m1 <- glm(survived ~ age + sex + factor(pclass) + fare, data = titanic,
          family = 'binomial')
summary(m1)    # ALL variables significant except 'fare'
exp(coef(m1)) # age = 0.966 odds, sex(male) = 0.08 odds, pclass = 2: 0.28 odds,
              # pclass = 3: 0.104 odds, fare: 1.000 odds

# pclass = 2: -1.257, pclass = 3: -2.261
# age: -0.03, sex(male): -2.493
# (1-0.28)*100 = 72%, Odds of survival for passenger class = 2 is 72% lower than pclass = 1.
# (1-0.104)*100 = 90%, Odds of survival for pclass = 3 is 90% lower than pclass = 1!

m2 <- glm( survived ~ age + sex 
           + factor(pclass),
           data=titanic, 
           family='binomial')
summary(m2)


m5 <- glm(survived ~ age * factor(pclass), data = titanic,
          family = 'binomial')

m6 <- glm(survived ~ sex * factor(pclass), data = titanic,
          family = 'binomial')

m7 <- glm(survived ~ sex + age + factor(pclass) + sex * age + sex * factor(pclass),
          data = titanic,
          family = 'binomial')

exp(m7$coe) # sexmale*pclass = 3: compared to BASE female pclass = 1, odds of survival
            # for male pclass = 3 is LOWER

library(binomTools)


m.rsq1 <- Rsq(m1)
m.rsq2 <- Rsq(m2)


## Ordinal Logistic Regression ##
View(cd)
cd$height <- factor(cd$height, levels = c(1, 2, 3),
                    labels = c("bottom", "middle", "top"))
cd$height
table(cd$height)  # bottom: 764 flies, middle: 505 flies, top: 752 flies reach
table(cd$genomic) # foxonull: 992 flies, wt: 1029 flies
table(cd$genomic, cd$height)

addmargins(prop.table(table(cd$genomic, cd$height), 1) * 100)
# foxonull genotype = LESS likely to climb high compared to wt genotype
library(ordinal)

m <- clm(height ~ genomic, data = cd)

summary(m)
# genomic = wt: 1.52
# odds of wt being in higher category = exp(1.52) = 4.57 (CI: 3.86, 5.43) times higher than for foxonull (p < 2e^-16)
exp(m$coef)  # genomic = wt: 4.57 (as above) to climb higher compared to foxonull
exp(confint(m)) # CI: (3.85, 5.43)
# Threshold: 
# Bottom: 552/(234+206) == 1.25, similar to bottom|middle = 1.234

predict(m, type = 'prob')
predict(m, type = 'prob', data.frame(genomic = 'wt')) # predict for genomic = wt
# bottom: 0.21    middle: 0.25    top: 0.534
predict(m, type = 'prob', data.frame(genomic = 'foxonull')) # predict for genomic = foxonull
# bottom: 0.552   middle: 0.246   top: 0.200

nominal_test(m)
# loglikelihood.     insignificant. FAIL TO REJECT NULL. proportionality holds.

# Exercise #3:
table(cd$RU)   # treatment: 1033, NO treatment: 988
table(cd$day)

addmargins(prop.table(table(cd$RU, cd$height), 1) * 100)
# with Treatment more likely to climb to the top!

m2 <- clm(height ~ RU, data = cd)
summary(m2)

exp(m2$coef[3]) # odds of fly with Treatment being in higher category is 1.62 times higher than NO treatment

m3 <- clm(height ~ factor(day), data = cd)
summary(m3)
exp(m3$coef[3])   # odds of fly with ^ day being in higher category is 0.95 times higher than less days

