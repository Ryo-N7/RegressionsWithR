###########
# No need to run the following two commands 
# but it's good to have here for quick
# reference if you want to clear the memory
# or list all defined objects

ls()
rm(list=ls())

###########
# check the working directory
getwd()

# change/set it to where you have saved the data
# for example
setwd("K:\\StatsAdmin\\Lecture Notes-Presentations\\Regressions with R\\datasets\\Final datasets")
setwd("S:\\ICH\\StatsAdmin\\Lecture Notes-Presentations\\Regressions with R\\datasets\\Final datasets")
getwd()

###########
# Load all the datasets in R by just dragging and dropping the .Rdata file
# on the R console or by just running the following line
load("Regressions with R datasets.Rdata")
# Alternatively, you can read in one dataset at a time using the 
# read.csv functions seen below

###########
# Install required packages - these commands 
# are here for quick reference, but each will
# be repeated at the point where it's needed

install.packages("car")
install.packages("binomTools")
install.packages("ordinal")
install.packages("survival")
install.packages("nlme")
install.packages('rms')

require(MASS) #pre-installed
require(car)
require(binomTools)
require(ordinal)
require(survival)
require(nlme)
require(rms)

### EXERCISE 1 ###########
##########################
########### TANGO ANALYSIS
##########################
# Read in the dataset
tango <- read.csv("QPCR tango.csv")
names(tango)
dim(tango)

# Linear Regression with intercept only 	
# default choice of family 'gaussian'
m <- glm(Tango ~ 1 , data=tango)
summary(m)
#### Interpretation
Notice that the value of the intercept matches
the overall mean values of the tango variable
mean(tango$Tango)
####

# Linear Regression with ONE binary predictor 	
m <- glm(Tango ~ RU , data=tango)
summary(m)
#### Interpretation
Repeating the analysis using the glm as opposed
to the lm function leads to the same results 
as shown on the notes. The glm function does not
allow the calculation of the R squared.
#### 

# To investigate the relationship between Tango
# and Genotype, we should start with some
# descriptives statistics:
# summaries by genecode 
table(tango$Genotype)
# table(tango$genecode)
tapply (tango$Tango,tango$Genotype,summary)
tapply (tango$Tango,tango$Genotype,sd)
# Similarly, if we want to look at Genotype
# and treatment against the tango expression
# at the same time
tapply (tango$Tango,
	list(tango$Genotype,tango$RU),mean)
tapply (tango$Tango,
	list(tango$Genotype,tango$RU),sd)

# Dot plot
trF <- table(tango$RU)
trF
par(pty="s")
plot(as.numeric(tango$Genotype),tango$Tango,
	ylab="Tango",xlab="Genecode",
	cex=2,xaxt="n",
	col=c("red","green")[tango$RU],
	cex.axis=2,
	cex.lab=1.5,
	pch=c("X","O")[tango$RU])
axis(1,1:3,labels=levels(tango$Genotype),
	cex.axis=2)
legend("top",legend=levels(tango$RU)[c(1,2)],
	col=c("red","green")[c(1,2)],
	pch=c("X","O")[c(1,2)],box.lty=2)

# Linear regression with ONE nominal predictor
m <- glm(Tango ~ Genotype, data=tango)
summary(m)
#### Interpretation
The average tango values for the Abeta gene
is equal to the value of the intercept, 0.9486.
The tango values for the RNAi gene are on average
lower by 0.38 units compared to Abeta
The UAS_Arm gene has average tango values higher
by 0.11 units compared to Abeta
####

# Try changing the reference category of the
# categorical predictor
Genotype2 <- relevel(tango$Genotype,ref='RNAi')
m <- glm(Tango ~ Genotype2, 
	data=tango)
summary(m)
# or 
Genotype3 <- relevel(tango$Genotype,ref='UAS Arm')
m <- glm(Tango ~ Genotype3, data=tango)
summary(m)

# Linear Regression with two predictors 	
m <- glm(Tango ~ RU + Genotype,
	data=tango)
summary(m)
#### Interpretation
The average tango value for the Abeta gene 
when treatment is absent is equal to the 
intercept, 0.7985.
The tango values for the RNAi gene are on average
lower by 0.38 units compared to Abeta with absent 
treatment
The UAS_Arm gene has average tango values higher
by 0.11 units compared to Abeta with absent 
treatment
The average tango expression is 0.3 higher 
for when the treatment was present compared to absent 
regardless of the gene variable
Adjusting for the treatment/genotype has not 
caused much of a change on the coefficients 
of genotype/treatment
####

# Linear regression with interactions
m <- glm(Tango ~ RU + Genotype + 
	RU * Genotype, data=tango)
summary(m)
#### Interpretation
The interaction term is not statistically 
significant meaning that there is no evidence
to suggest that treatment affects the tango 
expression differently for different groups of
gene
####

# Superimpose regression 'line' on data
par(pty="s")
plot(as.numeric(tango$Genotype),tango$Tango,
	ylab="Tango",xlab="Genecode",
	cex=2,xaxt="n",
	col=c("red","green")[tango$RU],
	cex.axis=2,
	cex.lab=1.5,
	pch=c("X","O")[tango$RU])
axis(1,1:3,labels=levels(tango$Genotype),
	cex.axis=2)
legend("top",legend=levels(tango$RU),
	col=c("red","green")[1:2],
	pch=c("X","O")[1:2])

points( tango$Genotype ,fitted(m),
	col=c("red","green")[tango$RU],
	cex=3,pch="-")
	
# Comparison of models can be achieved by 
# comparing their BICs and/or deviances
# or even the adjusted R squared values if the
# models have been fitted with the lm function
# and not glm	
m1 <- lm(Tango ~ 1,	data=tango)
BIC(m1)
deviance(m1)
summary(m1)$adj.r.squared
m2 <- lm(Tango ~ RU, data=tango)
BIC(m2)
deviance(m2)
summary(m2)$adj.r.squared
m3 <- lm(Tango ~ Genotype,data=tango)
BIC(m3)
deviance(m3)
summary(m3)$adj.r.squared
m4 <- lm(Tango ~ RU+Genotype,data=tango)
BIC(m4)
deviance(m4)
summary(m4)$adj.r.squared
m5 <- lm(Tango ~ RU*Genotype2,data=tango)
BIC(m5)
deviance(m5)
summary(m5)$adj.r.squared

##########################
# Developmental assay data
##########################
# Read in the dataset
dev <- read.csv("dev assay MR original.csv")
dim(dev)
names(dev)

table(dev$Trial)
# split the dataset by trial
# Trial A
devA <- dev[dev$Trial=="A",]
dim(devA)
# Trial B
# not needed for practical
devB <- dev[dev$Trial=="B",]
dim(devB)
# Trial C
# not needed for practical
devC <- dev[dev$Trial=="C",]
dim(devC)

par(pty="s")
plot(devA$Day,devA$Length,
	xlab="Days",
	ylab="Length")

# add a line on the existing plot connecting 
# the actual mean lengths at each of the days

lines(sort(unique(devA$Day)),
	tapply(devA$Length,devA$Day,mean),
	col='red',pch='X',lwd=2)
	
# depending on which trial it is that you want 
# to do analysis for, define a new object
# dev_sub (as in dev sub-dataset) to be equal 
# to the specifc sub-dataset
dev_sub <- devA
dev_sub <- devB
dev_sub <- devC
# the following command sets sub_dev to the original
# big dataset
dev_sub <- dev
# for each of them define the day squared
Day2 <- (dev_sub$Day)^2

# Trial A
# Linear regression with numeric predictor, Day
m1 <- glm(Length ~ Day, data=dev_sub)
# m1 <- lm(Length ~ Day, data=dev_sub)
# summary(m1)$adj.r.squared
summary(m1)	
plot(dev_sub$Day,dev_sub$Length,
	xlab="Days",
	ylab="Length")
lines(dev_sub$Day,fitted.values(m1),lwd=5)
# or 
points(unique(dev_sub$Day),
	unique(fitted.values(m1)),
	col="red",
	type="l",cex=2,lwd=5)

# Linear regression with two numeric predictors
# Recall we have defined Day2 as: 
# Day2 <- (dev_sub$Day)^2
m2 <- glm(Length ~ Day + Day2, data=dev_sub)
# m2 <- lm(Length ~ Day + Day2, data=dev_sub)
# summary(m2)$adj.r.squared
summary(m2)	
# plot(dev_sub$Day,dev_sub$Length,
#	xlab="Days",
#	ylab="Length",
# 	ylim=c(0,2000))
lines(unique(dev_sub$Day),unique(fitted(m2)),
	col="blue",
	cex=2,lwd=5)

plot(fitted(m2),residuals(m2))

# Linear regression with Group 
with (dev_sub, plot(as.numeric(Group),Length))
m3 <- glm(Length ~ Group, data=dev_sub)
summary(m3)
	
# Linear regression with Group and Day
with (dev_sub, plot(Day,Length,col=c(1,2,3,4,5,6,7,8)[Group]))
legend('topleft',legend=unique(dev_sub$Group),
	col=1:8,pch=rep(1,8),box.lty=2,cex=0.7)
m4 <- glm(Length ~ Day + Group, data=dev_sub)
summary(m4)	

# Linear regression with Group and Day
m5 <- glm(Length ~ Day + Day2 + Group, data=dev_sub)
summary(m5)
	
# Adding interactions with Day	
m6 <- glm(Length ~ Day + Day2 + Group 
	+ Day*Group, data=dev_sub)
summary(m6)	

# Adding interactions with Day squared 	
m7 <- glm(Length ~ Day + Day2 + Group 
	+ Day2*Group, data=dev_sub)
summary(m7)	

# Adding interactions with both Day and Day squared	
m8 <- glm(Length ~ Day + Day2 + Group 
	+ Day*Group + Day2*Group, data=dev_sub)
summary(m8)	

plot(dev_sub$Day,dev_sub$Length,ylim=c(0,1750),
		col=c(1,2,3,4,5,6,7,8)[dev_sub$Group])
for (i in 1:8){
lines(unique(dev_sub$Day[dev_sub$Group==levels(dev_sub$Group)[i]]),
	unique(predict(m8)[dev_sub$Group==levels(dev_sub$Group)[i]]),
	col=i,
	type="l",cex=3,
	lwd=2
	)
}
legend("topleft",legend=unique(dev_sub$Group),
	col=1:8,lty=rep(1,8),box.lty=2,cex=0.7)

#### Interpretation
Day seems to be consistently significantly 
for all models hence it would be wise to keep
it in the final model. The same goes for the
squared Day variable. The Group variable 
seems to come out as statistically significant
when accounted with Day and Day2 but not 
when the interactions are added too.

Looking at model 5
(Length ~ Day + Day2 + Group), we conclude
that for every extra Day the length value 
increases on average by 400 units, but this 
increase slows down by 27 units the higher the
days get whilst keeping the Group variable 
constant. Each of the Group coefficients shows 
the change in length compared to the baseline 
Group DR1567 whist keeping the Days constant. 
####

	
### EXERCISE 2 ###########
##########################
################## TITANIC
##########################
# Read in the dataset
titanic <- read.csv('titanic_final.csv')
dim(titanic)
names(titanic)
head(titanic)
 
table(titanic$survived)
addmargins(table(titanic$survived))
prop.table(table(titanic$survived))*100
prop.table(table(titanic$pclass))*100
prop.table(table(titanic$sex))*100
summary(titanic$fare)
table(titanic$sex,titanic$survived)
prop.table(table(titanic$sex,titanic$survived),1)*100

pairs(~ survived + age + fare,
	data=titanic)
	
plot (titanic$age, jitter(titanic$survived),
	ylab="Survived",xlab="Age",
	cex=1,yaxt="n",
	cex.axis=2,
	cex.lab=1.5,
	col=c(1,2)[titanic$sex],
	pty='s')
axis(2,0:1,labels=c('no','yes'),
	cex.axis=2)
legend('center',legend=c('female','male'),col=c(1,2),
	pch=c(1,1))

#### Interpretation
More mare were in the non survival group
####
	
m1 <- glm( survived ~ age + sex 
			+ fare + 
			factor(pclass),
		data=titanic, 
		family='binomial')
summary(m1)
#### Interpretation
All predictors are statistically significant, 
apart from fare.
####

m2 <- glm( survived ~ age + sex 
			+ factor(pclass),
		data=titanic, 
		family='binomial')
summary(m2)

exp(m2$coef)
exp(confint(m2))
#### Interpretation
The odds of survival for a male are 0.08 those
of a female after accounting for age,fare+class.
In other words, the odds of survival for males 
are lower by (1-0.08)*100=92% compared to 
females. Age was statistically significant but
its effect was not large, the OR for every
unit change in age was approximately equal to
1. Also, the ORs for the two passenges class 
levels are lower than 1, indicating that the both
classes 2 and 3 had lower odds of survival 
compared to class 1.
####

install.packages('binomTools')
require(binomTools)
m.rsq1 <- Rsq (m1)
HLtest(m.rsq1)
m.rsq2 <- Rsq (m2)
HLtest(m.rsq2)

# predicted probabilities
predict(m2,type='response')
hist(predict(m2,type='response'))
predict(m2,type='response',
	data.frame(age=20,sex='female',fare=33,pclass=2))
predict(m2,type='response',
	data.frame(age=50,sex='female',fare=33,pclass=2))

m3 <- glm( survived ~ age * sex * factor(pclass),
	data=titanic, family='binomial')	
summary(m3)

m4 <- glm( survived ~ age * sex ,
	data=titanic, family='binomial')
summary(m4)

m5 <- glm( survived ~ age * factor(pclass) ,
	data=titanic, family='binomial')
summary(m5)

m6 <- glm( survived ~ sex * factor(pclass) ,
	data=titanic, family='binomial')
summary(m6)

m7 <- glm( survived ~ sex + age + factor(pclass) +
			sex*age + sex*factor(pclass),
	data=titanic, family='binomial')
summary(m7)
exp(m7$coef)
exp(confint(m7))
m.rsq7 <- Rsq (m7)
HLtest(m.rsq7) # non significant result, which
	# indicates a good fit
#### Interpretation
The main effect of the three predictors has not
changed much after including the significant 
interactions in model 7. 
####

# Understanding the intercept
predict(m7,type='link',
	data.frame(age=0,
	sex='female',pclass=1))
	# this is the value of the raw intercept 
	# of m7
predict(m7,type='response',
	data.frame(age=0,
	sex='female',pclass=1))
	# this is the probability of survival for 
	# a female of age 0 and pclass 1

	
### EXERCISE 3 ###########
##########################
################# CLIMBING
##########################
# Read in the dataset
cd <- read.csv("climbing.csv")
names(cd)
dim(cd)	
cd$height <- factor(cd$height,levels=c(1,2,3),
	labels=c('bottom','middle','top'))

install.packages("ordinal")
require(ordinal)

m1 <- clm (height ~ RU ,data=cd)
summary(m1)
exp(m1$coef[3])
exp(confint.default(m1))[3,]

m2 <- clm (height ~ day ,data=cd)
summary(m2)
exp(m2$coef[3])
exp(confint.default(m2))[3,]

m3 <- clm (height ~ genomic + RU + day,data=cd)
summary(m3)
exp(m3$coef[3:5])
exp(confint.default(m3))[3:5,]

m4 <- clm (height ~ genomic*RU*day ,
	data=cd,x=TRUE,y=TRUE)
summary(m4)

m5 <- clm (height ~ genomic*RU + day ,
	data=cd,x=TRUE,y=TRUE)
summary(m5)
	# no significant interactions

#### Interpretation
Final model is m3 (height~genomic+RU+day)
	summary(m3)
	exp(m3$coef[3:5])
	exp(confint.default(m3))[3:5,]
The odds of wt gene to climb higher are 6.12 
times to odds of foxonull whilst day and 
treatment are kept constant.
####
	
# for diagnostics
require(rms)
m3.2 <- lrm (height ~ genomic+day+RU ,
	data=cd,x=TRUE,y=TRUE)
m3.2
vif(m3.2) 
	# only relevant to models with more 
	# than 1 predictors
hist(residuals(m3.2))
hist(resid(m3.2,type='dfbetas'))
hist(resid(m3.2,type='dffit'))


### EXERCISE 4 #######
##########################
################## FEEDING
##########################
# Read in the dataset
fa <- read.csv("feeding assay.csv")
dim(fa)
names(fa)

table(fa$Genotype)
tapply(fa$Sumfeeding,fa$Genotype,mean)
tapply(fa$Sumfeeding,fa$Genotype,var)
tapply(fa$Sumfeeding,fa$Genotype,sd)
tapply(fa$Sumfeeding,
	list(c(fa$Food),c(fa$Genotype)),mean)
tapply(fa$Sumfeeding,
	list(c(fa$Food),c(fa$Genotype)),var)
	
# from the notes
m.pois <- glm(Sumfeeding ~ Food,
	data=fa,family=poisson)
summary(m.pois)
BIC(m.pois)

# GLM Poisson model with foodfactor + genotype
m1 <- glm(Sumfeeding ~ Food + Genotype,
	data=fa,family=poisson)
summary(m1)
exp(m1$coef)
exp(confint(m1))
deviance(m1)
BIC(m1)
	# All predictors are significant and 
	# the BIC has dropped
	
predict(m1,data.frame(Genotype='G69',Food='1Y'),type='response')
predict(m1,data.frame(Genotype='G69',Food='2Y'),type='response')
predict(m1,data.frame(Genotype='G69',Food='0.2Y'),type='response')
predict(m1,data.frame(Genotype='W',Food='1Y'),type='response')
predict(m1,data.frame(Genotype='W',Food='2Y'),type='response')
predict(m1,data.frame(Genotype='W',Food='0.2Y'),type='response')

# GLM Poisson model with interaction
m2 <- glm(Sumfeeding ~ Genotype*Food,
	data=fa,family=poisson)
summary(m2)
BIC(m2)
		# BIC is even lower compared to before
par(pty="s")
plot(as.numeric(fa$Food),fa$Sumfeeding,
	ylab="",xlab="Food type",
	cex=1,xaxt="n",col=c("blue","red")[fa$Genotype],
	cex.axis=2,
	cex.lab=1.5,
	pch= c("X","O")[fa$Genotype])
axis(1,1:3,labels=levels(fa$Food),
	cex.axis=2)
legend("bottom",legend=levels(fa$Genotype)[1:2],
	col=c("blue","red")[1:2],
	pch=c("X","O")[1:2],
	box.lty=2)

points(1:3, predict(m,type='response',
	data.frame('Food'=c('0.2Y','1Y','2Y'),
				'Genotype'='W')),
	pch='-',cex=4, col='red')

points(1:3, predict(m,type='response',
	data.frame('Food'=c('0.2Y','1Y','2Y'),
				'Genotype'='G69')),
	pch='-',cex=4, col='blue')

	
### EXERCISE 5 ###########
##########################
################# LIFESPAN
##########################
lfs <- read.csv("lifespan.csv")
names(lfs)
dim(lfs)	
require(survival)	
lfsO <- Surv(lfs$day,lfs$event,type="right")	

table(lfs$trial)
table(lfs$RU)	
prop.table(table(lfs$trial,lfs$event),2)
prop.table(table(lfs$RU,lfs$event),2)

m1 <- coxph (lfsO ~ trial,	data =lfs)	
summary(m1)
exp(m1$coef)
exp(confint(m1))
	# There is a significant difference of
	# lifespan between the two trials.
	# Trial B has lower hazard ratio by 
	# (1-0.63)*100=37%

m2 <- coxph (lfsO ~ RU, data =lfs)	
summary(m2)
exp(m2$coef)
	# Lower HR for present treatment by
	# (1-0.77)*100=23%, significant
	
m3 <- coxph (lfsO ~ genotype + trial, data =lfs)	
summary(m3)
exp(m3$coef)
	# Accounting for both gene and trial, 
	# both variables remain significant  

m4 <- coxph (lfsO ~	genotype + RU, data =lfs)	
summary(m4)

m5 <- coxph (lfsO ~	genotype + RU + trial,
	data =lfs)	
summary(m5)

m6 <- coxph (lfsO ~	genotype * RU * trial,
	data =lfs)	
summary(m6)

require(rms)
m6.2 <- cph(lfsO ~	genotype * RU * trial,
	data =lfs,x=TRUE,y=TRUE)	
vif(m6.2)


### EXERCISE 6 ###########
##########################
################ DEV-ASSAY
##########################

# same as on the notes





