install.packages('plot3D')

#Load libraries
install.packages("copula")
install.packages('gamlss')
install.packages('TTR')
install.packages('quantmod')
install.packages('fitdistrplus')
install.packages('VineCopula')
install.packages("VC2copula")
# install.packages('scatterplot3d')

library('plot3D')
library("copula")
library('gamlss')
library('TTR')
library('quantmod')
library('fitdistrplus')
library('VineCopula')
# library('scatterplot3d')
library('plot3D')
library('ggplot2')
library('ggExtra')
library('VC2copula')




############## Data Import and Preprocessing ##############


raw_data <- read.csv('../input/covid-dataset/owid-covid-data.csv')
#Remove continents from location
raw_data <- raw_data[raw_data$continent != "",]

#Retain Important Columns
raw_data_imp_cols <- raw_data[,c(3,11,14,51,56,57,61,63)]

#Deal with NAs (Since all data columns are 0 or above, replacing NAs with 0 will not affect the max function)
raw_data_imp_cols[is.na(raw_data_imp_cols)] <- 0

#Aggregate Data
agg_data <- aggregate(x=raw_data_imp_cols, by= list(raw_data_imp_cols$location), FUN=max)

#Removing group 1 col
agg_data <- agg_data[,c(2,3,4,5,6,7,8,9)]

#Adding Mortality Rate
agg_data$Mortality <- agg_data[,3]/agg_data[,2]

#Check for missing data
#Remove countries where 5 columns are empty
agg_data_rmissing_data <- agg_data[agg_data$median_age != 0 & agg_data$cardiovasc_death_rate != 0 & agg_data$diabetes_prevalence != 0 & agg_data$hospital_beds_per_thousand != 0 & agg_data$human_development_index != 0,]

#Remove countries for which mortality rate is not available
agg_data_rmissing_data <- agg_data_rmissing_data[agg_data_rmissing_data$Mortality != 'NaN',]
head(agg_data_rmissing_data)



############## Original Data Dsitributions ##############

hist(agg_data_rmissing_data$median_age)
plot(agg_data_rmissing_data$median_age, agg_data_rmissing_data$Mortality)
plot(agg_data_rmissing_data$median_age, agg_data_rmissing_data$total_cases_per_million)
plot(agg_data_rmissing_data$median_age, agg_data_rmissing_data$total_deaths_per_million)
plot(agg_data_rmissing_data$human_development_index, agg_data_rmissing_data$total_deaths_per_million)






############## Identify Distributions ##############

identify_distr <- function(dist_col)
{
  
  fit <- fitDist(dist_col, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)
  
  summary(fit)
  fit
  
}

#Mortality
fit <- identify_distr(agg_data_rmissing_data$total_cases_per_million)



histDist(agg_data_rmissing_data$total_cases_per_million,"BCPEo",nbins = 100)
plotdist(agg_data_rmissing_data$total_cases_per_million,"BCPEo",para = list(11.3241588,0.3717137,0.3218911,3.2268296))



pobs_Mortality <- pobs(agg_data_rmissing_data$total_cases_per_million)



#Median Age
fit <- identify_distr(agg_data_rmissing_data$median_age)



histDist(agg_data_rmissing_data$median_age,"BCPE")
plotdist(agg_data_rmissing_data$median_age,"BCPE",para = list(30.9,0.283,0.816,11.5))



pobs_medianage <- pobs(agg_data_rmissing_data$median_age)



#Human Development Index (Can't Identify distribution properly)
fit <- identify_distr(agg_data_rmissing_data$human_development_index)



histDist(agg_data_rmissing_data$human_development_index,"BCPEo",nbins = 10)
#plotdist(agg_data_rmissing_data$human_development_index,"BCPE",para = list(30.9,0.283,0.816,11.5))



pobs_HDI <- pobs(agg_data_rmissing_data$human_development_index)



#Diabetes prevelance
fit <- identify_distr(agg_data_rmissing_data$diabetes_prevalence)



histDist(agg_data_rmissing_data$diabetes_prevalence,"GA",nbins = 10)
plotdist(agg_data_rmissing_data$diabetes_prevalence,"GA",para = list(8.07,0.495))



pobs_DP <- pobs(agg_data_rmissing_data$diabetes_prevalence)



#Hospital Beds per thousand
fit <- identify_distr(agg_data_rmissing_data$hospital_beds_per_thousand)



histDist(agg_data_rmissing_data$hospital_beds_per_thousand,"GA",nbins = 10)
plotdist(agg_data_rmissing_data$hospital_beds_per_thousand,"GA",para = list(2.95,00.768))



pobs_HBPT <- pobs(agg_data_rmissing_data$hospital_beds_per_thousand)



#Cardiovascular death rate
fit <- identify_distr(agg_data_rmissing_data$cardiovasc_death_rate)



histDist(agg_data_rmissing_data$cardiovasc_death_rate,"IG",nbins = 10)
plotdist(agg_data_rmissing_data$cardiovasc_death_rate,"IG",para = list(257,0.031))



pobs_CDR <- pobs(agg_data_rmissing_data$cardiovasc_death_rate)




pobs_Mortality <- pobs(agg_data_rmissing_data$Mortality)
pobs_death_per_mill <- pobs(agg_data_rmissing_data$total_deaths_per_million)


pobs_cardio <- pobs(agg_data_rmissing_data$cardiovasc_death_rate)
pobs_age <- pobs(agg_data_rmissing_data$median_age)
pobs_diab <- pobs(agg_data_rmissing_data$diabetes_prevalence)
pobs_hosp <- pobs(agg_data_rmissing_data$hospital_beds_per_thousand)
pobs_hdi <- pobs(agg_data_rmissing_data$human_development_index)




############## Compare Correlations ##############

cor.test(pobs_Mortality,pobs_cardio ,method = c("pearson", "kendall", "spearman"))
cor.test(pobs_Mortality,pobs_age, method = c("pearson", "kendall", "spearman"))
cor.test(pobs_Mortality,pobs_diab ,method = c("pearson", "kendall", "spearman"))


cor.test(pobs_Mortality,pobs_hosp ,method = c("pearson", "kendall", "spearman"))
cor.test(pobs_Mortality,pobs_hdi ,method = c("pearson", "kendall", "spearman"))
cor.test(pobs_Mortality,pobs_Mortality ,method = c("pearson", "kendall", "spearman"))


cor.test(agg_data_rmissing_data$total_deaths_per_million,agg_data_rmissing_data$cardiovasc_death_rate ,method = c("pearson", "kendall", "spearman"))
cor.test(agg_data_rmissing_data$total_deaths_per_million,agg_data_rmissing_data$median_age, method = c("pearson", "kendall", "spearman"))
cor.test(agg_data_rmissing_data$total_deaths_per_million,agg_data_rmissing_data$diabetes_prevalence ,method = c("pearson", "kendall", "spearman"))


cor.test(agg_data_rmissing_data$total_deaths_per_million,agg_data_rmissing_data$hospital_beds_per_thousand ,method = c("pearson", "kendall", "spearman"))
cor.test(agg_data_rmissing_data$total_deaths_per_million,agg_data_rmissing_data$human_development_index ,method = c("pearson", "kendall", "spearman"))
cor.test(agg_data_rmissing_data$total_deaths_per_million,pobs_Mortality ,method = c("pearson", "kendall", "spearman"))


cor.test(pobs_death_per_mill,pobs_cardio ,method = c("pearson"))
cor.test(pobs_death_per_mill,pobs_age, method = c("pearson"))
cor.test(pobs_death_per_mill,pobs_diab ,method = c("pearson"))


cor.test(pobs_death_per_mill,pobs_hosp ,method = c("pearson"))
cor.test(pobs_death_per_mill,pobs_hdi ,method = c("pearson"))
cor.test(pobs_Mortality,pobs_Mortality ,method = c("pearson"))




############## Create Plots ##############


cuts = c(0, 0.2, 0.4, 0.6, 0.8, 1)

x_cut <- cut(pobs_age, 20)
y_cut <- cut(pobs_death_per_mill, 20)

xy_table = table(x_cut, y_cut)
hist3D(z=xy_table, border="black", col = "lightsteelblue1", xlab = "CDF of Median Age", ylab = "CDF of Deaths per mill")


cuts = c(0, 0.2, 0.4, 0.6, 0.8, 1)

x_cut <- cut(pobs_hosp, 20)
y_cut <- cut(pobs_death_per_mill, 20)

xy_table = table(x_cut, y_cut)
hist3D(z=xy_table, border="black", col = "lightsteelblue1", xlab = "CDF of Hospital Beds per million", ylab = "CDF of Deaths per mill")

cuts = c(0, 0.2, 0.4, 0.6, 0.8, 1)

x_cut <- cut(pobs_hdi, 20)
y_cut <- cut(pobs_death_per_mill, 20)

xy_table = table(x_cut, y_cut)
hist3D(z=xy_table, border="black", col = "lightsteelblue1", 
       xlab = "CDF of Human Development Index", ylab = "CDF of Deaths per mill")



p <- ggplot(agg_data_rmissing_data, aes(x=median_age, y=total_deaths_per_million, color="blue")) +
  geom_point() +
  theme(legend.position="none")
p1 <- ggMarginal(p, type="histogram", size=3, fill = "lightskyblue2", color = "lightskyblue2")
print(p1)



data <- data.frame(pobs_age, pobs_death_per_mill)
p <- ggplot(data, aes(x=pobs_age, y=pobs_death_per_mill, color="blue")) +
  geom_point() +
  xlab("U(median_age)") + ylab("V(total_deaths_per_million)") + 
  theme(legend.position="none")
p1 <- ggMarginal(p, type="histogram", size=3, fill = "lightskyblue2", color = "lightskyblue2")
print(p1)


p <- ggplot(agg_data_rmissing_data, aes(x=hospital_beds_per_thousand, y=total_deaths_per_million, color="blue")) +
  geom_point() +
  theme(legend.position="none")
p1 <- ggMarginal(p, type="histogram", size=3, fill = "lightskyblue2", color = "lightskyblue2")
print(p1)


data <- data.frame(pobs_hosp, pobs_death_per_mill)
p <- ggplot(data, aes(x=pobs_hosp, y=pobs_death_per_mill, color="blue")) +
  geom_point() +
  xlab("U(hospital_beds)") + ylab("V(total_deaths_per_million)") +
  theme(legend.position="none")
p1 <- ggMarginal(p, type="histogram", size=3, fill = "lightskyblue2", color = "lightskyblue2")
print(p1)


p <- ggplot(agg_data_rmissing_data, aes(x=human_development_index, y=total_deaths_per_million, color="blue")) +
  geom_point() +
  theme(legend.position="none")
p1 <- ggMarginal(p, type="histogram", size=3, fill = "lightskyblue2", color = "lightskyblue2")
print(p1)


data <- data.frame(pobs_hdi, pobs_death_per_mill)
p <- ggplot(data, aes(x=pobs_hdi, y=pobs_death_per_mill, color="blue")) +
  geom_point() +
  xlab("U(HDI)") + ylab("V(total_deaths_per_million)") +
  theme(legend.position="none")
p1 <- ggMarginal(p, type="histogram", size=3, fill = "lightskyblue2", color = "lightskyblue2")
print(p1)




############## Fit Copulas and create plots ##############


# pobs_Mortality,pobs_age

selectedCopula <- BiCopSelect(pobs_age , pobs_death_per_mill)
selectedCopula

cp <- frankCopula(par = 4.28, dim = 2)
persp(cp, dCopula, col = "lightskyblue2", xlab = "U(median age)", ylab = 'V(death per million)')

u <- rCopula(3965,cp)
cor(u, method = c("pearson"))


# pobs_Mortality,pobs_age

selectedCopula <- BiCopSelect(pobs_hosp , pobs_death_per_mill)
selectedCopula

# cp <- BB8Copula(param = c(3.27, 0.84))
cp <- copulaFromFamilyIndex(10, par = 3.27, par2 = 0.84)
persp(cp, dCopula, col = "lightskyblue2", xlab = "U(hospital beds per thousands)", ylab = 'V(death per million)')

u <- rCopula(3965,cp)
cor(u, method = c("pearson"))


# pobs_Mortality,pobs_age

selectedCopula <- BiCopSelect(pobs_hdi , pobs_death_per_mill)
selectedCopula

# cp <- BB8Copula(param = c(3.04, 0.95))
cp <- copulaFromFamilyIndex(10, par = 3.04, par2 = 0.95)


persp(cp, dCopula, col = "lightskyblue2", xlab = "U(HDI)", ylab = 'V(death per million)')

u <- rCopula(3965,cp)
cor(u, method = c("pearson"))