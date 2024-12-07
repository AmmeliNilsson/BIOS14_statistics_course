  dat = read.table("../Data/exam2022_part2-1.txt", header=T) #reading the data
library(psych)
library(PerformanceAnalytics)
numeric_data <- dat[, sapply(dat, is.numeric)]
pairs.panels(numeric_data)
chart.Correlation(numeric_data, histogram=TRUE, pch=19) #many histograms to compare
#model w all variables, remove stepwise nonsignificant
Saturated_model <- lm(mass ~ sex + hornR + month + day + yr + daynr + age + cohort + hornL + density, data = dat)
summary(Saturated_model)
cor(dat[, sapply(dat, is.numeric)], use = "pairwise.complete.obs") #cohort is correlated w season, remove one of them
m1 <- lm(mass ~ sex + hornR + month + day + yr + daynr + cohort + hornL + density, data = dat)
summary(m1)
m2 <- lm(mass ~ sex + hornR + month + yr + daynr + cohort + hornL + density, data = dat)
summary(m2) #all significant!
m <- lm(mass ~ sex + month + yr + daynr + cohort + hornL + density, data = dat)
summary(m) #all varaibles significant
model <- lm(mass ~ sex + hornL + density + cohort, na.rm=T, data = dat) #choosing to remove varaibles that don not seem biologically important
summary(model)
#checking which model is best based on AIC and Weight
model1 <- lm(mass ~ sex + hornL + density + cohort, data = dat) 
model2 <- lm(mass ~ sex + hornL + cohort, data = dat)        
model3 <- lm(mass ~ sex + density + cohort, data = dat)    
model4 <- lm(mass ~ hornL + density + cohort, data = dat)     
model5 <- lm(mass ~ sex + hornL + density, data = dat)       
model6 <- lm(mass ~ sex + hornL, data = dat)                 
model7 <- lm(mass ~ sex + density, data = dat)           
model8 <- lm(mass ~ sex + cohort, data = dat)             
model9 <- lm(mass ~ hornL + density, data = dat)         
model10 <- lm(mass ~ hornL, data = dat)                         
model11 <- lm(mass ~ cohort, data = dat)                        
model12 <- lm(mass ~ 1, data = dat)                              
model_list <- list(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12)
AIC_table <- AIC(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12)
AIC_table$logLik <- unlist(lapply(model_list, logLik))
AIC_table <- AIC_table[order(AIC_table$AIC), ]
AIC_table$delta <- round(AIC_table$AIC - min(AIC_table$AIC), 2)
likelihood <- exp(-0.5 * AIC_table$delta)
AIC_table$weight <- round(likelihood / sum(likelihood), 2)
print(AIC_table)

library(glmmTMB)
model_glmmTMB <- glmmTMB(mass ~ sex + hornL + density + (1 | cohort), data = dat) #general linear mixed model
#comparing mixed models:
m1_glmmTMB <- glmmTMB(mass ~ hornL * sex + density + (1 | cohort), data = dat)
m2_glmmTMB <- glmmTMB(mass ~ hornL + sex + density + (1 | cohort), data = dat)
m3_glmmTMB <- glmmTMB(mass ~ hornL + density + (1 | cohort), data = dat)
m4_glmmTMB <- glmmTMB(mass ~ sex + density + (1 | cohort), data = dat)
m5_glmmTMB <- glmmTMB(mass ~ density + (1 | cohort), data = dat)
m6_glmmTMB <- glmmTMB(mass ~ hornL + (1 | cohort), data = dat)
m7_glmmTMB <- glmmTMB(mass ~ sex + (1 | cohort), data = dat)
m8_glmmTMB <- glmmTMB(mass ~ 1 + (1 | cohort), data = dat) 
mlist_glmmTMB <- list(m1_glmmTMB, m2_glmmTMB, m3_glmmTMB, m4_glmmTMB, m5_glmmTMB, m6_glmmTMB, m7_glmmTMB, m8_glmmTMB)
AIC_glmmTMB_table <- AIC(m1_glmmTMB, m2_glmmTMB, m3_glmmTMB, m4_glmmTMB, m5_glmmTMB, m6_glmmTMB, m7_glmmTMB, m8_glmmTMB)
AIC_glmmTMB_table$delta <- AIC_glmmTMB_table$AIC - min(AIC_glmmTMB_table$AIC)
AIC_glmmTMB_table$weight <- exp(-0.5 * AIC_glmmTMB_table$delta) / sum(exp(-0.5 * AIC_glmmTMB_table$delta))
print(AIC_glmmTMB_table)
#Comparing with REML=F
m1 <- glmmTMB(mass ~ hornL * sex + density + (1 | cohort), data = dat, REML = FALSE)  
m2 <- glmmTMB(mass ~ hornL + sex + density + (1 | cohort), data = dat, REML = FALSE) 
m3 <- glmmTMB(mass ~ hornL + density + (1 | cohort), data = dat, REML = FALSE)      
m4 <- glmmTMB(mass ~ sex + density + (1 | cohort), data = dat, REML = FALSE)    
m5 <- glmmTMB(mass ~ density + (1 | cohort), data = dat, REML = FALSE)            
m6 <- glmmTMB(mass ~ hornL + (1 | cohort), data = dat, REML = FALSE)        
m7 <- glmmTMB(mass ~ sex + (1 | cohort), data = dat, REML = FALSE)          
m8 <- glmmTMB(mass ~ 1 + (1 | cohort), data = dat, REML = FALSE)         
mlist <- list(m1, m2, m3, m4, m5, m6, m7, m8)
AIC_fixed_table <- AIC(m1, m2, m3, m4, m5, m6, m7, m8)
AIC_fixed_table$delta <- AIC_fixed_table$AIC - min(AIC_fixed_table$AIC)
AIC_fixed_table$weight <- exp(-0.5 * AIC_fixed_table$delta) / sum(exp(-0.5 * AIC_fixed_table$delta))
print(AIC_fixed_table)

library(Matrix)
library(glmmTMB)
str(dat)
Mixed_Model <- glmmTMB(mass ~ hornL * sex + density + (1 | cohort), data = dat, REML = TRUE) #Final model
summary(Mixed_Model)

coef(Mixed_Model) #all years

library(MuMIn)
r_squared_Mixed_Model <- r.squaredGLMM(Mixed_Model) # Marginal R²: The proportion of variance explained by the fixed effects:0.536
#Conditional R²: The proportion of variance explained by both the fixed effects and the random effects: 0.547
r_squared_Mixed_Model
#Do this for every varaible, removing one factor at a time, take minus to know the variance explained by it
M_density <- glmmTMB(mass ~ hornL * sex + (1 | cohort), data = dat, REML = TRUE)
r_squared_M_density <- r.squaredGLMM(M_density) # 0.524%
r_squared_M_density
M_sex <- glmmTMB(mass ~ hornL * density + (1 | cohort), data = dat, REML = TRUE)
r_squared_M_sex <- r.squaredGLMM(M_sex) #2.47%
r_squared_M_sex
M_hornL <- glmmTMB(mass ~ sex + density + (1 | cohort), data = dat, REML = TRUE)
r_squared_M_hornL <- r.squaredGLMM(M_hornL) #25%
r_squared_M_hornL

#Check that residuals are normally distrtributed, Yes!
qqnorm(residuals(Mixed_Model))
qqline(residuals(Mixed_Model))

#mass both sexes
library(dplyr) 
dat <- dat %>% mutate(sex = factor(sex), density = factor(density))
n_mass <- length(na.omit(dat$mass)) #Obs: n= 4394
mean_mass <- mean(dat$mass, na.rm = TRUE) #22.3
min_mass <- min(dat$mass, na.rm = TRUE)#43
max_mass <- max(dat$mass, na.rm = TRUE) #5
sd_mass <- sd(dat$mass, na.rm = TRUE) #5.39
sem_mass <- sd_mass / sqrt(length(na.omit(dat$mass))) #standard error: 0.0813
#mass females
dat_females <- dat %>% filter(sex == "F")
n_female_mass <- length(na.omit(dat_females$mass)) #1955 
mean_female_mass <- mean(dat_females$mass, na.rm = TRUE) #20.6
min_female_mass <- min(dat_females$mass, na.rm = TRUE) #6.5 
max_female_mass <- max(dat_females$mass, na.rm = TRUE)#35
sd_female_mass <- sd(dat_females$mass, na.rm = TRUE) #4.05
sem_female_mass <- sd_female_mass / sqrt(n_female_mass) #0.0915
#mass males
dat_males <- dat %>% filter(sex == "M")
n_male_mass <- length(na.omit(dat_males$mass)) #2439 
mean_male_mass <- mean(dat_males$mass, na.rm = TRUE) #23.7
min_male_mass <- min(dat_males$mass, na.rm = TRUE)#5
max_male_mass <- max(dat_males$mass, na.rm = TRUE)#43 
sd_male_mass <- sd(dat_males$mass, na.rm = TRUE) #5.90
sem_male_mass <- sd_male_mass / sqrt(n_male_mass)#0.1196

#horn length both sexes
mean_hornL <- mean(dat$hornL, na.rm = TRUE) #181.7
min_hornL <- min(dat$hornL, na.rm = TRUE) #10
max_hornL <- max(dat$hornL, na.rm = TRUE) #295
sd_hornL <- sd(dat$hornL, na.rm = TRUE) #44.0
sem_hornL <- sd_hornL / sqrt(length(na.omit(dat$hornL))) #0.664
#horn length females
mean_female_hornL <- mean(dat_females$hornL, na.rm = TRUE)#169
min_female_hornL <- min(dat_females$hornL, na.rm = TRUE) #11
max_female_hornL <- max(dat_females$hornL, na.rm = TRUE) #280
sd_female_hornL <- sd(dat_females$hornL, na.rm = TRUE) #41.3
sem_female_hornL <- sd_female_hornL / sqrt(length(na.omit(dat_females$hornL)))#0.935
#horn lenth males
mean_male_hornL <- mean(dat_males$hornL, na.rm = TRUE) #191.9
min_male_hornL <- min(dat_males$hornL, na.rm = TRUE) #10
max_male_hornL <- max(dat_males$hornL, na.rm = TRUE) #295 
sd_male_hornL <- sd(dat_males$hornL, na.rm = TRUE) # 43.5
sem_male_hornL <- sd_male_hornL / sqrt(length(na.omit(dat_males$hornL))) #0.880

# Density
dat$density <- ifelse(dat$density == "low", 0, ifelse(dat$density == "high", 1, NA)) #low=0 high=1
summary(dat$density)
mean_density <- mean(dat$density, na.rm = TRUE) #0.514
min_density <- min(dat$density, na.rm = TRUE)#0
max_density <- max(dat$density, na.rm = TRUE)#1
sd_density <- sd(dat$density, na.rm = TRUE) #0.4999
sem_density <- sd_density / sqrt(length(na.omit(dat$density)))#0.007541

library(ggplot2)
#Figure: Number of chamois collected for each cohort by sex
dat <- read.table("../Data/exam2022_part2-1.txt", header = TRUE)
summary_data <- aggregate(dat$sex, by = list(Cohort = dat$cohort, Sex = dat$sex), FUN = length)
colnames(summary_data) <- c("Cohort", "Sex", "Count")
library(ggplot2)
ggplot(summary_data, aes(x = factor(Cohort), y = Count, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Number of Chamois collected for each cohort by sex",
       x = "Cohort",
       y = "Number of Chamois",
       fill = "Sex") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("F" = "salmon", "M" = "skyblue"))

#Figure: Horn length & body mass by sex
ggplot(dat, aes(x = hornL, y = mass, color = sex)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE, aes(color = sex)) + 
  labs(title = "Horn length & body mass by sex", 
       x = "Horn length (mm)", 
       y = "Body mass (kg)") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("F" = "salmon", "M" = "skyblue"))

#Figure: Horn length & body mass by density
dat$density <- as.factor(dat$density)
ggplot(dat, aes(x = hornL, y = mass, color = density)) +
  geom_point(alpha = 0.7) +               
  geom_smooth(method = "lm", se = TRUE) +  
  scale_color_manual(values = c("darkblue", "6fa3d3")) +  
  labs(
    title = "Horn length & body mass by density",
    x = "Horn length (mm)",
    y = "Body mass (kg)",
    color = "Density"
  ) +
  theme_minimal()

