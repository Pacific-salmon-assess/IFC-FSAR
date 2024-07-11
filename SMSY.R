
library(tidyverse)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")



setwd(cohoDir)

SRDat <- read.csv("./DataIn/IFCoho_SRbyCU.csv")

muLSurv <- SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
muLSurv <- muLSurv$muLSurv[1]

2.88 + 0.42*muLSurv
A <- exp(2.88 + 0.42*muLSurv)

(1 - gsl::lambert_W0(exp(1 - log(A)))) / 0.000121

ricksurv <- read.csv("./SamSimInputs/SR_IndivRicker_Surv_mcmc.csv")
rickpc <- read.csv("./SamSimInputs/SR_IndivRicker_SurvCap_mcmc.csv")

rickcombo <- rbind(ricksurv,rickpc)
str(rickcombo)

rickcombo <- rickcombo %>% 
  mutate(logA = alpha + gamma*muLSurv,
         A = exp(logA),
         SMSY = (1 - gsl::lambert_W0(exp(1 - logA))) / beta,
         Smsy80 = SMSY*0.80)

benchmarks <-rickcombo %>% 
  group_by(stk) %>% 
  summarise(A.mean = mean(A),
            A.SD = sd(A),
            A.05 = quantile(A, 0.05),
            A.10 = quantile(A, 0.10),
            A.25 = quantile(A, 0.25),
            A.50 = quantile(A, 0.50),
            A.75 = quantile(A, 0.75),
            A.90 = quantile(A, 0.90),
            A.95 = quantile(A, 0.95),
            Sgen.mean = mean(Sgen),
            Smsy80.mean = mean(Smsy80),
            Sgen.SD = sd(Sgen),
            Sgen.05 = quantile(Sgen, 0.05),
            Sgen.10 = quantile(Sgen, 0.10),
            Sgen.25 = quantile(Sgen, 0.25),
            Sgen.50 = quantile(Sgen, 0.50),
            Sgen.75 = quantile(Sgen, 0.75),
            Sgen.90 = quantile(Sgen, 0.90),
            Sgen.95 = quantile(Sgen, 0.95),
            Smsy80.SD = sd(Smsy80),
            Smsy80.05 = quantile(Smsy80, 0.05),
            Smsy80.10 = quantile(Smsy80, 0.10),
            Smsy80.25 = quantile(Smsy80, 0.25),
            Smsy80.50 = quantile(Smsy80, 0.50),
            Smsy80.75 = quantile(Smsy80, 0.75),
            Smsy80.90 = quantile(Smsy80, 0.90),
            Smsy80.95 = quantile(Smsy80, 0.95))
benchmarks

SRDat <- rename(SRDat, stk = CU_ID)
benchmarks <- unique(select(SRDat, stk, CU_Name)) %>% left_join(benchmarks)

mean(c(1646, 2515))
mean(c(314, 304))
mean(c(1977, 2781))
mean(c(2482, 3171))
mean(c(2573, 4050))

bench1 <- select(benchmarks, CU_Name, Sgen.mean, Sgen.SD:Sgen.95) %>% 
  mutate(Parameter = "Sgen")
bench2 <- select(benchmarks, CU_Name, Smsy80.mean, Smsy80.SD:Smsy80.95) %>% 
  mutate(Parameter = "80%Smsy")
bench3 <- select(benchmarks, CU_Name, A.mean:A.95) %>% 
  mutate(Parameter = "Alpha.Prime")

colnames(bench1) <- c("CU_Name", "Mean", "SD", "P05", "P10", "P25", "P50", 
                      "P75", "P90", "P95", "Parameter")
colnames(bench2) <- c("CU_Name", "Mean", "SD", "P05", "P10", "P25", "P50", 
                      "P75", "P90", "P95", "Parameter")
colnames(bench3) <- c("CU_Name", "Mean", "SD", "P05", "P10", "P25", "P50", 
                      "P75", "P90", "P95", "Parameter")

source("C:/Users/arbeiderm/Documents/R/Round R Function.R")

bench <- rbind(bench1, bench2, bench3) %>% 
  select(CU_Name, Parameter, Mean:P95) %>% 
  arrange(CU_Name)

bench[ , 3:11] <- round2(bench[ , 3:11], 2)

write.csv(bench, "./DataOut/Updated_IFR_Coho_WSP_Benchmarks.csv")

write.csv(bench, "./DataOut/Updated_IFR_Coho_WSP_Benchmarks with Alpha Prime.csv")
