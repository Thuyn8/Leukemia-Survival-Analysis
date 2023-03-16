#Packages

library(survival)

library(dplyr)

library(tidyverse)

library(survminer)

library(gtsummary)

library(sjPlot)

library(readr)



bmt <- read.csv("~/Documents/survival analysis/project/bmt.csv")



#loading the file

bmt <- read_csv("bmt-2.csv", show_col_types = FALSE)



#Q1----

#"disease-free survival"



s.bmt<-with(bmt,Surv(time=tdfs, event = deltadfs))

survfit.bmt<-survfit(s.bmt~1,data=bmt,conf.type="log-log")



# the median disease free survival for the enrolled patients is 481 days 95 CI 363-748 days.



plot(survfit.bmt, conf.type = "log-log",
     
     main = "Disease Free Survival Curve From Bone Marrow Transplant ",
     
     col = c(2,4),
     
     lty = 1,
     
     lwd = 2,
     
     conf.int = F,
     
     xlab = "Time (in days)",
     
     ylab = "Disease Free Survival Probability",
     
     cex.main=0.65)



#Q2---- table1 (2 tbls)

table1 <- bmt %>% select("age","male","cmv","donorage","donormale","donorcmv","waittime","fab","hospital", "disgroup", "mtx")

table1 %>% tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"))

table1 %>% tbl_summary(by = fab, statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)")) %>% add_p()%>%as_flex_table()



table1 %>% tbl_summary(by = disgroup, statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)")) %>% add_p()%>%
  
  as_flex_table()



#Q3 "baseline risk factors"

cox.bmt.age <- coxph(s.bmt ~ age , data = bmt)

summary(cox.bmt.age)



cox.bmt.male <- coxph(s.bmt ~ as.factor(male) , data = bmt)

summary(cox.bmt.male)



cox.bmt.cmv <- coxph(s.bmt ~ as.factor(cmv), data = bmt)

summary(cox.bmt.cmv)



cox.bmt.donorage <- coxph(s.bmt ~ donorage , data = bmt)

summary(cox.bmt.donorage)



cox.bmt.donormale <- coxph(s.bmt ~ as.factor(donormale) , data = bmt)

summary(cox.bmt.donormale)



cox.bmt.donorcmv <- coxph(s.bmt ~ as.factor(donorcmv) , data = bmt)

summary(cox.bmt.donorcmv)



cox.bmt.waittime <- coxph(s.bmt ~ waittime, data = bmt)

summary(cox.bmt.waittime)



cox.bmt.fab <- coxph(s.bmt ~ as.factor(fab) , data = bmt)

summary(cox.bmt.fab)



cox.bmt.dis <- coxph(s.bmt ~ as.factor(disgroup) , data = bmt)

summary(cox.bmt.dis)



cox.bmt.hospital <- coxph(s.bmt ~ as.factor(hospital), data = bmt)

cox.bmt.hospital



cox.bmt.mtx <- coxph(s.bmt ~ as.factor(mtx) , data = bmt)

summary(cox.bmt.mtx)



cox.bmt <- coxph(s.bmt ~ age + as.factor(male) + as.factor(cmv) +
                   
                   donorage + as.factor(donormale) + as.factor(donorcmv) +
                   
                   waittime + as.factor(fab) + as.factor(disgroup)+as.factor(mtx), data = bmt)

summary(cox.bmt)



tab_model(cox.bmt,transform = "exp", auto.label = FALSE, collapse.ci=TRUE, file = "Cox.dfs.xls")





# using an alpha of 0.10 to compensate for the small sample size, it appears that FAB and hospital

# were associated with different disease free survival.



#Q4

## disease-free survival

bmt.del <- tmerge(data1 = bmt,
                  
                  data2 = bmt,
                  
                  id = id,
                  
                  death = event(tdfs, deltadfs),
                  
                  tx.tv = tdc(ta))

bmt.del$same_sex <- ifelse(bmt.del$male == bmt.del$donormale,1,0)



surv.bmt.del <- with(bmt.del, Surv(tstart, tstop, death))

cox.bmt.del <- coxph(surv.bmt.del~ tx.tv + 
                       as.factor(disgroup) + 
                       age + 
                       as.factor(male) +
                       donorage +
                       as.factor(donormale) +
                       as.factor(cmv)+
                       as.factor(donorcmv) +
                       as.factor(mtx) +
                       as.factor(hospital) +
                       waittime,  
                     data = bmt.del) 

summary(cox.bmt.del)



## risk of relapse

bmt.del.relapse <- tmerge(data1 = bmt,
                          
                          data2 = bmt,
                          
                          id = id,
                          
                          death = event(tdfs, deltar),
                          
                          tx.tv = tdc(ta))



bmt.del.relapse$same_sex <- ifelse(bmt.del.relapse$male == bmt.del.relapse$donormale,1,0)



surv.bmt.relapse <- with(bmt.del.relapse, Surv(tstart, tstop, death))

cox.bmt.relapse <- coxph(surv.bmt.relapse~ tx.tv + 
                           as.factor(disgroup) + 
                           age + 
                           as.factor(male) +
                           donorage +
                           as.factor(donormale) +
                           as.factor(cmv)+
                           as.factor(donorcmv) +
                           as.factor(mtx) +
                           as.factor(hospital) +
                           waittime, data= bmt.del.relapse) 
summary(cox.bmt.relapse)



#Q5:

s.object1<-with(bmt, Surv(tdfs,deltadfs))



mod5<-coxph(s.object1 ~as.factor(disgroup)+ age + as.factor(male) +donorage+
              
              as.factor(donormale)+as.factor(cmv)+as.factor(donorcmv)
            
            +as.factor(mtx)+as.factor(hospital)+waittime, data = bmt, subset=(deltaa==1))

summary(mod5)



## mtx survival

test<-survfit(s.bmt~mtx,data=bmt,conf.type="log-log")



ggsurvplot(
  
  test,                    
  
  data = bmt,            
  
  risk.table = TRUE,      
  
  pval = TRUE,            
  
  conf.int = TRUE,        
  
  xlim = c(0,365),
  
  xlab = "Time in days",
  
  ylab = "Disease free survival",
  
  title= "Post transplant disease free survival",
  
  break.time.by = 90,   
  
  ggtheme = theme_light(),
  
  risk.table.y.text.col = T,
  
  risk.table.y.text = FALSE,
  
  legend.labs =
    
    c("no MTX", "MTX")   
  
)



table(bmt$mtx)



##

# QUESTION 6

# Survival object

s.object3<-with(bmt, Surv(ta,deltaa))



# Model

mod7<-coxph(s.object1 ~as.factor(mtx)+as.factor(disgroup)+age+as.factor(male)
            
            +donorage+as.factor(donormale)+as.factor(cmv)+as.factor(donorcmv) +as.factor(hospital),data = bmt)

summary(mod7)



## KM

plot(survfit(s.object3 ~as.factor(mtx) , data=bmt, conf.type = "log-log"),
     
     main = "Kaplan-Meier Estimators Comparing the Survival Times From Transplant Until

onset of Acute Graff Versus Host Disease",col = c(2,4), lty = 1, lwd = 2,conf.int = F,xlab = "Time (in days)", ylab = "Survival Probability",cex.main=0.65)

legend("topright",legend=c("Received Methotrexate","Didn't Receive Methotrexate"),fill=c(2,4),col=c(2,4),cex=0.7,bty="n")
