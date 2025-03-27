library(tidyverse)
library(here)
library(rio)
library(janitor)
library(survival)
library(survminer)
library(gtsummary)
library(naniar)
library(ggpubr)
library(glmnet)
library(mice)
library(rstatix)
data=import("Data.xlsx",which = "Clean") %>% 
  clean_names() %>% 
  mutate(v2_Fib_recat=ifelse(v2_fibro_scr_e2_c2%in%c("F0","F1","F2"),"F0-F2",
                             ifelse(v2_fibro_scr_e2_c2%in%c("F3","F4"),"F3-F4",
                                    ifelse (v2_fibro_scr_e2_c2=="NA","NA",NA))),
         v1_Fib_recat=ifelse(v1_fibro_scr_e1_c1%in%c("F0","F1","F2"),"F0-F2",
                             ifelse(v1_fibro_scr_e1_c1%in%c("F3","F4"),"F3-F4",
                                    ifelse (v1_fibro_scr_e1_c1=="NA","NA",NA))),
         v3_Fib_recat=ifelse(v3_fibro_scr_e3_c3%in%c("F0","F1","F2"),"F0-F2",
                             ifelse(v3_fibro_scr_e3_c3%in%c("F3","F4"),"F3-F4",
                                    ifelse (v3_fibro_scr_e3_c3=="NA","NA",NA))),
         CPC_RECAT=ifelse(v1_pbc_cpc_e1_c1>=5 & v1_pbc_cpc_e1_c1<=6,"Class A",
                          ifelse(v1_pbc_cpc_e1_c1>6 ,"Class B","Normal")),
         v2_decomp_scr=rowSums(across(c(v2_ascites,v2_varices,v2_variceal_bleeding,v2_encephalopathy,v2_lt,v2_hcc))),
         v3_decomp_scr=rowSums(across(c(v3_ascites,v3_varices,v3_variceal_bleeding,v3_encephalopathy,v3_lt,v3_hcc))),
         v2_decomp_stat=ifelse(v2_decomp_scr>0,1,0),
         v3_decomp_stat=ifelse(v3_decomp_scr>0,1,0),
         AMA_RECAT=ifelse(str_detect(v1_lab9_u_e1_c1,"Positive"),"Positive","Negative")
         )
         

# Age at diagnosis post hoc:
age=data %>% 
  select(new_subject_id,contains("AGE"),contains("COB")) %>% 
  mutate(country=as.factor(substr(new_subject_id,1,2)))

pairwise.t.test(age[["v1_pbc_dia_age_e1_c1"]], age[["country"]], p.adj = "bonferroni")

#### Table1####
Tab1=data %>% 
  select(v2_ttt_resp_e2_c2_1,v1_pbc_dia_age_e1_c1,sex_e1_c1,race_e1_c1,v1_pbc_sym_e1_c1,
         v1_lab1_e1_c1,v1_lab2_e1_c1,v1_lab3_e1_c1,v1_lab4_e1_c1,v1_lab5_e1_c1,v1_lab6_e1_c1,v1_lab7_e1_c1,v1_lab8_e1_c1,
         v1_lab9_e1_c1,v1_lab10_e1_c1,v1_lab11_e1_c1,v1_lab12_e1_c1,v1_lab13_e1_c1,v1_lab14_e1_c1,
         v1_pbc_cpc_e1_c1,v1_Fib_recat,v1_hri_class_e1_c1,v1_lab9_u_e1_c1,CPC_RECAT,v1_pbc_sym_e1_c1,v2_decomp_stat)
#continuous vars:
Tab1 %>% 
  select(v2_ttt_resp_e2_c2_1,v1_pbc_dia_age_e1_c1,starts_with("v1_lab")) %>% 
  select(-c(v1_lab9_e1_c1,v1_lab10_e1_c1,v1_lab11_e1_c1)) %>% 
  tbl_summary(
    by = v2_ttt_resp_e2_c2_1,
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{mean} ({sd})",
      "{median} ({p25}, {p75})"
    ),
    digits = all_continuous() ~ 2,
  ) %>% 
  add_n() %>% 
  add_p()

#Categorical vars:
Tab1_cat=data %>% 
  select(v2_ttt_resp_e2_c2_1,sex_e1_c1,race_e1_c1,v1_pbc_sym_e1_c1,v1_Fib_recat,
         v1_hri_class_e1_c1,v1_lab9_u_e1_c1,CPC_RECAT)

Tab1_cat %>% 
  tbl_summary(
    by=v2_ttt_resp_e2_c2_1
  ) %>% 
  add_n() %>% 
  add_overall() %>% 
  add_p()

##### Table2 inquiries####
tab2=data %>% 
  select(new_subject_id,v2_ascites,v2_varices,v2_variceal_bleeding,v2_lt,v2_encephalopathy,v2_hcc,
         v3_ascites,v3_varices,v3_variceal_bleeding,v3_lt,v3_encephalopathy,v3_hcc) %>% 
  mutate_all(~replace_na(., 0))
  

tab2 %>% 
  select(-new_subject_id) %>% 
  tbl_summary()

chk_ascites=tab2 %>% 
  select(new_subject_id,v2_ascites,v3_ascites) %>% 
  filter(v2_ascites==1 |v3_ascites==1 )

chk_varices=tab2 %>% 
  select(new_subject_id,v2_varices,v3_varices) %>% 
  filter(v2_varices==1 |v3_varices==1 )

chk_encephalopathy=tab2 %>% 
  select(new_subject_id,v2_encephalopathy,v3_encephalopathy) %>% 
  filter(v2_encephalopathy==1 |v3_encephalopathy==1)

chk_lt=tab2 %>% 
  select(new_subject_id,v2_lt,v3_lt) %>% 
  filter(v2_lt==1 |v3_lt==1)

chk_hcc=tab2 %>% 
  select(new_subject_id,v2_hcc,v3_hcc) %>% 
  filter(v2_hcc==1 |v3_hcc==1)

export(chk_ascites,"chk_ascites.xlsx")
export(chk_varices,"chk_varices.xlsx")
export(chk_encephalopathy,"chk_encephalopathy.xlsx")
export(chk_lt,"chk_lt.xlsx")
export(chk_hcc,"chk_hcc.xlsx")


## end of follow up clinical outcome by treatment response:
data1=data %>% 
  select (new_subject_id,v3_ttt_e3_c3_1,v3_ttt_resp_e3_c3_1) %>% 
  filter(!v3_ttt_e3_c3_1 %in% c("Obeticholic Acid","Obeticholic Acid 5m PO Od x","tacrolimus","prescribed","spironolactone","Paracetamol: 1,000 mg","omprazole",NA))

table(data1$v3_ttt_resp_e3_c3_1)

## Median (IQR) duration of follow up for the overall cohort :
data2=data %>% 
  select(new_subject_id,v1_dt_e1_c1_10,v3_out_stat_e3_c3,v3_out_fu_dt_e3_c3,v3_out_dod_e3_c3) %>% 
  mutate(
    vis1dt=as.Date(v1_dt_e1_c1_10,origin="1899-12-30"),
    vis3fudt=as.Date(v3_out_fu_dt_e3_c3,origin="1899-12-30"),
    dod=as.Date(v3_out_dod_e3_c3,origin="1899-12-30"),
    vis3fudtn=coalesce(dod,vis3fudt),
    FUDUR=vis3fudtn-vis1dt,
    DURN=as.numeric(FUDUR)
  )
fivenum(data2$DURN)

###### OCA patients #######
OCA_V2=data %>% 
  filter (v2_ttt_e2_c2_2 %in% c("Obeticholic Acid","obeticholic acid")) %>% 
  select(v2_rri1_e2_c2,v2_rri2_e2_c2,v2_ttt_resp_e2_c2_2,v1_pbc_dia_age_e1_c1,sex_e1_c1,race_e1_c1,race_sp_e1_c1,v1_pbc_sym_e1_c1,v2_lab1_e2_c2,v2_lab2_e2_c2,v2_lab3_e2_c2,
         v2_lab4_e2_c2,v2_lab5_e2_c2,v2_lab6_e2_c2,v2_lab7_e2_c2,v2_lab8_e2_c2,v2_fibro_scr_e2_c2,v2_hri_class_e2_c2) %>% 
  mutate(v2_lab8_e2_c2=as.numeric(v2_lab8_e2_c2))
 
OCA_V2 %>% 
  tbl_summary(
    by=v2_ttt_resp_e2_c2_2,
    type = c(v1_pbc_dia_age_e1_c1,v2_lab1_e2_c2,v2_lab2_e2_c2,v2_lab3_e2_c2,v2_lab4_e2_c2,v2_lab5_e2_c2,v2_lab6_e2_c2,v2_lab7_e2_c2,v2_lab8_e2_c2) ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{mean} ± ({sd})",
      "{median} ({p25}, {p75})"
    ),
    missing = "no"
  ) %>% 
  add_n()




OCA_V3=data %>% 
  filter (v3_ttt_e3_c3_2 %in% c("Obeticholic Acid","obeticholic acid","Obeticholic acid"))%>% 
  select( v3_rri1_e3_c3 ,v3_rri2_e3_c3,v3_ttt_resp_e3_c3_2,v1_pbc_dia_age_e1_c1,sex_e1_c1,race_e1_c1,race_sp_e1_c1,v1_pbc_sym_e1_c1,v3_lab1_e3_c3,v3_lab2_e3_c3,v3_lab3_e3_c3,
         v3_lab4_e3_c3,v3_lab5_e3_c3,v3_lab6_e3_c3,v3_lab7_e3_c3,v3_lab8_e3_c3,v3_fibro_scr_e3_c3,v3_hri_class_e3_c3)%>% 
  mutate(v3_lab5_e3_c3=as.numeric(v3_lab5_e3_c3),
         v3_lab6_e3_c3=as.numeric(v3_lab6_e3_c3))

OCA_V3 %>% 
  tbl_summary(
    by=v3_ttt_resp_e3_c3_2,
    type = c(v1_pbc_dia_age_e1_c1,v3_lab1_e3_c3,v3_lab2_e3_c3,v3_lab3_e3_c3,
             v3_lab4_e3_c3,v3_lab5_e3_c3,v3_lab6_e3_c3,v3_lab7_e3_c3,v3_lab8_e3_c3) ~ "continuous2",
    digits = v3_lab4_e3_c3 ~2,
    statistic = all_continuous() ~ c(
      "{mean} ± ({sd})",
      "{median} ({p25}, {p75})"
    ),
    missing = "ifany"
  ) %>% 
  add_n()

#### survival analysis #####
data3=data %>% 
  select(new_subject_id,v1_fibro_find_e1_c1,v1_fh_pbc_e1_c1,race_e1_c1, v1_surg_e1_c1,v1_pbc_criteria_e1_c1,v1_hri_class_e1_c1,sex_e1_c1,
         age_e1_c1,v1_pbc_dia_age_e1_c1,v1_dt_e1_c1_10,v3_out_stat_e3_c3,v3_out_fu_dt_e3_c3,v3_out_dod_e3_c3,
         v3_ascites,v3_varices,v3_lt,v3_variceal_bleeding,v3_hcc,v3_encephalopathy,
         v3_ttt_resp_e3_c3_1,v3_decomp_stat,v2_decomp_stat,v1_Fib_recat,CPC_RECAT,AMA_RECAT,
         v1_fibro_find_e1_c1,sex_e1_c1,v1_pbc_dia_age_e1_c1,v1_fh_pbc_e1_c1,race_e1_c1,v1_pbc_cpc_e1_c1,v1_fibro_scr_e1_c1,v1_hri_class_e1_c1,v2_ttt_resp_e2_c2_1,v2_dt_e2_c2_196,
         v1_lab1_e1_c1,v1_lab2_e1_c1,v1_lab3_e1_c1,v1_lab4_e1_c1,v1_lab5_e1_c1,v1_lab6_e1_c1,v1_lab7_e1_c1,v1_lab8_e1_c1,v1_lab9_e1_c1,v1_lab10_e1_c1,v1_lab11_e1_c1,v1_lab12_e1_c1,v1_lab13_e1_c1,v1_lab14_e1_c1) %>% 
  mutate(
    vis1dt=as.Date(v1_dt_e1_c1_10,origin="1899-12-30"),
    vis2dt=as.Date(v2_dt_e2_c2_196,origin="1899-12-30"),
    vis3fudt=as.Date(v3_out_fu_dt_e3_c3,origin="1899-12-30"),
    dod=as.Date(v3_out_dod_e3_c3,origin="1899-12-30"),
    vis3fudtn=coalesce(dod,vis3fudt),
    FUDUR=vis3fudtn-vis1dt,
    FU2DUR=vis2dt-vis1dt,
    DURN=ceiling(as.numeric(FUDUR)/30),
    DURN2=ceiling(as.numeric(FU2DUR)/30),
    Death=ifelse(v3_out_stat_e3_c3=="Alive",0,1)
  ) %>% 
  mutate_at(vars(v3_ascites, v3_varices,v3_lt,v3_variceal_bleeding,v3_hcc,v3_encephalopathy), ~replace_na(., 0)) %>% 
  mutate(across(.cols=c(v3_ascites, v3_varices,v3_lt,v3_variceal_bleeding,v3_hcc,v3_encephalopathy,Death),.fns=as.numeric),
         sex_e1_c1=as.factor(sex_e1_c1),
         v1_fibro_find_e1_c1=as.factor(v1_fibro_find_e1_c1))


## Decompensation as composite endpoint:
kms_decomp=survfit(Surv(DURN,v3_decomp_stat)~v3_ttt_resp_e3_c3_1,data=data3)

kms_decomp_p=ggsurvplot(kms_decomp,
                    pval=T,pval.method = T,pval.method.coord = c(0, 0.1),pval.coord = c(15, 0.1),
                    conf.int=F,risk.table = T,risk.table.col="strata",
                    cumevents = TRUE,
                    legend.labs=c("Not Responded","Responded"),
                    ncensor.plot=F,
                    surv.median.line = "hv",
                    xlab = "Month",
                    tables.theme = theme_cleantable())
kms_decomp_p

### All in one KM curve:
install.packages("Polychrome")
library(Polychrome)

# build-in color palette
library(RColorBrewer)
display.brewer.all()

fit <- list(Ascites= kms_ascites, Variceal_bleeding = kms_vb,HCC=kms_hcc,
            Encephalopathy=kms_encephal,LT=kms_LT,Death=kms_death)

ggsurvplot(fit,
           data=data3,combine = TRUE,
           risk.table = T,risk.table.col="strata",
           legend.labs=c("Asc:Not Responded","Asc:Responded","VB:Not Responded","VB:Responded",
                         "HCC:Not Responded","HCC:Responded","ENC:Not Responded","ENC:Responded",
                         "LT:Not Responded","LT:Responded","DT:Not Responded","DT:Responded"),
           xlab = "Day",
           tables.theme = theme_cleantable(),
           tables.height = 0.4,
           censor = FALSE,
           ncensor.plot=T,
           palette = "Paired",
           surv.median.line = "hv")

##### Logistic regressions for responders vs non responders #####
LRV2=data %>% 
  select(new_subject_id,sex_e1_c1,v1_pbc_dia_age_e1_c1,v1_fh_pbc_e1_c1,race_e1_c1,v1_Fib_recat,CPC_RECAT,
         v1_hri_class_e1_c1,v2_ttt_resp_e2_c2_1,v2_decomp_stat,AMA_RECAT,
         v1_lab1_e1_c1,v1_lab2_e1_c1,v1_lab3_e1_c1,v1_lab4_e1_c1,v1_lab5_e1_c1,v1_lab6_e1_c1,v1_lab7_e1_c1,v1_lab8_e1_c1,v1_lab12_e1_c1,v1_lab13_e1_c1,v1_lab14_e1_c1) %>% 
  mutate(
         Resp=ifelse(v2_ttt_resp_e2_c2_1=="Not Responded",0,1)) %>% 
  mutate(across(.cols = c(CPC_RECAT,v1_hri_class_e1_c1,v1_fh_pbc_e1_c1,v1_Fib_recat,Resp,sex_e1_c1,race_e1_c1,v2_ttt_resp_e2_c2_1),.fns=as.factor))


LRV2$CPC_RECAT=relevel(LRV2$CPC_RECAT,ref="Normal")
LRV2$v1_Fib_recat=relevel(LRV2$v1_Fib_recat,ref="NA")

uvreg=LRV2 %>% 
  select(v2_ttt_resp_e2_c2_1,v1_pbc_dia_age_e1_c1,sex_e1_c1,CPC_RECAT,
         v1_Fib_recat,v1_hri_class_e1_c1,v2_decomp_stat,AMA_RECAT,
         starts_with("v1_lab")) %>% 
  tbl_uvregression(
    method = glm,
    y=v2_ttt_resp_e2_c2_1,
    method.args = list(family=binomial),
    exponentiate = T
  ) %>% 
  add_global_p() %>% 
  bold_p(t=0.2)

uvreg

uni_reg=glm(v2_ttt_resp_e2_c2_1~v1_lab7_e1_c1,data = LRV2,family = binomial(link="logit"))
exp(uni_reg$coefficients)
exp(confint(uni_reg))

mreg=glm(v2_ttt_resp_e2_c2_1~v1_pbc_dia_age_e1_c1+sex_e1_c1+
           v1_lab1_e1_c1+v1_lab2_e1_c1+v1_lab3_e1_c1+v1_lab4_e1_c1+v1_lab5_e1_c1+v1_lab12_e1_c1+
           v1_Fib_recat+v1_hri_class_e1_c1+v2_decomp_stat
         
         ,data = LRV2,family = binomial(link="logit"))
#summary(mreg)
exp(mreg$coefficients)
exp(confint(mreg))

mreg1=tbl_regression(mreg,
                    exponentiate = T
) %>% 
  add_global_p()

mreg1

Merged_tables= tbl_merge(
  list(uvreg,mreg1),
  tab_spanner=c("**Univariable**","**Multivariable**")
)
Merged_tables


##### Cox regressions #####

data4=data3 %>% 
  mutate_at(c("v1_fh_pbc_e1_c1","race_e1_c1","v1_surg_e1_c1","v1_pbc_criteria_e1_c1","v1_hri_class_e1_c1","v1_fibro_find_e1_c1"),factor) %>% 
  mutate(CPC=ifelse(v1_pbc_cpc_e1_c1>=5 & v1_pbc_cpc_e1_c1<=6,"Class A",
                    ifelse(v1_pbc_cpc_e1_c1>6 & v1_pbc_cpc_e1_c1<=9,"Class B",
                           ifelse(v1_pbc_cpc_e1_c1>9 ,"Class C","Class"))),
         Resp=ifelse(v2_ttt_resp_e2_c2_1=="Not Responded",0,1)) %>% 
  mutate(across(.cols = c(v1_pbc_cpc_e1_c1,v1_fibro_find_e1_c1,v1_hri_class_e1_c1,v1_fh_pbc_e1_c1,CPC,Resp,sex_e1_c1,race_e1_c1,v1_pbc_cpc_e1_c1,v1_fibro_scr_e1_c1,v2_ttt_resp_e2_c2_1),.fns=as.factor))

cox_ascites=coxph(Surv(DURN,v3_ascites)~v1_pbc_dia_age_e1_c1+
                    v3_ttt_resp_e3_c3_1+v1_fibro_find_e1_c1+
                    v1_fh_pbc_e1_c1+
                    v1_lab1_e1_c1+v1_lab2_e1_c1+v1_lab3_e1_c1+v1_lab4_e1_c1+v1_lab5_e1_c1+
                    v1_lab6_e1_c1+v1_lab7_e1_c1+v1_lab8_e1_c1+
                    v1_lab12_e1_c1
                  
                    ,data = data4)


tbl_regression(cox_ascites,
               exponentiate = T
) %>% add_global_p()

#### CoX decompensated status at v3 ####
data_cox=data3 %>% 
  select(new_subject_id,v3_decomp_stat,DURN,v1_pbc_dia_age_e1_c1,sex_e1_c1,v1_Fib_recat,v1_hri_class_e1_c1,
         v1_lab1_e1_c1,v1_lab2_e1_c1,v1_lab3_e1_c1,v1_lab4_e1_c1,v1_lab6_e1_c1,v1_lab7_e1_c1,v1_lab5_e1_c1,v1_lab8_e1_c1,
         v1_lab12_e1_c1,v1_lab13_e1_c1,v1_lab14_e1_c1,
         CPC_RECAT,AMA_RECAT) %>% 
  filter(CPC_RECAT!="Class B") %>% 
  mutate(Age_recat=ifelse(v1_pbc_dia_age_e1_c1<40,"<40",">=40"))

data_cox$v1_Fib_recat=relevel(as.factor(data_cox$v1_Fib_recat),ref="NA")
data_cox$CPC_RECAT=relevel(as.factor(data_cox$CPC_RECAT),ref="Normal")


cox_decomp=coxph(Surv(DURN,v3_decomp_stat)~v1_pbc_dia_age_e1_c1+Age_recat+
                   v1_Fib_recat+CPC_RECAT+AMA_RECAT+
                    v1_lab1_e1_c1+v1_lab2_e1_c1+v1_lab3_e1_c1+v1_lab4_e1_c1+v1_lab5_e1_c1+
                   v1_lab6_e1_c1+v1_lab7_e1_c1+v1_lab8_e1_c1+v1_lab12_e1_c1+v1_lab13_e1_c1+v1_lab14_e1_c1
                  
                  ,data = data_cox)

extracted_data <- model.frame(cox_decomp) 
export(extracted_data,"extracted_data.xlsx")
descr=import("extracted_data.xlsx")
tbl_uvregression(
  data_cox[,-1],
  method = coxph,
  y = Surv(DURN, v3_decomp_stat),
  exponentiate = TRUE,
  pvalue_fun = function(x) style_pvalue(x, digits = 2)
) %>% 
  add_global_p()

descr %>% 
  select(Status,Age_recat,v1_Fib_recat,CPC_RECAT,AMA_RECAT) %>% 
  tbl_summary(
    by = Status,
    statistic = list(
      all_categorical() ~ "{n}({p}%)"
  )) %>% 
  add_overall()



#### CoX decompensated status at v2 ####
data_cox2=data3 %>% 
  select(new_subject_id,v2_decomp_stat,DURN2,v1_pbc_dia_age_e1_c1,sex_e1_c1,v1_Fib_recat,v1_hri_class_e1_c1,
         v1_lab1_e1_c1,v1_lab2_e1_c1,v1_lab3_e1_c1,v1_lab4_e1_c1,v1_lab6_e1_c1,v1_lab7_e1_c1,v1_lab5_e1_c1,v1_lab8_e1_c1,
         v1_lab12_e1_c1,v1_lab13_e1_c1,v1_lab14_e1_c1,
         CPC_RECAT,AMA_RECAT) %>% 
  filter(CPC_RECAT!="Class B") %>% 
  mutate(Age_recat=ifelse(v1_pbc_dia_age_e1_c1<40,"<40",">=40"))

data_cox2$v1_Fib_recat=relevel(as.factor(data_cox$v1_Fib_recat),ref="NA")
data_cox2$CPC_RECAT=relevel(as.factor(data_cox$CPC_RECAT),ref="Normal")


cox_decomp2=coxph(Surv(DURN2,v2_decomp_stat)~v1_pbc_dia_age_e1_c1+Age_recat+
                   v1_Fib_recat+CPC_RECAT+AMA_RECAT+
                   v1_lab1_e1_c1+v1_lab2_e1_c1+v1_lab3_e1_c1+v1_lab4_e1_c1+v1_lab5_e1_c1+
                   v1_lab6_e1_c1+v1_lab7_e1_c1+v1_lab8_e1_c1+v1_lab12_e1_c1+v1_lab13_e1_c1+v1_lab14_e1_c1
                 
                 ,data = data_cox2)

tbl_uvregression(
  data_cox2[,-1],
  method = coxph,
  y = Surv(DURN2, v2_decomp_stat),
  exponentiate = TRUE,
  pvalue_fun = function(x) style_pvalue(x, digits = 2)
) %>% 
  add_global_p()


extracted_data2 <- model.frame(cox_decomp2) 
export(extracted_data2,"extracted_data2.xlsx")
descr2=import("extracted_data_v2.xlsx")
descr2 %>% 
  select(Status,Age_recat,v1_Fib_recat,CPC_RECAT,AMA_RECAT) %>% 
  tbl_summary(
    by = Status,
    statistic = list(
      all_categorical() ~ "{n}({p}%)"
    )) %>% 
  add_overall()

##### Response based on total serum bilirubin 0.6 x ULN ####
Respv2=data %>% 
  select(new_subject_id,v2_ttt_resp_e2_c2_1,v2_lab4_e2_c2) %>% 
  mutate(Total_serum_bilirubin=ifelse(v2_lab4_e2_c2>=0.6,">=0.6","<0.6"))


Respv2 %>% 
  tbl_cross(
  row=Total_serum_bilirubin,
  col=v2_ttt_resp_e2_c2_1,
  percent = "row"
) %>% 
  add_p()

Respv3=data %>% 
  select(new_subject_id,v3_ttt_resp_e3_c3_1,v3_lab4_e3_c3) %>% 
  mutate(Total_serum_bilirubin=ifelse(v3_lab4_e3_c3>=0.6,">=0.6","<0.6"))

Respv3 %>% 
  tbl_cross(
    row=Total_serum_bilirubin,
    col=v3_ttt_resp_e3_c3_1,
    percent = "row",
    missing = "no"
  ) %>% 
  add_p()  

#####Figure: Lab Mean(SD)/Medain IQR:
data5=data %>% 
  select(new_subject_id,101,104,107,110,113,116,119,122,134,137,140,
         198,201,203,206,208,211,214,217,313,316,319,322,325,328,331,334) %>% 
  pivot_longer(!new_subject_id,names_to = "Lab_Test",values_to = "Value") %>% 
  mutate(Visit=case_when(
    str_detect(Lab_Test,"v1_lab")~"At Baseline",
    str_detect(Lab_Test,"v2_lab")~"At 1 year follow-up",
    str_detect(Lab_Test,"v3_lab")~"End of follow-up"
    )
    ) %>% 
  mutate(Lab_Test=case_when(
    Lab_Test=="v1_lab1_e1_c1"~"ALT(U/L)",
    Lab_Test=="v1_lab2_e1_c1"~"AST(U/L)",
    Lab_Test=="v1_lab3_e1_c1"~"ALP(U/L)",
    Lab_Test=="v1_lab4_e1_c1"~"GGT(U/L)",
    Lab_Test=="v1_lab5_e1_c1"~"Total Bilirubin(mmol/L)",
    Lab_Test=="v1_lab6_e1_c1"~"Total Cholesterol(mmol/L)",
    Lab_Test=="v1_lab7_e1_c1"~"Platelets(10^9/L)",
    Lab_Test=="v1_lab8_e1_c1"~"Albumin(g/L)",
    Lab_Test=="v1_lab12_e1_c1"~"Prothrombin ratio(Seconds)",
    Lab_Test=="v1_lab13_e1_c1"~"IgM(g/L)",
    Lab_Test=="v1_lab14_e1_c1"~"IgG(g/L)",
    Lab_Test=="v2_lab1_e2_c2"~"ALT(U/L)",
    Lab_Test=="v2_lab2_e2_c2"~"AST(U/L)",
    Lab_Test=="v2_lab3_e2_c2"~"ALP(U/L)",
    Lab_Test=="v2_lab4_e2_c2"~"Total Bilirubin(mmol/L)",
    Lab_Test=="v2_lab5_e2_c2"~"Total Cholesterol(mmol/L)",
    Lab_Test=="v2_lab6_e2_c2"~"Platelets(10^9/L)",
    Lab_Test=="v2_lab7_e2_c2"~"Albumin(g/L)",
    Lab_Test=="v2_lab8_e2_c2"~"Prothrombin ratio(Seconds)",
    Lab_Test=="v3_lab1_e3_c3"~"ALT(U/L)",
    Lab_Test=="v3_lab2_e3_c3"~"AST(U/L)",
    Lab_Test=="v3_lab3_e3_c3"~"ALP(U/L)",
    Lab_Test=="v3_lab4_e3_c3"~"Total Bilirubin(mmol/L)",
    Lab_Test=="v3_lab5_e3_c3"~"Total Cholesterol(mmol/L)",
    Lab_Test=="v3_lab6_e3_c3"~"Platelets(10^9/L)",
    Lab_Test=="v3_lab7_e3_c3"~"Albumin(g/L)",
    Lab_Test=="v3_lab8_e3_c3"~"Prothrombin ratio(Seconds)"
  ))

#IQR as ranges:
data5 %>%
  filter(Lab_Test=="ALP(U/L)") %>% 
  select(Lab_Test,Value,Visit) %>% 
  tbl_summary(
    by=Visit,
    statistic =Value~"{median} ({p25}, {p75})",
    digits = all_continuous() ~ 0
  )

# Fibrosis score v1,v2,v3:
datafbr=data %>% 
  select(new_subject_id,v1_Fib_recat,v2_Fib_recat,v3_Fib_recat)

datafbr %>%
  select(-new_subject_id) %>% 
  tbl_summary(
  )

# plot1:
data5 %>% 
  ggboxplot(x="Lab_Test",y="Value",
            add = "mean_sd",facet.by = "Visit")+
  rotate_x_text(45)

#plot2:
data6=data5 %>% 
  group_by(Visit,Lab_Test) %>% 
  summarise(mean_lab=mean(Value,na.rm = T),
            sd_lab=sd(Value,na.rm=T),
            median=median(Value,na.rm=T),
            Q1=quantile(Value,0.25,na.rm=T),
            Q3=quantile(Value,0.75,na.rm=T)
            ) %>% 
  mutate(Visit=factor(Visit,levels=c("At Baseline","At 1 year follow-up","End of follow-up")),
         Lab_test=factor(Lab_Test,levels=c("ALT(U/L)","AST(U/L)","ALP(U/L)","GGT(U/L)","Total Bilirubin(mmol/L)","Total Cholesterol(mmol/L)","Platelets(10^9/L)","Albumin(g/L)","Prothrombin ratio(Seconds)","IgM(g/L)","IgG(g/L)")))

data6$Lab_Test=factor(data6$Lab_Test,levels=c("ALT(U/L)","AST(U/L)","ALP(U/L)","GGT(U/L)","Total Bilirubin(mmol/L)","Total Cholesterol(mmol/L)","Platelets(10^9/L)","Albumin(g/L)","Prothrombin ratio(Seconds)","IgM(g/L)","IgG(g/L)"))

#Mean±SD plot: 
data6 %>% 
  ggplot(aes(x=Lab_Test, y=mean_lab)) +  
  geom_point()+ 
  geom_errorbar(aes(ymin=mean_lab-sd_lab, ymax=mean_lab+sd_lab), width=.2, 
                position=position_dodge(0.05)) +
  facet_wrap(~Visit)+
  labs(x="Lab paraameter",y="Mean±SD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Median (IQR):
data6 %>% 
  ggplot(aes(x = Lab_Test)) +
  geom_pointrange(aes(y = median,
                      ymin = Q1,
                      ymax = Q3),
                  show.legend = FALSE) +
  labs(y = 'Median (IQR)',x="Lab paraameter")+
  facet_wrap(~Visit)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

data5wide=data %>% 
  select(new_subject_id,101,104,107,110,113,116,119,122,134,137,140,
         198,201,203,206,208,211,214,217,313,316,319,322,325,328,331,334)

export(data5wide,"FR.xlsx")
