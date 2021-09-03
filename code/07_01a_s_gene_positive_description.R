##########################################################
# Name of file: 07a_s_gene_dropout_description.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@nhs.net
# Original date: 06 August 2020
# Latest update author (if not using version control) - Chris Robertson chrisobertson@nhs.net
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort and merges in the risk groups 
#                         selects out the positive cases only and fits a prediction model
#                         reads in new positives and forecasts 28 day deaths + Hospitalisations
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
library(lubridate)
#Load data

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop
project_path <- paste0(Location,"EAVE/GPanalysis/progs/CR")

EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-07-28.rds"))
EAVE_cohort <- filter(EAVE_cohort, !duplicated(EAVE_LINKNO))

table(EAVE_cohort$death_covid, is.na(EAVE_cohort$NRS.Date.Death), exclude=NULL)
table(EAVE_cohort$icu_death, is.na(EAVE_cohort$date_icu_death), exclude=NULL)
table(EAVE_cohort$hosp_covid, is.na(EAVE_cohort$date_hosp_covid), exclude=NULL)

#z <- filter(EAVE_cohort, (icu_death==1 & is.na(date_icu_death)) )
#correct errors
EAVE_cohort <- EAVE_cohort %>% 
  mutate(death_covid = if_else(death_covid==1 & is.na(NRS.Date.Death), 0,death_covid)) %>% 
  mutate(icu_death = if_else((icu_death==1) & is.na(date_icu_death), 0, icu_death))

a_begin <- as.Date("2021-04-01")
#remove all who have died before the beginning
EAVE_cohort <- filter(EAVE_cohort, is.na(NRS.Date.Death) | (!is.na(NRS.Date.Death) & NRS.Date.Death > a_begin))

EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))
EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)

#adjust inconsistencies in the endpoints and times - all hosp have an admission date
z_max_date_death <- max(EAVE_cohort$NRS.Date.Death, na.rm=T)
z_max_date_icu <- max(EAVE_cohort$date_icu_death, na.rm=T)
#EAVE_cohort <- EAVE_cohort %>% mutate(NRS.Date.Death = case_when(death_covid==1 & is.na(NRS.Date.Death) ~ SpecimenDate + 21,
#                                                       TRUE ~ NRS.Date.Death),
#                            date_icu_death = case_when(icu_death==1 & is.na(date_icu_death) ~ SpecimenDate + 14,
#                                                       TRUE ~ date_icu_death ) ) %>% 
#  mutate(NRS.Date.Death = case_when(NRS.Date.Death > z_max_date_death  ~ z_max_date_death,
#                                    TRUE ~ NRS.Date.Death),
#         date_icu_death = case_when(date_icu_death > z_max_date_icu ~ z_max_date_icu,
#                                    TRUE ~ date_icu_death ) )
EAVE_cohort <- EAVE_cohort %>% mutate(death_covid = case_when(death_covid==1 & is.na(NRS.Date.Death) ~ 0,
                                                              TRUE ~ death_covid),
                                      icu_death = case_when(icu_death==1 & is.na(date_icu_death) ~ 0,
                                                            TRUE ~ icu_death ) )

#z <- readRDS(paste0(Location,"EAVE/GPanalysis/data/ECOSS_deduped_linked.rds"))

rg <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE.rds"))
rg <- filter(rg, !duplicated(EAVE_LINKNO))
rg <- rg %>% dplyr::select(-EAVE_CHRONIC_PANCREATITIS, -EAVE_PREGNANCY, -EAVE_TRANSPLANT)
z <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds"))
z <- filter(z, !duplicated(EAVE_LINKNO))
z <- z %>% dplyr::select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>% 
  dplyr::rename(EAVE_Smoke = EAVE_Smoking_Status_Worst)

rg <- rg %>% left_join(z, by="EAVE_LINKNO")

#read in the Previous Tests data
cdw_full  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/CDW_full.rds"))
cdw_full <- cdw_full %>% mutate(date_ecoss_specimen = as_date(date_ecoss_specimen))
cdw_full <- cdw_full %>% #dplyr::select(EAVE_LINKNO, date_ecoss_specimen, test_result) %>% 
  #  filter(subject_upi %in% z_id) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen, desc(test_result)) %>%
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen)))  #get one test per person per day - preferentially positive test 
  
summary(cdw_full)
z <- cdw_full %>% filter(date_ecoss_specimen < a_begin) %>% 
  group_by(EAVE_LINKNO) %>% dplyr::summarise(n_tests=n()) %>% 
  dplyr::select(EAVE_LINKNO, n_tests)
rg <- rg %>% left_join(z, by="EAVE_LINKNO")
rg <- rg %>% mutate(n_tests = if_else(is.na(n_tests), 0L, n_tests))

#read in the GP vaccination data
source("00_Read_GP_Vaccinations.R")

Positive_Tests <- cdw_full %>% filter(test_result=="POSITIVE") %>% 
  dplyr::select(EAVE_LINKNO, date_ecoss_specimen) %>% dplyr::rename(SpecimenDate=date_ecoss_specimen)
summary(Positive_Tests)

#use the EAVE hospitalisations
z <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO, SpecimenDate, hosp_covid, date_hosp_covid, NRS.Date.Death) %>% 
  filter(hosp_covid==1) %>% 
  filter(date_hosp_covid > a_begin) %>% 
  dplyr::rename(admission_date = date_hosp_covid) %>% 
  mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>% 
  dplyr::select(-hosp_covid, -NRS.Date.Death)
covid_hospitalisations <- z

#update with more timely data
z <- readRDS(paste0(Location,"EAVE/GPanalysis/data/covid_hospitalisations.rds"))
summary(z)
#z <- z %>% filter(!is.na(admission_date))
z1 <- z %>% filter(admission_date >= a_begin) %>% 
  filter(!(EAVE_LINKNO %in% covid_hospitalisations$EAVE_LINKNO)) # find id's not in current list
z1 <- z1 %>% left_join(Positive_Tests, by="EAVE_LINKNO") %>% 
#  dplyr::rename(SpecimenDate=specimen_date) %>% 
  dplyr::select("EAVE_LINKNO","SpecimenDate","admission_date")
covid_hospitalisations <- bind_rows(covid_hospitalisations, z1) %>% 
  filter(EAVE_LINKNO %in% EAVE_cohort$EAVE_LINKNO)

z <- covid_hospitalisations %>% group_by(admission_date) %>% dplyr::summarise(N=n())
z %>% ggplot(aes(x=admission_date, y=N)) + geom_point()
summary(covid_hospitalisations)

#The dates of these two below depend on the updated endpoint linkages
#use the EAVE severe cases
z <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO, SpecimenDate, icu_death, date_icu_death, NRS.Date.Death) %>% 
  filter(icu_death==1) %>% 
  filter(date_icu_death > a_begin) %>% 
  dplyr::rename(admission_date = date_icu_death) %>% 
  mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>% 
  dplyr::select(-icu_death, -NRS.Date.Death)
covid_icu_death <- z
summary(covid_icu_death)

#use the EAVE death cases
z <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO, SpecimenDate, death_covid, NRS.Date.Death) %>% 
  filter(death_covid==1) %>% 
  filter(NRS.Date.Death > a_begin) %>% 
  dplyr::rename(admission_date = NRS.Date.Death) %>% 
  dplyr::select(-death_covid)
covid_death <- z
summary(covid_death)
#use the EAVE death cases
#z <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO, SpecimenDate, death_covid, NRS.Date.Death) %>% 
#  filter(!is.na(NRS.Date.Death)) %>% 
#  filter(NRS.Date.Death > a_begin) %>% 
#  dplyr::rename(admission_date = NRS.Date.Death) %>% 
#  dplyr::select(-death_covid)
#any_death <- z

all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))
summary(all_deaths)
#get covid death certificate deaths
z <- all_deaths %>%  
  mutate(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~if_else(. %in% c("U071","U072"), 1,0)))
z <- z %>% rowwise() %>% mutate(rowsum = sum(c_across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9))) %>% 
  mutate(covid_death_cert = if_else(rowsum>=1,1,0)) %>% 
  dplyr::select(EAVE_LINKNO, NRS.Date.Death, covid_death_cert)
all_deaths_covid_dth_cert <- z


all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/any_hospitalisation_post_01022020.rds"))
summary(all_hospitalisations)
#all_hospitalisations <- left_join(all_hospitalisations, covid_hospitalisations, by="EAVE_LINKNO", suffix=c("","_covid")) %>% 
#  filter(is.na(admission_date_covid) | !is.na(admission_date_covid)&(admission_date_covid > admission_date + 14) ) %>% 
#  dplyr::select(-SpecimenDate, - admission_date_covid)

#cohort + risk groups
df_cohort <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO:ur6_2016_name, age_gp, eave_weight) %>% 
  left_join(rg, by="EAVE_LINKNO")

z <- Positive_Tests %>%  mutate(days = as.numeric(SpecimenDate - a_begin)) %>% 
  mutate(test_before_begin = cut(days, breaks = c((min(days)-1), -28, -21, -14, -7, 0, max(days)),
                                labels=c("1+m", "4w","3w","2w","0-6d","post-start")))
df_cohort <- df_cohort %>% left_join(dplyr::select(z, EAVE_LINKNO, test_before_begin), by="EAVE_LINKNO")
df_cohort <- df_cohort %>% mutate(test_before_begin = as.character(test_before_begin)) %>% 
  mutate(test_before_begin = if_else(is.na(test_before_begin), "no pos test",test_before_begin) )

#z <- readRDS(paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid.rds"))
#z <- z %>% dplyr::select(-(Sex:age_gp), -Q_BMI)
#z <- filter(z, !duplicated(EAVE_LINKNO))
#z1 <- df_cohort %>% dplyr::select(-(EAVE_ASTHMA:EAVE_HYPERTENSION), -EAVE_Smoke, -EAVE_BP, -n_risk_gps, -eave_weight) %>% 
#  left_join(z, by="EAVE_LINKNO")
#z1 <- z1 %>% mutate(n_tests = if_else(is.na(n_tests),0L,n_tests) )
#df_cohort <- filter(z1, !is.na(eave_weight)) #omit any who - need to fix -  do not match - omits under 18's
#df_cohort$n_tests_gp <- cut(df_cohort$n_tests, breaks = c(-1,0,1,2,3,9,100), labels=c("0","1","2","3","4-9","10+"))

df_cohort <- df_cohort %>% left_join(Vaccinations, by="EAVE_LINKNO") 
df_cohort <- df_cohort %>% mutate(flag_incon = if_else(is.na(flag_incon), 0,flag_incon))

z_ids <- c(Vaccinations$EAVE_LINKNO, all_deaths$EAVE_LINKNO, covid_hospitalisations$EAVE_LINKNO, 
           filter(EAVE_cohort, tested==1)$EAVE_LINKNO, all_hospitalisations$EAVE_LINKNO) %>% unique()
#summary(filter(EAVE_cohort, !(EAVE_LINKNO %in% z_ids))$eave_weight)
z_N <- round(sum(df_cohort$eave_weight) )
z_k <- sum(df_cohort$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(df_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))
z <- df_cohort %>% mutate(ew = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )
df_cohort <- z %>% dplyr::select(-eave_weight) %>% dplyr::rename(eave_weight=ew)

z <- read_csv(paste0(Location,"/EAVE/GPanalysis/data/restored/map_files/Datazone2011Lookup.csv")) %>% 
  dplyr::select(DataZone, InterZone, Council, HB)
df_cohort <- df_cohort %>% left_join(z, by="DataZone") %>% 
  mutate(HB = if_else(is.na(HB),"Unknown", HB),
         InterZone = if_else(is.na(InterZone),"Unknown", InterZone),
         Council = if_else(is.na(Council),"Unknown", Council))

#saveRDS(df_cohort, paste0(project_path,"/output/temp/df_cohort.rds"))
#df_cohort <- readRDS(paste0(project_path,"/output/temp/df_cohort.rds"))

s_gene <- readRDS(paste0(Location,"EAVE/GPanalysis/data/cases_s_gene.rds")) %>% 
  filter(SpecimenDate >= a_begin) %>% 
  arrange(EAVE_LINKNO, SpecimenDate) %>% 
  filter(!duplicated(EAVE_LINKNO))
summary(s_gene)



wgs <- readRDS(paste0(Location,"EAVE/GPanalysis/data/WGS_latest.rds")) %>% 
  mutate_at(c("Collection_Date","Sequencing_Date","Alignment_Date"), ~ as.Date(. , format="%d/%m/%Y")) %>% 
  filter(Collection_Date >= a_begin)

#z_id <- filter(wgs, duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO)
#z <- filter(wgs, EAVE_LINKNO %in% z_id) %>% arrange(EAVE_LINKNO)
#duplicates all have the same lineage
wgs <- wgs %>% arrange(EAVE_LINKNO, Collection_Date) %>% filter(!duplicated(EAVE_LINKNO))
a_end_wgs <- max(wgs$Collection_Date) - 2  #table(wgs$Collection_Date)  #check each time


#below not needed for nhs and lh data
s_df <- s_gene %>%  inner_join(df_cohort, by="EAVE_LINKNO")
s_df <- s_df %>%  mutate(ageYear = ifelse(ageYear >= 100, 100, ageYear)) 
s_df <- s_df %>% mutate(days = as.numeric(SpecimenDate -min(SpecimenDate)))
#s_df <- s_df %>% mutate(bmi_gp = cut(bmi_impute, c(12, 19, 24,29,34,39,100)))
#z_vars <- c("Q_HOME_CAT","Q_LEARN_CAT","Q_DIAG_CKD_LEVEL","EAVE_BP","EAVE_Smoke")
#s_df <- s_df %>% mutate_at(all_of(z_vars) , ~ as.factor(.))

s_df <- s_df %>%  mutate(vs1 = case_when(is.na(date_vacc_1) | date_vacc_1 > SpecimenDate ~ "uv",
                                         date_vacc_1 <= SpecimenDate &  date_vacc_1 > SpecimenDate - 28 ~ "v1_0:27",
                                         TRUE ~ "v1_28+"),
                         vs2 = case_when(is.na(date_vacc_2) | date_vacc_2 > SpecimenDate ~ "uv",
                                         date_vacc_2 <= SpecimenDate &  date_vacc_2 > SpecimenDate - 14 ~ "v2_0:13",
                                         TRUE ~ "v2_14+")) %>% 
  mutate(vs = if_else(vs2=="uv", vs1,vs2))

s_df <- s_df %>% left_join(dplyr::select(wgs, EAVE_LINKNO, VariantofInterest), by="EAVE_LINKNO")
s_df <- s_df %>% mutate(variant = case_when(is.na(VariantofInterest) ~ "not_sequenced",
                                            VariantofInterest=="VOC-20DEC-01" ~ "alpha",
                                            VariantofInterest=="VOC-21APR-02" ~ "delta",
                                            TRUE ~ "other"))

saveRDS(s_df, paste0(project_path,"/SGene.rds"))
#s_df <- readRDS(paste0(project_path,"/SGene.rds"))


z_names <- names(s_df)
z_names <- z_names[!(z_names %in% c("EAVE_LINKNO","SpecimenDate","true_s_gene_dropout","ageYear","DataZone"))]

#z_names_eave <- z_names[grepl("^Q_", z_names)]
#z_names_eave <- c(z_names_eave, c("EAVE_BP","EAVE_Smoke", "bmi_gp", "n_tests_gp", "test_before_begin"))
#z_names_rest <- z_names[!(z_names %in% z_names_eave)]
#z_names_rest <- z_names_rest[17:23]
z_names_rest <- c("age_gp","Sex","simd2020_sc_quintile", "n_risk_gps", "vs")

for (z_var in z_names_rest) { 
#z_var <- "vs"
z1 <- s_df %>% group_by_at(c("true_s_gene_dropout", z_var) ) %>% dplyr::summarise(N=n())
z2 <- s_df %>% group_by_at(c("true_s_gene_dropout") ) %>% dplyr::summarise(Total=n())
z.df <- left_join(z1,z2) %>% 
  mutate(P=N/Total*100) %>% 
  dplyr::rename(Var = 2) %>% 
  filter(true_s_gene_dropout != "Unknown")

gg <- z.df %>% ggplot(aes(x=Var, y=P, fill=true_s_gene_dropout)) + geom_col(position="dodge") +
  labs(y="Percentage", x= z_var, fill="S Gene")

print(gg)
}

table(s_df$true_s_gene_dropout, s_df$Sex)
table(s_df$variant, s_df$Sex)
z <- filter(s_df, SpecimenDate <= a_end_wgs)
table(z$variant, z$true_s_gene_dropout, exclude = NULL)
sum(z$variant != "not_sequenced")/nrow(z)  #proportion sequenced to a_end_wgs

for (z_var in z_names_eave) { 
  #z_var <- "Q_DIAG_AF"
  z1 <- s_df %>% group_by_at(c("true_s_gene_dropout", z_var) ) %>% dplyr::summarise(N=n())
  z2 <- s_df %>% group_by_at(c("true_s_gene_dropout") ) %>% dplyr::summarise(Total=n())
  z.df <- left_join(z1,z2) %>% 
    mutate(P=N/Total*100) %>% 
    dplyr::rename(Var = 2) %>% 
   # filter(Var=="Yes")%>% 
    filter(true_s_gene_dropout != "Unknown") %>% 
    mutate(Var=factor(Var))
  
  gg <- z.df %>% ggplot(aes(y=P, x=true_s_gene_dropout, fill=Var)) + geom_col() +
    labs(y="Percentage", title = z_var, x="S Gene")
  
  print(gg)
}

#z <- dplyr::select(EAVE_cohort, EAVE_LINKNO, hosp_covid , date_hosp_covid, Time.To.Hosp,
#                   death_covid, NRS.Date.Death, Time.To.Death )
s_df <- s_df %>% 
  left_join(dplyr::select(covid_hospitalisations, EAVE_LINKNO, admission_date),  by="EAVE_LINKNO") %>% 
  dplyr::rename(date_hosp_covid=admission_date)
a_end <- max(s_df$date_hosp_covid, na.rm=T)

s_df <- s_df %>% 
  left_join(dplyr::select(all_deaths, EAVE_LINKNO, NRS.Date.Death),  by="EAVE_LINKNO") 


#need to modify endpoints for deaths
z_df <- s_df %>%  filter(SpecimenDate <= a_end)
z_df <- z_df %>% mutate(Time.To.Hosp = if_else(!is.na(date_hosp_covid), 
                                            as.numeric(date_hosp_covid - SpecimenDate), 
                                            as.numeric(a_end - SpecimenDate)) ,
                     Time.To.Death = if_else(!is.na(NRS.Date.Death), 
                                             as.numeric(NRS.Date.Death - SpecimenDate), 
                                             as.numeric(a_end - SpecimenDate)) ) %>% 
  mutate(Time.To.Hosp = if_else(!is.na(NRS.Date.Death) & NRS.Date.Death < a_end & is.na(date_hosp_covid), 
                                as.numeric(NRS.Date.Death - SpecimenDate),
                                Time.To.Hosp ))
z_df <- z_df %>% mutate(Time.To.Hosp = if_else(Time.To.Hosp<0, 0, Time.To.Hosp)) %>% 
  mutate(hosp_covid = if_else(!is.na(date_hosp_covid),1,0))
  
library(survival)

z_names <- c(z_names_eave,z_names_rest[6:7])

#Response var
z.rv <- "hosp_covid" 
z.rv.time <- "Time.To.Hosp" 
#z.rv <- "covid_adm" 
#z.rv.time <- "time_test_covid_adm"
prediction_vars_h <- c("EAVE_CARE_HOME","EAVE_CHRONIC_HEART_DIS",             
 "EAVE_CHRONIC_KIDNEY_DIS", "EAVE_CHRONIC_PANCREATITIS","EAVE_CHRONIC_RESP_DIS","EAVE_DEMENTIA" ,                     
 "EAVE_DIABETES","EAVE_HAEMAT_MALIGNANCY" , "EAVE_HYPERTENSION","EAVE_MS_DEGEN_DIS",                  
 "EAVE_STROKE_TIA" ,"EAVE_TRANSPLANT", "EAVE_ULCER_DIS" ,"EAVE_BP", "EAVE_Smoke")

if (z.rv=="covid_adm") prediction_vars <- prediction_vars_h
if (z.rv=="covid_dth") prediction_vars <- z_names

z.df <- z_df
#z.df <- filter(df, true_s_gene_dropout != "Unknown")
#z.df <- filter(df, true_s_gene_dropout %in%  c("Positive S Gene","True S Gene Dropout"))
z.df$true_s_gene_dropout <- factor(z.df$true_s_gene_dropout, 
                levels=c("True S Gene Dropout", "Positive S Gene", "Weak Positive"), 
                labels=c("S Gene Negative","S Gene Positive","Weak S Positive"))
#z.df <- z.df %>% filter(ageYear >= 65 & ageYear <= 129)
z.df <- z.df %>% mutate(days = days-min(days))
z.df <- z.df %>% mutate(Time.To.Hosp = if_else(Time.To.Hosp>=28,28, Time.To.Hosp))
#z.df <- z.df %>% mutate(vacc = case_when(vs %in% c("v1_28+","v2_0:13","v2_14+") ~ "v1_28+v2",
#                                         TRUE ~ "uv_v1_0:27"))
z.df <- z.df %>% mutate(vacc = case_when(vs %in% c("v1_28+","v2_0:13","v2_14+") ~ "v1_28+v2",
                                         TRUE ~ vs))
z.df <- z.df %>% dplyr::rename(s_gene = true_s_gene_dropout)

z.df <- z.df %>% mutate(hosp_covid = if_else(Time.To.Hosp >= 15,0, hosp_covid),
  Time.To.Hosp = if_else(Time.To.Hosp >= 15,15,Time.To.Hosp))
z.df <- z.df %>% mutate(vacc = case_when(vs %in% c("v1_28+","v2_0:13","v2_14+") ~ "v1_28+v2",
                                         TRUE ~ vs))
#z.df %>% group_by(vacc_type, vacc) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()
z.df <- z.df %>% mutate(vt = vacc) %>% 
  mutate(vt = case_when(vt %in% c("v1_0:27","v1_28+v2") ~ paste(vacc_type,vt, sep="_"),
                        TRUE ~ vt)) %>% 
  mutate(vt = fct_relevel(vt, "uv")) 


z.df <- mutate(z.df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   s_gene"))
z <- pyears(fmla.final, data=z.df, data.frame = TRUE)
z$data %>% mutate(rate.month=event/pyears/12)



#fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  pspline(ageYear) + n_risk_gps + true_s_gene_dropout"))
#fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   true_s_gene_dropout + strata(days)"))
#fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   strata(days) + Sex + simd2020_sc_quintile + n_risk_gps + age_gp*true_s_gene_dropout"))
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) +pspline(days) + Sex + simd2020_sc_quintile +
                               s_gene +  n_risk_gps +vs "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) +pspline(days) + Sex + simd2020_sc_quintile +
                               s_gene +  n_risk_gps +vt "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) + strata(days) + Sex + simd2020_sc_quintile +
                               s_gene +  n_risk_gps +vs "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) + s_gene + vacc +
                               s_gene*vacc +  n_risk_gps + pspline(days) + Sex + simd2020_sc_quintile "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) + s_gene + vacc +
                               s_gene*vt +  n_risk_gps + pspline(days) + Sex + simd2020_sc_quintile "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) + pspline(days) + Sex + simd2020_sc_quintile +
             n_risk_gps + s_gene +  s_gene:vacc  "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) + pspline(days) + Sex + simd2020_sc_quintile +
             n_risk_gps + s_gene +  s_gene:vt  "))
z.fit <- coxph(fmla.final , data=z.df, subset = !(vt %in% c("Mo_v1_0:27","Mo_v1_28+v2")))
summary(z.fit)
z.df %>% group_by(vt) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()
z.df %>% group_by(s_gene) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()
z.df %>% group_by(vt,s_gene) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()




#fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) +pspline(days) + Sex + simd2020_sc_quintile +
#                                 n_risk_gps +vs "))
z.fit <- coxph(fmla.final , data=z.df)
summary(z.fit)
drop1(z.fit, test="Chisq")
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  pspline(ageYear) +pspline(days)  + Sex + simd2020_sc_quintile + ",
                               " true_s_gene_dropout + strata(Council) + ",
                               paste(prediction_vars, collapse= "+")))
z.fit <- coxph(fmla.final , data=z.df)
summary(z.fit)

#z.sub <- MASS::stepAIC(z.fit, direction="backward", k=2)


z1 <- summary(z.fit)$coefficients[5:32,]
z2 <- summary(z.fit)$conf.int[25:52,]
z_coef_hosp <- cbind(z2[,c(1,3,4)], z1[,6])

z_names_h <- prediction_vars

ptemp <- termplot(z.fit, se=TRUE, plot=FALSE)
ageterm <- ptemp$ageYear  # this will be a data frame
center <- with(ageterm, y[x==30])
ytemp <- ageterm$y + outer(ageterm$se, c(0, -1.96, 1.96),'*')
matplot(ageterm$x, exp(ytemp - center), log='y',type='l', lty=c(1,2,2), col=1,
        xlab="Age at Positive Test", ylab="Hospitalisation Hazard Ratio")
abline(h=1, lty=2)

ageterm <- ptemp$days  # this will be a data frame
center <- with(ageterm, y[x==0])
ytemp <- ageterm$y + outer(ageterm$se, c(0, -1.96, 1.96),'*')
matplot(ageterm$x, exp(ytemp - center), log='y',type='l', lty=c(1,2,2), col=1,
        xlab="Days from April 01", ylab="Hospitalisation Hazard Ratio")
abline(h=1, lty=2)

#z.df <- filter(df, true_s_gene_dropout != "Unknown")
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  s_gene "))
z.survfit <- survfit(fmla.plot, data=z.df)
plot(z.survfit, fun="event", col=1:3, xlab="Days from Test to Hospital Admission",
     ylab="Risk")
legend("topleft",col=1:3, lty=1, legend=levels(z.df$s_gene))

z_df <- df %>% filter(SpecimenDate<= Sys.Date()-15) %>% 
  filter(true_s_gene_dropout != "Unknown") %>% 
  select_at(c("age_gp", "Sex", "simd2020_sc_quintile", "true_s_gene_dropout",prediction_vars_h, "n_risk_gps", "covid_adm")) %>% 
  pivot_longer(cols=age_gp:n_risk_gps)

z_var <- "EAVE_Smoking_Status_Worst"

fun_1 <- function(df, z_resp) { 
  df %>% group_by(value) %>% dplyr::summarise(N=n(), Percent=mean(get(z_resp))*100) %>% 
  as.data.frame()
}

fun_1(filter(z_df, name=="Sex"), z_resp="covid_adm")
z.out <- plyr::ddply(z_df, .(name), fun_1, z_resp="covid_adm")

fun_2 <- function(df, z_resp, z_title) { 
z <-  df %>% group_by(name, value) %>% dplyr::summarise(N=n(), Percent=mean(get(z_resp))*100) %>% 
    ungroup() %>% as.data.frame()
z_name <- unique(z$name)
g <- z %>% ggplot(aes(x=value, y=Percent)) + geom_col() +
  labs(x=z_name, title=z_title)
print(g)
}

fun_2(filter(z_df, name=="Sex"), z_resp="covid_adm" ,z_title = "Hospitalisation")
plyr::d_ply(z_df, .(name), fun_2, z_resp="covid_adm",z_title = "Hospitalisation")

#logistic model hospital
z_df <- df %>% filter(SpecimenDate<= Sys.Date()-15) 
gam.model <- as.formula(paste(z.rv,"~  s(ageYear) + s(days) + Sex + simd2020_sc_quintile + ",
                               " true_s_gene_dropout + ",
                               paste(prediction_vars, collapse= "+")))
z.gam <- mgcv::gam(gam.model , data=z_df, family=binomial)
summary(z.gam)

z.gam <- mgcv::gam(covid_adm ~ s(ageYear) + s(days) + Sex + simd2020_sc_quintile + n_risk_gps + true_s_gene_dropout, data=z_df, family=binomial)
#z.gam <- mgcv::gam(covid_adm ~ true_s_gene_dropout, data=z_df, family=binomial)
summary(z.gam)
z1 <- table(z_df$true_s_gene_dropout, z_df$covid_adm)
z2 <- prop.table(z1 <- table(z_df$true_s_gene_dropout, z_df$covid_adm),1)
cbind(z1,z2)

#proportions
z1 <- table(z_df$true_s_gene_dropout, z_df$covid_dth)
z2 <- prop.table(z1 <- table(z_df$true_s_gene_dropout, z_df$covid_dth),1)
cbind(z1,z2)

#logistic model deaths
z_df <- df %>% filter(SpecimenDate<= Sys.Date()-29) 
gam.model <- as.formula(paste("covid_dth ~  s(ageYear) + s(days) + Sex + simd2020_sc_quintile + ",
                              " true_s_gene_dropout + n_risk_gps"))
z.gam <- mgcv::gam(gam.model , data=z_df, family=binomial)
summary(z.gam)


#cumulative risk plot
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  true_s_gene_dropout "))
z.survfit <- survfit(fmla.plot, data=z_df, subset=true_s_gene_dropout != "Unknown")
plot(z.survfit, fun="event", col=1:3, xlab="Days from Test to Death",
     ylab="Risk")
legend("topleft",col=1:3, lty=1, legend=levels(z.df$true_s_gene_dropout))


#trends in admissions among the different groups
#all have at least 14 days follow up
z_df <- df %>% filter(SpecimenDate < max(SpecimenDate)-14) %>% 
  group_by(SpecimenDate, true_s_gene_dropout) %>% 
  dplyr::summarise(N=n(), n_covid_adm = sum(covid_adm)) %>% ungroup() %>% 
  mutate(p_admit = n_covid_adm/N)
z_df %>% ggplot(aes(x=SpecimenDate, y=p_admit, colour=true_s_gene_dropout)) +geom_smooth()+
  labs(x="Date Sampled Taken", y= "Proportion Admitted", 
       title = "Admissions to hospital within 14 days of a positive test")

#trends
z_df <- df %>% filter(SpecimenDate < max(SpecimenDate)-2) %>% 
  filter(true_s_gene_dropout %in% c("True S Gene Dropout","Positive S Gene","Unknown_NHS")) %>% 
  group_by(SpecimenDate, true_s_gene_dropout) %>% 
  dplyr::summarise(N=n()) %>% ungroup()
g1 <- z_df %>% ggplot(aes(x=SpecimenDate, y=N, colour=true_s_gene_dropout)) +geom_point() +
  geom_smooth(se=FALSE, span=0.25) + geom_vline(xintercept=as.Date("2020-12-16"), linetype="dashed") +
  geom_vline(xintercept=as.Date("2020-12-31"), linetype="dashed") + geom_vline(xintercept=as.Date("2021-01-11"), linetype="dashed")+
  labs(x="Date Sampled Taken", y= "Number Positive", 
       title = "Positive cases per day")
print(g1)


#changes over time
z_df <- df %>% mutate(days_gp = cut(days, breaks=c(-1, 27, 55, 77), 
                                    labels=c("P1","P2","P3")),
                      age_gp2 = cut(ageYear, breaks=c(-1,49,64, 74, 110),
                                    labels=c("0-49","50-64","65-74","75+")))
z.rv <- "covid_dth" 
z.rv.time <- "time_test_covid_dth" 
#z.rv <- "covid_adm" 
#z.rv.time <- "time_test_covid_adm" 
z.df <- filter(z_df, true_s_gene_dropout %in% c("Positive S Gene","True S Gene Dropout"))
z.df$true_s_gene_dropout <- factor(z.df$true_s_gene_dropout)

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  days_gp + pspline(ageYear) + Sex +
        simd2020_sc_quintile + n_risk_gps + true_s_gene_dropout:days_gp"))
z.fit <- coxph(fmla.final , data=z.df)
#summary(z.fit)
#drop1(z.fit,test="Chisq")

#time period model
z2 <- summary(z.fit)$conf.int[c(25:27),]
z_coef_hosp <- z2[,c(1,3,4)]
dimnames(z_coef_hosp)[[2]] <- c("HR","LCL_95","UCL_95")
dimnames(z_coef_hosp)[[1]] <- c("Nov16-Dec13", "Dec14-Jan10", "Jan11-")
z_hr_period <- z_coef_hosp

##age gp
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  strata(days_gp) + age_gp + Sex +
        simd2020_sc_quintile + n_risk_gps + true_s_gene_dropout:age_gp2"))
z.fit <- coxph(fmla.final , data=z.df)
#summary(z.fit)
#drop1(z.fit,test="Chisq")

#age gp model
z2 <- summary(z.fit)$conf.int[c(29,31,33,35),]
z_coef_hosp <- 1/z2[,c(1,4,3)]
dimnames(z_coef_hosp)[[2]] <- c("HR","LCL_95","UCL_95")
dimnames(z_coef_hosp)[[1]] <- paste("Age",levels(z.df$age_gp2))
z_hr_age <- z_coef_hosp


z.df %>% group_by(age_gp2, true_s_gene_dropout) %>% 
  dplyr::summarise(N=n(), R=sum(covid_dth), days=sum(time_test_covid_dth))


#########################################################
# variant


fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   variant"))
z <- pyears(fmla.final, data=z.df, data.frame = TRUE, subset= variant %in% c("alpha","delta"))
z$data %>% mutate(rate.month=event/pyears/12)

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) +pspline(days) + Sex + simd2020_sc_quintile +
                               variant +  n_risk_gps +vs "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) + s_gene + vacc +
                               variant*vacc +  n_risk_gps "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(ageYear) + pspline(days) + Sex + simd2020_sc_quintile +
                               n_risk_gps + variant +  variant:vacc  "))


#z.fit <- coxph(fmla.final , data=z.df, subset=true_s_gene_dropout != "Weak Positive")
z.fit <- coxph(fmla.final , data=z.df, subset= variant %in% c("alpha","delta"))
summary(z.fit)
drop1(z.fit, test="Chisq")
