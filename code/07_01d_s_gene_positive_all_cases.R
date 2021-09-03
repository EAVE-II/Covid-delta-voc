##########################################################
# Name of file: 07_01d_s_gene_positive_all_cases.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@nhs.net
# Original date: 06 August 2020
# Latest update author (if not using version control) - Chris Robertson chrisobertson@nhs.net
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: sets up analysis for all positive tests
#                         run 07_01a_s_gene_positive_description.R
#                         to read in the data - line 291
#                         then get first positive test post 01 April
# Approximate run time: Unknown
##########################################################

#Need to get Positive Tests going right back

z_df <- cdw_full %>% filter(test_result == "POSITIVE") %>% filter(date_ecoss_specimen >= a_begin)
#get the first positive test in the period
z_df <- z_df %>% arrange(EAVE_LINKNO, date_ecoss_specimen) %>% 
  filter(!duplicated(EAVE_LINKNO))
#Ecoss is the nhs labs
z_df <- z_df %>% mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh") )
#  dplyr::select(EAVE_LINKNO, ecossid, sex, age_year, specimen_date, lab)
z_df <- z_df %>%  
  dplyr::select(EAVE_LINKNO, subject_sex, age, date_ecoss_specimen, lab, test_result_record_source) %>% 
  dplyr::rename(specimen_date=date_ecoss_specimen, age_year=age)

#get the first positive test for anyone before a_begin
z <- EAVE_cohort %>% filter(result==1) %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate) %>% 
  arrange(EAVE_LINKNO, SpecimenDate) %>% 
  filter(!duplicated(EAVE_LINKNO)) %>% 
  filter(SpecimenDate < a_begin)

z_df <- z_df %>% left_join(z, by="EAVE_LINKNO")
z_df <- z_df %>% mutate(days_from_first_pos = as.numeric(specimen_date - SpecimenDate)) %>% 
  dplyr::select(-SpecimenDate)

z <- s_gene %>% dplyr::select(EAVE_LINKNO, SpecimenDate, true_s_gene_dropout)
z_df <- z_df %>% left_join(z, by="EAVE_LINKNO")
#table(z_df$specimen_date == z_df$SpecimenDate, exclude=NULL) #virtually all the same
z_df <- z_df %>% mutate(true_s_gene_dropout = if_else(specimen_date != SpecimenDate, "unknown", true_s_gene_dropout))
z_df <- z_df %>%  mutate(s_gene = if_else(is.na(true_s_gene_dropout), "unknown", true_s_gene_dropout))
z_df <- z_df %>% mutate(s_gene = factor(s_gene, 
                  levels=c("True S Gene Dropout","Positive S Gene","Weak Positive", "unknown"),
                  labels=c("S_Neg","S_Pos","Weak_S_Pos","Unknown")))
z_df$true_s_gene_dropout <- NULL
z_df$SpecimenDate <- NULL

#link in hospitalisations - use all as some will be in hospital at time of test
#keep all admissions post a_begin and those in hospital at a_begin
z_h <- all_hospitalisations %>% dplyr::select(-validchi) %>% 
  filter(!(!is.na(discharge_date) & discharge_date <= a_begin))
#summary(filter(z, admission_date < a_begin))

z <- z_df %>% left_join(z_h, by="EAVE_LINKNO")
#discharges before a_begin
z <- z %>% mutate(ad=admission_date, dd=discharge_date, em = emergency) %>% 
  mutate(ad = as.Date(ifelse(!is.na(dd) & dd <= specimen_date, NA, ad), origin="1970-01-01"),
         em = ifelse(!is.na(dd) & dd <= specimen_date, NA, em), 
         dd = as.Date(ifelse(!is.na(dd) & dd <= specimen_date, NA, dd), origin="1970-01-01") )
#multiple admissions
z_id <- z %>% filter(duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO) %>% unique()
z_01 <- z %>% filter(!(EAVE_LINKNO %in% z_id)) # 0 or 1 admission
z_m <- z %>% filter((EAVE_LINKNO %in% z_id))  #multiple admissions
z_m <- z_m %>%  arrange(EAVE_LINKNO, ad) %>% 
  filter(!duplicated(EAVE_LINKNO))  # NA on ad go to the end, pick the first admission post specimen

z <- bind_rows(z_01, z_m) %>% 
  dplyr::select(-admission_date, -discharge_date, -emergency) %>% 
  dplyr::rename(admission_date=ad, discharge_date=dd, emergency=em)

z_df <- z
a_end <- max(c(max(z_df$admission_date, na.rm=T), max(z_df$discharge_date, na.rm=T)))
z_df <- z_df %>% mutate(Time.To.Hosp = if_else(!is.na(admission_date), as.numeric(admission_date - specimen_date),
                                               as.numeric(a_end - specimen_date))) %>% 
  mutate(In_Hosp_At_Test = if_else(Time.To.Hosp <= -3,  1L, 0L ) ,
         hosp_covid = if_else(!is.na(admission_date) & Time.To.Hosp <= 14,  1L, 0L ),
         Time.To.Hosp = case_when(Time.To.Hosp < 0 ~ 0,
                                  Time.To.Hosp >= 15 ~ 15,
                                  TRUE ~ Time.To.Hosp))
z_df <- z_df %>% mutate(hosp_covid_emerg = if_else(!is.na(emergency) & emergency ,hosp_covid, 0L ) )
z_df <- z_df %>% mutate(days = as.numeric(specimen_date- min(specimen_date)) )

#add in the vaccinations and risk groups
z_df <- z_df %>%  
  left_join(dplyr::select(df_cohort, -Sex, -ageYear, -age_gp ), by="EAVE_LINKNO" )
z_df <- z_df %>% mutate(in_eave = if_else(is.na(simd2020_sc_quintile), 0L,1L))

z_df <- z_df %>%  mutate(vs1 = case_when(is.na(date_vacc_1) | date_vacc_1 > specimen_date ~ "uv",
                                         date_vacc_1 <= specimen_date &  date_vacc_1 > specimen_date - 28 ~ "v1_0:27",
                                         TRUE ~ "v1_28+"),
                         vs2 = case_when(is.na(date_vacc_2) | date_vacc_2 > specimen_date ~ "uv",
                                         date_vacc_2 <= specimen_date &  date_vacc_2 > specimen_date - 14 ~ "v2_0:13",
                                         TRUE ~ "v2_14+")) %>% 
  mutate(vs = if_else(vs2=="uv", vs1,vs2))

z_df <- z_df %>% mutate(age_gp = cut(age_year, breaks = c(-1, 15, 39, 59, 74, 120),
                                     labels=c("0-15","16-39","40-59","60-74","75+")))
z_df <- z_df %>% dplyr::rename(sex=subject_sex)
z_df <- z_df %>% filter(!is.na(age_year))

#add in icu/death  - icu dates depend upon the endpoints linkage
z <- covid_icu_death %>% dplyr::select(-SpecimenDate) %>% 
  dplyr::rename(date_icu_death = admission_date)
z <- z_df %>% left_join(z, by="EAVE_LINKNO") %>% 
  mutate(covid_icu_death = if_else(is.na(date_icu_death), 0L, 1L)) %>% 
  mutate(Time.To.ICU.Death = if_else(is.na(date_icu_death), as.numeric(a_end - specimen_date),
                                    as.numeric(date_icu_death - specimen_date))) %>% 
  mutate(Time.To.ICU.Death = if_else(Time.To.ICU.Death < 0, 0, Time.To.ICU.Death))
z_df <- z

#covid death derived from all deaths
z <- z_df %>% left_join(all_deaths_covid_dth_cert, by="EAVE_LINKNO") %>% 
  mutate(death = if_else(is.na(NRS.Date.Death), 0L, 1L)) %>% 
  mutate(Time.To.Death = if_else(is.na(NRS.Date.Death), as.numeric(a_end - specimen_date),
                                       as.numeric(NRS.Date.Death - specimen_date))) %>% 
  mutate(Time.To.Death = if_else(Time.To.Death < 0, 0, Time.To.Death)) %>% 
  mutate(covid_death = if_else(is.na(covid_death_cert), 0, covid_death_cert)) %>% 
  mutate(covid_death = if_else(!is.na(NRS.Date.Death) & (NRS.Date.Death - specimen_date <= 28) , 1, covid_death))
z_df <- z

z_df <- z_df %>% mutate(prev_pos = if_else(is.na(days_from_first_pos), "not_prev_pos","prev_pos"))

z_df <- z_df %>% mutate(sg_lab = case_when(s_gene != "Unknown" ~ as.character(s_gene),
                                        TRUE ~ paste(s_gene, lab, sep="_"))) %>% 
  mutate(sg_lab=factor(sg_lab, levels =c("S_Neg","S_Pos","Weak_S_Pos","Unknown_lh","Unknown_nhs")))

z.df <- z_df
#z.df <- z.df %>% mutate(vacc = case_when(vs %in% c("v1_28+","v2_0:13","v2_14+") ~ "v1_28+v2",
#                                         TRUE ~ "uv_v1_0:27"))
z.df <- z.df %>% mutate(vacc = case_when(vs %in% c("v1_28+","v2_0:13","v2_14+") ~ "v1_28+v2",
                                         TRUE ~ vs))
#z.df %>% group_by(vacc_type, vacc) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()
z.df <- z.df %>% mutate(vt = vacc) %>% 
  mutate(vt = case_when(vt %in% c("v1_0:27","v1_28+v2") ~ paste(vacc_type,vt, sep="_"),
                        TRUE ~ vt)) %>% 
  mutate(vt = fct_relevel(vt, "uv")) 

z.df <- mutate(z.df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))
z.df <- mutate(z.df, Time.To.Death = if_else(Time.To.Death==0,0.1,Time.To.Death))

z.df <- z.df %>% left_join(dplyr::select(wgs, EAVE_LINKNO, VariantofInterest), by="EAVE_LINKNO")
z.df <- z.df %>% mutate(variant = case_when(is.na(VariantofInterest) ~ "not_sequenced",
                                            VariantofInterest=="VOC-20DEC-01" ~ "alpha",
                                            VariantofInterest=="VOC-21APR-02" ~ "delta",
                                            TRUE ~ "other")) %>% 
  dplyr::select(-VariantofInterest)

z_df <- z.df

z.df <- z_df %>% mutate(vs_phe = as.character(vs)) %>% 
  mutate(vs_phe = case_when(vs_phe %in% c("v1_28+","v2_0:13") ~ "v1_28v2_13",
                            TRUE ~ vs_phe)) %>% 
  mutate(vs_phe=factor(vs_phe))
z.df <- z.df %>% mutate(vt_phe = fct_cross(vacc_type, vs_phe, sep="_")) %>% 
  mutate(vt_phe = fct_explicit_na(vt_phe, na_level = "Missing")) %>% 
  mutate(vt_phe = fct_recode(vt_phe, uv="AZ_uv", uv="Mo_uv", uv="PB_uv",uv="Missing")) %>% 
  filter(vt_phe != "UNK_v2_14+") %>% 
  mutate(vt_phe = fct_drop(vt_phe))

z_df <- z.df

saveRDS(z_df, paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Sgene_all_positive.rds") )


z_df <- z_df %>%  filter(flag_incon == 0) # remove inconsistent vaccination records
########################################################
#Modelling

#Response var
z.rv <- "hosp_covid" 
z.rv <- "hosp_covid_emerg" 
z.rv.time <- "Time.To.Hosp" 




fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   s_gene"))
z <- pyears(fmla.final, data=z.df, data.frame = TRUE, subset= lab=="lh" & in_eave==1)
z$data %>% mutate(rate.month=event/pyears/12)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  s_gene "))
z.survfit <- survfit(fmla.plot, data=z.df, subset= lab=="lh" & in_eave==1)
plot(z.survfit, fun="event", col=1:4, xlab="Days from Test to Hospital Admission",
     ylab="Risk")
legend("topleft",col=1:4, lty=1, legend=levels(z.df$s_gene))


fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year) +pspline(days) + sex + simd2020_sc_quintile +
                               s_gene +  n_risk_gps +vs "))
z.fit <- coxph(fmla.final , data=z.df, subset = In_Hosp_At_Test==0 & lab=="lh" & in_eave==1)
drop1(z.fit, test="Chisq")

fun.extract <- function(z.fit) {
  #takes a coxph fit using penalised splines and drops off the ps terms
  #make sure no variable begins ps
  z <- summary(z.fit)
  z <- data.frame(z$conf.int)
  z <- z %>% mutate(names = row.names(z)) %>% 
   filter(!(grepl("^ps", names))) %>% 
   dplyr::relocate(names, .before=1) %>% 
   dplyr::select(-exp..coef.)
  names(z) <- c("names","HR","LCL","UCL")
 z
}

#tables
z.rv <- "hosp_covid_emerg" 
z.rv.time <- "Time.To.Hosp" 

#community
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1 & age_year >= 18 & s_gene != "Unknown")

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  s_gene +vt "))
z <- pyears(fmla.final, data=z.df.fit, data.frame=TRUE)
z$data

# nhs labs
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  vs "))
z.survfit <- survfit(fmla.plot, data=z.df, subset= In_Hosp_At_Test==0 & lab=="nhs" & in_eave==1)
plot(z.survfit, fun="event", col=1:5, xlab="Days from Test to Hospital Admission",
     ylab="Risk")
legend("bottomright",col=1:5, lty=1, legend=unique(z.df$vs))


fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year) +pspline(days) + sex + simd2020_sc_quintile +
                                 n_risk_gps +vs "))
z.fit <- coxph(fmla.final , data=z.df, subset = In_Hosp_At_Test==0 & lab=="nhs" & in_eave==1)
summary(z.fit)
drop1(z.fit, test="Chisq")



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

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year) + pspline(days) + sex + simd2020_sc_quintile +
             n_risk_gps + s_gene +  s_gene:vt  "))
z.df.fit <- z.df %>% filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1) %>% 
  filter(s_gene != "Unknown") %>% 
  mutate(s_gene = factor(s_gene)) %>% 
  filter(!(vt %in% c("Mo_v1_0:27" , "Mo_v1_28+v2"))) %>% mutate(vt = fct_drop(vt))
z.fit <- coxph(fmla.final , data=z.df.fit )
summary(z.fit)
z.df %>% group_by(vt) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()
z.df %>% group_by(s_gene) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()
z.df.fit %>% group_by(vt,s_gene) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% as.data.frame()

########################################################################################################

z_df %>% group_by(lab, age_gp) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% 
  mutate(P=R/N*100) %>% ungroup() %>% as.data.frame()
z_df %>% group_by(lab, n_risk_gps) %>% dplyr::summarise(N=n(), R=sum(hosp_covid)) %>% 
  mutate(P=R/N*100) %>% ungroup() %>% as.data.frame()

ptemp <- termplot(z.fit, se=TRUE, plot=FALSE)
ageterm <- ptemp$age_year  # this will be a data frame
center <- with(ageterm, y[x==31])
ytemp <- ageterm$y + outer(ageterm$se, c(0, -1.96, 1.96),'*')
matplot(ageterm$x, exp(ytemp - center), log='y',type='l', lty=c(1,2,2), col=1,
        xlab="Age at Positive Test", ylab="Death Hazard Ratio")
abline(h=1, lty=2)

ageterm <- ptemp$days  # this will be a data frame
center <- with(ageterm, y[x==1])
ytemp <- ageterm$y + outer(ageterm$se, c(0, -1.96, 1.96),'*')
matplot(ageterm$x, exp(ytemp - center), log='y',type='l', lty=c(1,2,2), col=1,
        xlab=paste0("Days from ", format(a_begin, "%d %B" )), ylab="Hazard Ratio of Discharge")
abline(h=1, lty=2)

##############################################################################
z.rv <- "covid_death" 
z.rv.time <- "Time.To.Death" 

#community
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1)
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1 & s_gene %in% c("S_Neg","S_Pos"))
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1  & age_year >= 18)
z.df.fit$s_gene <- fct_drop(z.df.fit$s_gene)
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1 & variant %in% c("alpha","delta"))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + simd2020_sc_quintile +
                               s_gene +  n_risk_gps +vs "))
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + simd2020_sc_quintile +
                               s_gene +  n_risk_gps + s_gene:vacc "))
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  +
                               variant +  n_risk_gps +vs "))
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  +
                               s_gene +  n_risk_gps + s_gene:vs_phe "))
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + simd2020_sc_quintile +
                                n_risk_gps + vt_phe "))

z.fit <- coxph(fmla.final , data=z.df.fit)
summary(z.fit)
drop1(z.fit, test="Chisq")

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  vt_phe "))
z <- pyears(fmla.final, data=z.df.fit, data.frame=TRUE)
z$data

#nhs testing
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="nhs" & in_eave==1 & variant %in% c("alpha","delta") )

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + 
                                +  n_risk_gps + variant + vs "))
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  +
                               variant +  n_risk_gps + variant:vacc "))
z.fit <- coxph(fmla.final , data=z.df.fit)
summary(z.fit)
drop1(z.fit, test="Chisq")
