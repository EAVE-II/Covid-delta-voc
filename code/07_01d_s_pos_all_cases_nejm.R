##########################################################
# Name of file: 07_01d_s_pos_all_cases_nejm.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@phs.net
# Original date: 31 August 2021
# Latest update author (if not using version control) - Chris Robertson chrisobertson@phs.net
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: sets up analysis for all positive tests
#                         run 07_01a_s_gene_positive_description.R
#                         to read in the data - line 291
#                         then get first positive test post 01 April
#                         run 07_01d_s_gene_positive_all_cases.R to get the data
# Approximate run time: Unknown
##########################################################

z.df <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Sgene_all_positive.rds"))
a_begin <- as.Date("2021-04-01")
z.df <- z.df %>%  filter(age_year >= 18) # remove inconsistent vaccination records - need to correct this
z.df <- z.df %>%  filter(flag_incon == 0) # remove inconsistent vaccination records - need to correct this

z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1)

table(z.df.fit$s_gene,z.df.fit$covid_death)

z <- z.df.fit %>% group_by(NRS.Date.Death) %>% dplyr::summarise(N=sum(covid_death))
g1 <- z %>% ggplot(aes(x=NRS.Date.Death, y=N)) + geom_point() +
  labs(y="Number of Deaths", x="Date of Death", title="Daily Covid Deaths in Scotland", subtitle="Among adults tested in the community")

z <- z.df.fit %>% group_by(specimen_date) %>% 
  dplyr::summarise(N=n(), N_S_Pos = sum(s_gene=="S_Pos"), N_Seq = sum(variant != "not_sequenced"), N_delta=sum(variant=="delta")) %>% 
  mutate(N_Seq = if_else(specimen_date > as.Date("2021-08-06"), NA_integer_, N_Seq),
         N_delta = if_else(specimen_date > as.Date("2021-08-06"), NA_integer_, N_delta)) %>% 
  mutate(P_S_Pos = N_S_Pos/N*100, P_delta = N_delta/N_Seq*100) %>% 
  pivot_longer(cols=N:P_delta)

g2 <- z %>% filter(name %in% c("N","N_Seq")) %>% 
  ggplot(aes(x=specimen_date, y=value, colour=name)) + geom_point() + 
  labs(y="Number", x= "Date of Specimen", title="Cases and Sequenced in Scotland") +
  scale_colour_manual(name="", values= c("black","blue"), breaks=c("N","N_Seq"), labels=c("Positive","Sequenced"))

g3 <-z %>% filter(name %in% c("P_S_Pos","P_delta")) %>% 
  ggplot(aes(x=specimen_date, y=value, colour=name)) + geom_point() +
  labs(y="Percentage", x= "Date of Specimen", title="Percentage S positive and delta in Scotland",
       subtitle="Denominators - all cases and all sequenced, respectively") +
  scale_colour_manual(name="", values= c("black","blue"), breaks=c("P_S_Pos","P_delta"), 
                      labels=c("S Positive","Sequenced delta"))

gridExtra::grid.arrange(g2,g3, g1, nrow=2, ncol=2)


#models
#Table S2
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1 & s_gene %in% c("S_Neg","S_Pos"))
z.df.fit$s_gene <- fct_drop(z.df.fit$s_gene)

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + simd2020_sc_quintile +
                               s_gene +  n_risk_gps +vs "))
z.fit <- coxph(fmla.final , data=z.df.fit)
summary(z.fit)

#Table S3
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1 & s_gene %in% c("S_Neg","S_Pos"))
z.df.fit$s_gene <- fct_drop(z.df.fit$s_gene)

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + simd2020_sc_quintile +
                               s_gene +  n_risk_gps + s_gene:vs_phe "))

z.fit <- coxph(fmla.final , data=z.df.fit)
summary(z.fit)

#Table S4
z.df.fit <- z.df %>%  filter(In_Hosp_At_Test==0 & lab=="lh" & in_eave==1 & s_gene %in% c("S_Pos"))
z.df.fit$s_gene <- fct_drop(z.df.fit$s_gene)

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + simd2020_sc_quintile +
                                n_risk_gps + vs_phe "))

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(age_year)  + sex  + simd2020_sc_quintile +
                                n_risk_gps + vt_phe "))

z.fit <- coxph(fmla.final , data=z.df.fit)
summary(z.fit)

