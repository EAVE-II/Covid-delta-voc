z <- filter(z.df, hosp_covid==1)
by(z$Time.To.Hosp, list(z$ageYear <= 65), summary)

table(z$vs, z$vacc_type, exclude=NULL)
df_cohort %>% group_by(vacc_type) %>% dplyr::summarise(R1=sum(eave_weight), 
                                                       Dose_2=sum(eave_weight*!is.na(date_vacc_2))) %>% 
  mutate(Dose_1 = R1-Dose_2) %>% 
  dplyr::select(vacc_type, Dose_1, Dose_2)

table(z$s_gene, z$vs)
table(z$age_gp, z$vs)

hist(z$ageYear)

z %>% ggplot(aes(x=ageYear, fill=vs)) + geom_histogram() + facet_wrap(~s_gene) +
  labs(y="Count", x="Age", fill = "Vaccination")

z %>% ggplot(aes(x=n_risk_gps, fill=vs)) + geom_bar() + facet_wrap(~s_gene) +
  labs(y="Count", x="Number of Risk Groups", fill = "Vaccination")


z <- z.df %>% filter(hosp_covid==1) %>% 
  group_by(date_hosp_covid, s_gene) %>% 
  dplyr::summarise(N=n())

z %>%  ggplot(aes(x=date_hosp_covid, y=N, fill=s_gene)) + geom_col() +
  labs(x="Date Admitted to Hospital", y="Number Covid Admissions", fill= "s Gene Status",
       title="Hospital Admissions (S Gene status linked to EAVE-II)")

z <- z.df %>% filter(hosp_covid==1) %>% 
  group_by(date_hosp_covid, s_gene, vacc) %>% 
  dplyr::summarise(N=n()) %>% 
  mutate(s_gene = case_when(s_gene=="S Gene Negative" ~ "S-",
                            s_gene=="S Gene Positive" ~ "S+",
                            TRUE ~ "WP")) %>% 
  mutate(var=paste(s_gene, vacc))

z %>%  ggplot(aes(x=date_hosp_covid, y=N, fill=var)) + geom_col() +
  labs(x="Date Admitted to Hospital", y="Number Covid Admissions", 
       fill= "S Gene + Vaccine",
       title="Hospital Admissions (known S Gene status linked to EAVE-II cohort)")


########################################
#Vaccination by April 01

z <- df_cohort %>% dplyr::select(EAVE_LINKNO, ageYear, age_gp, simd2020_sc_quintile, date_vacc_1, date_vacc_2, vacc_type, eave_weight) %>% 
  mutate(vacc_1 = if_else(!is.na(date_vacc_1) & date_vacc_1 <= as.Date("2021-04-01"), 1, 0),
         vacc_2 = if_else(!is.na(date_vacc_2) & date_vacc_2 <= as.Date("2021-04-01"), 1, 0) )

z_tab <- z %>% #group_by(ageYear<= 64) %>% 
  dplyr::summarise(v1 = sum(vacc_1), v2=sum(vacc_2), N=sum(eave_weight)) %>% 
  mutate(p1=v1/N, p2=v2/N)

z_tab <- z %>% filter(vacc_1==1 & ageYear >= 16) %>% group_by(age_gp, vacc_type) %>% dplyr::summarise(N=n())
z_tab %>% ggplot(aes(x=age_gp, y=N, fill=vacc_type)) + geom_col() + labs(x="Age Group", y="Count", fill="Vaccine")


z_tab <- z %>% filter(vacc_1==1 & ageYear >= 16 & simd2020_sc_quintile != "NA") %>%
  group_by(simd2020_sc_quintile, vacc_type) %>% dplyr::summarise(N=n())
z_tab %>% ggplot(aes(x=simd2020_sc_quintile, y=N, fill=vacc_type)) + geom_col(position="dodge") +
  labs(x="Deprivation Quintile", y="Count", fill="Vaccine")
