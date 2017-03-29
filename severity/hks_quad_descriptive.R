source("hks_quad_setup.r")

#source("hks_quad_setup.R")

# Characteristics of the cohort
# Which groups am I going to use and how should I summerize them?

# Some descriptive stats for the quality of the cohort for the variables we are interested in
# Need to make sure that the date variables are being loaded correctly
max(allDat$blood_1_date,na.rm=TRUE)
min(allDat$blood_2_date,na.rm=TRUE)
max(allDat$blood_2_date,na.rm=TRUE)
min(allDat$blood_3_date,na.rm=TRUE)
max(allDat$blood_3_date,na.rm=TRUE)
min(allDat$blood_4_date,na.rm=TRUE)
max(allDat$blood_4_date,na.rm=TRUE)

# What were the characteristics of the study
names(allDat)

# Total number of households
#length(table(substr(allDat$New_id,4,7)))

# Total number of people in households in which at least one person participated
#overallN <- dim(allDat)[1]

#sum(table(allDat$ag1)) / overallN

# Currently being used for Table 1

  # Create work status
 # allDat$occstat.r1[allDat$occ.r1 %in% c(1,2,3,4,5,6,7,8,9) ] <- "atwork"
 # allDat$occstat.r1[allDat$occ.r1 %in% c(10,11,13) ] <- "athome"
 # allDat$occstat.r1[allDat$occ.r1 %in% c(12) ] <- "student"
 # allDatVac$occstat.r1[allDatVac$occ.r1 %in% c(1,2,3,4,5,6,7,8,9) ] <- "atwork"
 # allDatVac$occstat.r1[allDatVac$occ.r1 %in% c(10,11,13) ] <- "athome"
 # allDatVac$occstat.r1[allDatVac$occ.r1 %in% c(12) ] <- "student"
 # allDatNoVac$occstat.r1[allDatNoVac$occ.r1 %in% c(1,2,3,4,5,6,7,8,9) ] <- "atwork"
 # allDatNoVac$occstat.r1[allDatNoVac$occ.r1 %in% c(10,11,13) ] <- "athome"
 # allDatNoVac$occstat.r1[allDatNoVac$occ.r1 %in% c(12) ] <- "student"

#cbind(table(allDat$ag1) ,  table(allDatVac$ag1) , table(allDatNoVac$ag1) )
#cbind(table(allDat$sex.r1) ,  table(allDatVac$sex.r1) , table(allDatNoVac$sex.r1) )
#cbind(table(allDat$occstat.r1[allDat$Age_1>=15]) ,  table(allDatVac$occstat.r1[allDatVac$Age_1>=15]) , table(allDatNoVac$occstat.r1[allDatNoVac$Age_1>=15]) )
#cbind(table(allDat$occstat.r1[allDat$Age_1<15]) ,  table(allDatVac$occstat.r1[allDatVac$Age_1<15]) , table(allDatNoVac$occstat.r1[allDatNoVac$Age_1<15]) )

#cbind(table(allDat$edu_a.r1[allDat$occ.r1 !=12]) ,  table(allDatVac$edu_a.r1[allDatVac$occ.r1 !=12]) , table(allDatNoVac$edu_a.r1[allDatNoVac$occ.r1 !=12]) )
  # 1 = HKI ; 2,3 = KLN ; 4,5,6 = NT ; 7 = missing
#table(allDat$district.r1)
#table(allDatVac$district.r1)
#table(allDatNoVac$district.r1)

  # Compare vaccinated and non-vaccinated

#  prop.test(x=c(11,31), n=c(113 , 306))
#  prop.test(x=c(19,82), n=c(113 , 306))
#  prop.test(x=c(44,177), n=c(113 , 306))
#  prop.test(x=c(39,16), n=c(113 , 306))

#  prop.test(x=c(62,186), n=c(113 , 306))
#  prop.test(x=c(51,120), n=c(113 , 306))

#  prop.test(x=c(55,145), n=c(113 , 306))
#  prop.test(x=c(40,118), n=c(113 , 306))
#  prop.test(x=c(12,26), n=c(113 , 306))
#  prop.test(x=c(6,17), n=c(113 , 306))

#  fisher.test(  matrix( c(1,112,4,302      )   ,ncol=2) )
#  prop.test(x=c(27,38), n=c(113 , 306))
#  prop.test(x=c(39,158), n=c(113 , 306))
#  prop.test(x=c(6,21), n=c(113 , 306))
#  prop.test(x=c(22,40), n=c(113 , 306))
#  prop.test(x=c(18,43), n=c(113 , 306))

#  prop.test(x=c(16,61), n=c(113 , 306))
#  prop.test(x=c(39,95), n=c(113 , 306))
#  prop.test(x=c(56,143), n=c(113 , 306))
#  fisher.test(  matrix( c(0,113,0,306      )   ,ncol=2) )


# Make plots of seroconversion rates and generate text for the latex table
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r12_ph1n1,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_ph1n1_r12.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r13_ph1n1,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_ph1n1_r13.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r14_ph1n1,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_ph1n1_r14.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r23_ph1n1,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_ph1n1_r23.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r24_ph1n1,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_ph1n1_r24.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r34_ph1n1,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_ph1n1_r34.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r12_sh3n2,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_sh3n2_r12.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r13_sh3n2,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_sh3n2_r13.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r14_sh3n2,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_sh3n2_r14.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r23_sh3n2,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_sh3n2_r23.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r24_sh3n2,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_sh3n2_r24.pdf")
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r34_sh3n2,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_sh3n2_r34.pdf")

# Some confidence intervals among the positives, i.e. for the distribtion
d_r12_h1 <- c(10,7,5,0)
d_r24_h1 <- c(11,18,35,2)
d_r23_h3 <- c(2,4,16,2)

chisq.test(as.matrix(rbind(c(10,7,5),c(11,18,37))))
chisq.test(as.matrix(rbind(c(11,18,37),c(2,4,18))))

sum(d_r12_h1)
sum(d_r24_h1)
tmpcis <- function(x){
	cat(format(100*srg.ci.binom(rep(sum(x),length(x)),x,0.975),digits=3))
	cat("\n")
	cat(format(100*x/rep(sum(x),length(x)),digits=3))
	cat("\n")
	cat(format(100*srg.ci.binom(rep(sum(x),length(x)),x,0.275),digits=3))
}
tmpcis(d_r12_h1)
tmpcis(d_r24_h1)

# Lines below are for the supporitng figure, not the table


hks.round.CARs(vecAg1Study,table(qDatNoVac$PH1N1_com_result_tbase,qDatNoVac$ag1)[2,],vecAg1HK, filename="bar_micro_neut_ph1n1_r12.pdf",ymax=0.6)
hks.round.CARs(vecAg1Study,table(qDatNoVac$ff_r12_ph1n1,qDatNoVac$ag1)["TRUE",],vecAg1HK, filename="bar_ph1n1_r12_comp_micro.pdf",ymax=0.6)

# Generate an arrow plot of the unvaccinated individuals who showed a 4-fold rise or
# greater to two strains in the same time step
# Think about doing other versions of this for different vaccination sub groups
#hks.coinf.arrows(qDatNoVac,"hks_coinf_arrows.pdf")


# Make some age specific and overall infection rate charts
# The auxcases numbers come from the docx in the data dir
hks.wave.strain.severity(
		"sev_chart_all_cause_death.pdf",
		"sev_table_all_cause_death.tex",
		qDatNoVac,
		vecAg1HK,
		auxcases = c(61,160,359,73,409,37))

hks.wave.strain.severity(
		"sev_chart_resp_death.pdf",
		"sev_table_resp_death.tex",
		qDatNoVac,
		vecAg1HK,
		auxcases=c(57,130,289,34,190,17))

hks.wave.strain.severity(
		"sev_chart_resp_hosp.pdf",
		"sev_table_resp_hosp.tex",
		qDatNoVac,
		vecAg1HK,
		auxcases=c(3472,2480,4589,942,5299,476))

# Plot long timeseries of chp data and also round / wave plot for Fig 1
#hks.plot.long.gp("hks_plot_long_chp.pdf",aDatGP)
#hks.plot.waves.aux(aDatCHP,qDatNoVac,"data_and_aux_raw.pdf")

# Remember the CHP data and the possible survival function approach
#hks.fig.qmh.wave.survival("fig_QMH_Wave_Survival.pdf","QMH_surveillance.csv")

# Visualize the raw data with lines and lasagne
# Matierals for Fig 2
# Some of these are broken
hks.fig.cards("house_of_cards_",qDatNoVac,strain="H1N1",change="Change")
hks.fig.cards("house_of_cards_",qDatNoVac,strain="H1N1",change="NoChange")
hks.fig.cards("house_of_cards_",qDatNoVac,strain="H3N2",change="Change")
hks.fig.cards("house_of_cards_",qDatNoVac,strain="H3N2",change="NoChange")
hks.lasagne.plot(qDatNoVac,"H1N1","lasagne_")
hks.lasagne.plot(qDatNoVac,"H3N2","lasagne_")

# Inter subtype specific innumity
chisq.test(table(qDatNoVac$ff_r12_ph1n1,qDatNoVac$ff_r24_ph1n1))
