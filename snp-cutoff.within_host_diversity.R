# If you find this code useful, please cite:
# Coll et al. Definition of a genetic relatedness cutoff to exclude recent transmission of meticillin-resistant Staphylococcus aureus: a genomic epidemiology analysis. Lancet Microbe. 2020 Dec;1(8):e328-e335. doi: 10.1016/S2666-5247(20)30149-X. PMID: 33313577; PMCID: PMC7721685. (https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30149-X/fulltext)

# Intructions on using this R script:
#   1. Change working_dir variable to include the full path of your working directory
#   2. Change full path to SupplementaryData1.xlsx
#   3. The R package gdata needs to be installed, which is used to load .xlsx files
#   4. The R package ggplot2 needs to be installed, which is used for plotting
#   5. The R package lme4 needs to be installed, which is used to run linear mixed models
#   6. The R package svglite needs to be installed, which is used to save plots as SVG files


##########################################################################################################
###                                         1. INPUT FILES                                            ####
##########################################################################################################

working_dir = "";

setwd(working_dir)

dataS1_file = "SupplementaryData1.xlsx"

require(gdata)
require(ggplot2)
require(svglite)

dataS1 = read.xls(dataS1_file, sheet = 1, header = T); # SNP distances derived from mapping to the whole chromosome of the ST22 strain HO 5096 0412 reference genome (and removing MGEs)
dataS1 = read.xls(dataS1_file, sheet = 2, header = T); # SNP distances derived from mapping to the core genome portion of the ST22 strain HO 5096 0412 reference genome (and removing MGEs)
dataS1 = read.xls(dataS1_file, sheet = 3, header = T); # SNP distances derived from mapping to the whole genome portion of the ST30 strain MRSA252 reference genome (and removing MGEs)
dataS1 = read.xls(dataS1_file, sheet = 4, header = T); # SNP distances derived from mapping to the core genome portion of the SST30 strain MRSA252 reference genome (and removing MGEs)

output_sufix = "st22_ref.whole";
output_sufix = "st22_ref.core";
output_sufix = "st30_ref.whole";
output_sufix = "st30_ref.core";

dim(dataS1)
# [1] 1557   12

# Total number of MRSA isolates used
length(unique(c(as.vector(dataS1$SequencingTag1),as.vector(dataS1$SequencingTag2))))
# [1] 1276


##########################################################################################################
###                                     2. PERCENTAGE OF MIXED STRAINS                                ####
##########################################################################################################

# Number of individuals with more than one isolate

length(unique(dataS1$AnonymisedPatientId))
# [1] 459

# Number of individuals with more than one strain (as defined by having isolates from different CCs or
# isolates from the same CC but different clades, labelled as outliers)

length(unique(dataS1$AnonymisedPatientId[which(grepl("outlier",dataS1$Note)==TRUE)]))
# [1] 23

# Removing outliers

dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
dim(dataS1)
# [1] 1510   12

##########################################################################################################
###                           3. SNP DISTANCES AMONG ISOLATES COLLECTED ON THE SAME DAY               ####
##########################################################################################################

# Isolates from the same patient collected on the same day will be used to calculate the cloud of diversity

dataS1sd = dataS1[which(dataS1$TimeGap==0),];

percentiles_95 = vector()

for(r in 1:100)
{

keepInd = vector()

individuals = unique(as.vector(dataS1sd$AnonymisedPatientId))

for(i in 1:length(individuals))
{
  tmp = which(dataS1sd$AnonymisedPatientId==individuals[i])
  # Extracting total number of collection dates available per host
  dates_host = unique(c(as.vector(dataS1sd$CollectionDate1[tmp]), as.vector(dataS1sd$CollectionDate2[tmp])))
  
  # If multiple collection dates are available, select earliest one
  if(length(dates_host)>1)
  {
    earliest_date = min(as.Date(dates_host, format="%Y-%m-%d"))
    tmp = which(dataS1sd$AnonymisedPatientId==individuals[i] & as.Date(dataS1sd$CollectionDate1) == earliest_date)
  }
  # Extracting only one comparison per host
  if(length(tmp)==1)
  {
    keepInd = c(keepInd,tmp)
  } else
  {
    keepInd = c(keepInd,  sample(tmp, 1))
  }
}

# Keeping one isolate pair per patient (the one with the maximum SNP distance)
dataS1sd_max = dataS1sd[keepInd,]
dim(dataS1sd_max)
# [1] 82 12
percentiles_95 = c(percentiles_95, quantile(dataS1sd_max$SNPs, probs = 0.95))

}

quantile(percentiles_95)
# 0%   25%   50%   75%  100% 
# 12.00 12.95 13.90 13.95 13.95  > st22_ref.whole
# 6.90 7.95 8.90 9.85 9.90  > st22_ref.core
# 15.00 20.70 20.70 22.60 22.65 > st30_ref.whole
# 6.9500 7.9500 9.8500 9.8625 9.9000  > st30_ref.core









##########################################################################################################
###                               4. EMPIRICAL DISTRIBUTION OF CLOUD OF DIVERSITY                     ####
##########################################################################################################

# The "cloud of diversity" follows an exponential distribution  

quantile(dataS1sd$SNPs)
# 0%  25%  50%  75% 100% 
# 0    1    3    5   82  > st22_ref.whole

quantile(dataS1sd_max$SNPs)
# 0%  25%  50%  75% 100% 
# 0    1    3    5   82 

quantile(dataS1sd_max$SNPs, probs = 0.95)
# 95% 
# 13.95 

quantile(dataS1sd$SNPs, probs = 0.95)
# 95% 
# 14.85 

## Plots

# Empirical distribution of the cloud of diversity across all CCs

plot_width = 6; plot_height = 5;

plot_cloud_of_diversity = function(data, text_x_offset, plot_title)
{
  size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
  axis_text_size = 15; axis_title_size = 20; ann_text_size = 5;
  
  co_y = round(as.numeric(quantile(data$SNPs, probs = 0.95)));
  
  co_x = which(data$SNPs <= co_y); co_x = co_x[length(co_x)];
  
  g1 <- ggplot(data, aes(x=seq(1,nrow(data),1), y=SNPs)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
    geom_point(shape = 21, colour = dot_color, fill = dot_color, size = size_dot) + 
    ylim(0, 100) + 
    geom_segment(aes(x= 0, y = co_y, xend = co_x, yend = co_y), linetype="dashed", size=size_axis_lines) + 
    geom_segment(aes(x= co_x, y = 0, xend = co_x, yend = co_y), linetype="dashed", size=size_axis_lines) + 
    annotate("text", x = 0 + text_x_offset, y = co_y + text_y_offset, label = paste("95 percentile =",co_y,"SNPs",sep=" "), family=font, size=ann_text_size) + 
    ylab("Number of SNPs") + 
    xlab("Patients") +
    ggtitle(plot_title) +
    theme(text = element_text(family = font)) +
    theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
  
  return(g1)
}

text_x_offset = 15;
dataS1sd = dataS1sd_max[order(as.numeric(dataS1sd_max$SNPs)),];
plot_title = "Empirical Cloud of Diversity"
g1 = plot_cloud_of_diversity(dataS1sd,text_x_offset,plot_title)
plot_file = paste("empirical_clould_of_diversity.allCCs.",output_sufix,".pdf",sep="");
ggsave(plot_file, plot = g1, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")




##########################################################################################################
###                                   6. APPLYING LINEAR MIXED MODELS                                 ####
##########################################################################################################

# Linear mixed models are applied to calculate the SNP accumulation rate and to model the "cloud of diversity"
# The number of SNPs between MRSA isolates (SNPs) is modelled as a function of the time gap (TimeGap) between isolates
# The intercept (that is, number of SNPs at time 0) is interpreted as the "cloud of diversity" at time 0 and assumbed to 
# vary by patient (AnonymisedPatientId, random variable)

library(lme4)
lmer_all  =  lmer(SNPs  ~ TimeGap   + (1|AnonymisedPatientId), data=dataS1)
summary(lmer_all)


##########################################################################################################
###                                    5. PLOT SNP DISTANCES OVER TIME                                ####
##########################################################################################################

# Binning data point by TimeGap in months

bins_from = seq(0,330,30); bins_to = seq(30,360,30); bin_month = seq(1,12,1);
dataS1$bin = NA;

for(b in 1:length(bin_month))
{
  tmp = which(dataS1$TimeGap>=bins_from[b] & dataS1$TimeGap<bins_to[b])
  if(length(tmp)>0){ dataS1$bin[tmp] = as.character(bin_month[b]); }
}

# Converting month bin label to factor
dataS1$bin = factor(dataS1$bin,seq(1,12,1))

# Creating X labels
xlab <- paste(levels(dataS1$bin),"\n(N=",table(dataS1$bin),")",sep="")
size_axis_lines = 0.3; axis_x_text_size = 8; axis_y_text_size = 15; axis_title_size = 20;
plot_title = "Number of SNPs over time"

# Boxplot of binned SNP distances per month
boxplot = ggplot(dataS1,aes(x = bin, y = SNPs)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",
        size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_boxplot(varwidth = TRUE, notch=FALSE, outlier.size = 0.5) + 
  scale_y_continuous(limits=c(0,25)) +
  scale_x_discrete(labels=xlab) + 
  ylab("Number of SNPs") + 
  xlab("Time distance in months\n(Pairwise Comparisons)") +
  ggtitle(plot_title) +
  theme(text = element_text(family = font)) +
  theme(axis.text.x = element_text(size=axis_x_text_size, color="black"), axis.text.y = element_text(size=axis_y_text_size, color="black"), 
        axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))

boxplot_file = paste("snp_distances_over_time.allCCs.",output_sufix,".svg",sep="");
ggsave(boxplot_file, plot = boxplot, device = "svg", width = 6, height = 5, dpi = 300, units = "in")


##########################################################################################################
###                                 6. CALCULATION OF THE MRSA SUBSTITUTION RATE                      ####
##########################################################################################################

## Converting units to substitutions per site per year

chromosome_length = 2832299; # ST22 chromosome reference 
mge_regions_length = 159127; # total length of MGEs on the ST22 chromosome
chromosome_length = 1759534; # ST22 core-genome reference 
mge_regions_length = 2809; # MGEs on the ST22 core-genome reference

### Using MRSA isolates from all CCs
timegap_coefficient_all = coef(summary(lmer_all))[2,1]
timegap_coefficient_all
# [1] 0.01120942
# [1] 0.007115821

## Converting units to SNPs per genome per year
substitution_rate_all = timegap_coefficient_all*365
substitution_rate_all
# [1] 4.091438
# [1] 2.597275

## 95% confidence interval
CI = confint(lmer_all, "TimeGap", level = 0.95)
CI_lower_bound = CI[1]*365
CI_lower_bound
# [1] 2.409668
# [1] 1.475942
CI_upper_bound = CI[2]*365
CI_upper_bound
# [1] 5.774858
# [1] 3.719656

## Converting units to substitutions per site per year
substitution_rate_all_ps = substitution_rate_all/(chromosome_length-mge_regions_length);
substitution_rate_all_ps
# [1] 1.530556e-06
# [1] 1.478475e-06

CI_lower_bound_ps = CI_lower_bound/(chromosome_length-mge_regions_length)
CI_lower_bound_ps
# [1] 9.014263e-07
# [1] 8.401667e-07

CI_upper_bound_ps = CI_upper_bound/(chromosome_length-mge_regions_length)
CI_upper_bound_ps
# [1] 2.160302e-06
# [1] 2.117381e-06


##########################################################################################################
###                      APPLYING LINEAR MIXED MODELS WITH SUB-SAMPLED DATA                          ####
##########################################################################################################

# Because patients differ in the number of isolate genomes available per sample (colonies sequenced)
# the dataset needs to be deduplicated to keep only isolate per sample This is done by:
#   1. Keeping all different samples (i.e. collection dates) per patient
#   2. Keeping 1 isolate genome per sample
#   3. If a patient has a single sample (i.e. single collection date), two isolates from this sample are kept
# Isolates are selected randomly

library(lme4)
number_iterations = 100
mutation_rates_iterations = vector()
mutation_rates_lci_iterations = vector()
mutation_rates_uci_iterations = vector()
beta_zero_iterations = vector()

host_ids = unique(as.vector(dataS1$AnonymisedPatientId))
length(host_ids)
# [1] 445

for(i in 1:number_iterations)
{
  print(paste("Iteration number: ",i,sep=""))
  # Vector to store pairwise isolate comparisons to keep in each iteration
  host_isolates_kept = vector()
  
  # Across all hosts, sub-sample randomly to de-duplicate dataset
  for(h in 1:length(host_ids))
  {
    hhh = which(dataS1$AnonymisedPatientId == host_ids[h])
    # Sample only one CC
    host_clonal_complexes = as.character(unique(c(as.vector(dataS1$ClonalComplex1[hhh]), as.vector(dataS1$ClonalComplex2[hhh]))))
    host_clonal_complex = sample(host_clonal_complexes,1)
    # Select comparison of sampled CC
    hhh = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$ClonalComplex1 == host_clonal_complex)
    # Get all available collection dates
    host_collection_dates = unique(c(as.vector(dataS1$CollectionDate1[hhh]),as.vector(dataS1$CollectionDate2[hhh])))
    # If only one available collection data/sample > keep two random isolates
    if(length(host_collection_dates)==1)
    {
      host_isolates = vector()
      hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate1 == host_collection_dates[1])
      if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(dataS1$SequencingTag1[hhhd1])); }
      hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate2 == host_collection_dates[1])
      if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(dataS1$SequencingTag2[hhhd1])); }
      host_isolates = unique(host_isolates)
      host_isolates = sample(host_isolates, 2);
      host_isolates_kept = c(host_isolates_kept, host_isolates)
    } else
    {
      # Else, for each collection date, randomly sample one isolate
      for(d in 1:length(host_collection_dates))
      {
        host_isolates = vector()
        hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate1 == host_collection_dates[d])
        if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(dataS1$SequencingTag1[hhhd1])); }
        hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate2 == host_collection_dates[d])
        if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(dataS1$SequencingTag2[hhhd1])); }
        host_isolates = unique(host_isolates)
        host_isolates = sample(host_isolates, 1);
        host_isolates_kept = c(host_isolates_kept, host_isolates)
      }
    }
  }

  print(paste("Number of isolates sub-sampled: ",length(host_isolates_kept), sep=""))
  
  # Keeping comparisons including isolates sub-sampled
  iii1 = which(!is.na(match(dataS1$SequencingTag1, host_isolates_kept)))
  iii2 = which(!is.na(match(dataS1$SequencingTag2, host_isolates_kept)))
  iii = iii1[which(!is.na(match(iii1,iii2)))]
  print(paste("Number of pairwise comparisons sub-sampled: ",length(iii), sep=""))
  print(paste("Number of patients sub-sampled: ",length(unique(dataS1$AnonymisedPatientId[iii])), sep=""))
  tmp = match(host_ids, dataS1$AnonymisedPatientId[iii])
  print(paste("Missing sub-sampled patients: ",paste(host_ids[which(is.na(tmp))], collapse = ";"), sep=""))
  
  ### Running linear mixed model
  dataS1_sub = dataS1[iii,]
  dataS1_sub <- droplevels(dataS1_sub)
  dataS1_sub$AnonymisedPatientId=as.factor(dataS1_sub$AnonymisedPatientId)
  # dataS1_sub = rbind(dataS1_sub, dataS1_sub)
  lmer_all  =  lmer(SNPs  ~ TimeGap   + (1|AnonymisedPatientId), data=dataS1_sub)
  timegap_coefficient_all = coef(summary(lmer_all))[2,1]
  mutation_rates_iterations = c(mutation_rates_iterations, timegap_coefficient_all)
  CI = confint(lmer_all, "TimeGap", level = 0.95)
  CI_lower_bound = CI[1]; mutation_rates_lci_iterations = c(mutation_rates_lci_iterations, CI_lower_bound);
  CI_upper_bound = CI[2]; mutation_rates_uci_iterations = c(mutation_rates_uci_iterations, CI_upper_bound);
  beta0_all = as.vector(unlist(coef(lmer_all)$AnonymisedPatientId["(Intercept)"]))
  beta0_95per = quantile(beta0_all,prob=0.95)
  beta_zero_iterations = c(beta_zero_iterations, beta0_95per)
}

##### Mutation rates and 95% cloud of diversity across all 100 iterations
quantile(mutation_rates_iterations*365)
# 0%         25%         50%         75%        100% 
# 4.431645 4.607433 4.707698 4.817113 4.913721  > st22_ref.whole
# 2.727992 2.852357 2.935401 2.988398 3.095441  > st22_ref.core
# 4.513096 4.749575 4.902059 5.119699 5.441598  > st30_ref.whole
# 2.704350 2.839135 2.887466 2.952478 3.069959  > st30_ref.core

quantile(beta_zero_iterations)
# 0%      25%      50%      75%     100% 
# 18.51329 18.93635 19.22067 19.47765 19.72289  > st22_ref.whole
#  9.776901 10.004914 10.452521 10.822996 10.846875 > st22_ref.core
# 18.18067 19.29616 21.44925 22.16469 22.29634  > st30_ref.whole
# 10.02485 10.08829 10.36068 10.61637 10.66863  > st30_ref.core



### Extracting median substitution rate and 95% CI across iterations
if(output_sufix == "st22_ref.whole")
{
  chromosome_length = 2832299; # ST22 chromosome reference 
  mge_regions_length = 159127; # total length of MGEs on the ST22 chromosome
}
if(output_sufix == "st22_ref.core")
{
  chromosome_length = 1759534; # ST22 core-genome reference 
  mge_regions_length = 2809; # MGEs on the ST22 core-genome reference
}
if(output_sufix == "st30_ref.whole")
{
  chromosome_length = 2902619; # ST30 chromosome reference 
  mge_regions_length = 314479; # total length of MGEs on the ST30 chromosome
}
if(output_sufix == "st30_ref.core")
{
  chromosome_length = 1754228; # ST30 core-genome reference 
  mge_regions_length = 3866; # MGEs on the ST30 core-genome reference
}

mutation_rate_median = sort(mutation_rates_iterations)[50]
tmp = which(mutation_rates_iterations == mutation_rate_median)
mutation_rate_lci = mutation_rates_lci_iterations[tmp[1]]
mutation_rate_uci = mutation_rates_uci_iterations[tmp[1]]
print(paste(mutation_rate_median," (",mutation_rate_lci," - ",mutation_rate_uci,")",sep=""))
# [1] "0.0128724506584772 (0.00781873051982872 - 0.0179269008529202)"  > st22_ref.whole
# [1] "0.0134591253302896 (0.0070021257138012 - 0.0199192055155959)"  > st30_ref.whole
# [1] "0.00800455270067963 (0.00463717411557472 - 0.0113724084147899)" > st22_ref.core
# [1] "0.00791734559036188 (0.00448982538564908 - 0.0113475954839798)" > st30_ref.core

# Converting units to substitutions per site per year
print(paste(mutation_rate_median*365," (",mutation_rate_lci*365," - ",mutation_rate_uci*365,")",sep=""))
# [1] "4.69844449034419 (2.85383663973748 - 6.54331881131588)"  > st22_ref.whole
# [1] "4.9125807455557 (2.55577588553744 - 7.2705100131925)" > st30_ref.whole
# [1] "2.92166173574806 (1.69256855218477 - 4.15092907139832)" > st22_ref.core
# [1] "2.88983114048209 (1.63878626576192 - 4.14187235165264)"  > st30_ref.core
print(paste(mutation_rate_median*365/(chromosome_length-mge_regions_length)," (",mutation_rate_lci*365/(chromosome_length-mge_regions_length)," - ",mutation_rate_uci*365/(chromosome_length-mge_regions_length),")",sep=""))
# [1] "1.75762894806028e-06 (1.06758436783622e-06 - 2.44777321149401e-06)"  > st22_ref.whole
# [1] "1.89811244583203e-06 (9.87495222645389e-07 - 2.80916411523043e-06)" > st30_ref.whole
# [1] "1.66312982154183e-06 (9.63479515681038e-07 - 2.36287926192109e-06)" > st22_ref.core
# [1] "1.65099056108513e-06 (9.36255623557821e-07 - 2.36629471598026e-06)"  > st30_ref.core





##########################################################################################################
###                                 8. MODELED DISTRIBUTION OF CLOUD OF DIVERSITY                     ####
##########################################################################################################

# Extracting varying intercepts for all patients

beta0_all = as.vector(unlist(coef(lmer_all)$AnonymisedPatientId["(Intercept)"]));
quantile(beta0_all)
# 0%        25%        50%        75%       100% 
# -0.5326214  2.2241942  4.2607700  7.8025847 56.8797944 
# -0.6022719  1.2223954  2.0550299  3.8494356 34.6246406 
quantile(beta0_all,prob=0.95)
# 95% 
# 19.43245 
# 10.51029 


## Plots

# Modelled distribution of the cloud of diversity across all CCs
beta0_all_df = as.data.frame(cbind(seq(1,length(beta0_all),1),sort(beta0_all)))
colnames(beta0_all_df) = c("X","SNPs")
text_x_offset = 70;
plot_title = "Modelled Cloud of Diversity";
g3 = plot_cloud_of_diversity(beta0_all_df,text_x_offset,plot_title)
plot_file = paste("modelled_clould_of_diversity.allCCs.",output_sufix,".pdf",sep="");
ggsave(plot_file, plot = g3, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")




##########################################################################################################
####****                          9. SIMULATION MODEL OF TRANSMISSION                             ****####
##########################################################################################################
# Load in the required functions
# These are the distribution sampler / the curve for the number of SNPs over time / the simulation model
source("simu_transmission_fn.R")
require(reshape2)

## Simulation population 
npat = 459 # number of patient samples in cohort 1
nruns = 200000 # number of transmission samples - can be increased to reduce variation in number of SNPs

ndays = 180 # time between samples


## Choose which mapping: All of CC22/CC30 or core genome only
map <- "CC22" # ST22 strain HO 5096 0412 mapping data 
#map <- "CC30"
#map <- "CC22_core"
#map <- "CC30_core"

## mu = substitution_rate
if(map == "CC22"){mu = 4.7 / 365}
if(map == "CC30"){mu = 4.9 / 365}
if(map == "CC22_core"){mu = 2.9 / 365}
if(map == "CC30_core"){mu = 2.9 / 365}

## Parameters from model fit
param_general_fit <- read.csv(paste0("output/param_general_fit_",map,".csv"))[,2]

## Data for time zero transferred variability
## Which sheet? 
if(map == "CC22"){sheet_num = 1}
if(map == "CC30"){sheet_num = 2}
if(map == "CC22_core"){sheet_num = 3}
if(map == "CC30_core"){sheet_num = 4}
dataS1 = read.xls(dataS1_file, sheet = sheet_num, header = T) 
# Removing outliers
dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
# Same day 
dd0 = dataS1[which(dataS1$TimeGap == 0),] # 104 pairs

h <- hist(dd0$SNPs, breaks = seq(0,max(dd0$SNPs),1))
t0prob_dist <- h$counts/sum(h$counts)
w <- which(t0prob_dist <= 0.01) ## Remove those at < 1%: has big impact on rMyDist.
t0prob_dist[w] <- 0

nonz <- which(t0prob_dist > 0) ## Make the length of the vector equal to the max number of SNPs (remove those 0 values above this)
t0prob_dist <- t0prob_dist[1:max(nonz)]

###### Sample from the baseline variance in the source patient at transmission ##############################
## Run the simulation 10 times to give a range on the maximum, some random variation expected, depending on the number of runs used. 
m <- rep(0,10)
for(i in 1:10){
  ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit, t0prob_dist)
  
  ## Maximum number of SNPs needed to capture 95% or 99% of the transmission events
  m[i] <- max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]) # = below this, capture 95% of all transmissions
  
}


## Plot the output from the last run
# g5 = ggplot(ss$store_limits, aes(x=value, fill = variable)) + geom_histogram(aes(y=..density..), binwidth = 1, position = "identity") +   
#   facet_wrap(~variable) + guides(fill=FALSE) + scale_y_continuous(paste0("Density across ", nruns, " simulations")) + scale_x_continuous("Number of SNPs")
# 
# plot_file = paste0("simulation_model_distribution_of_SNPs_",map,".pdf");
# 
# ggsave(plot_file, plot = g5, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")


## Results: 
max(m)
range(m)

##            95% of transmission events (max [range])
## map = CC22      = 17 (16 - 17) 
## map = CC30      = 55 (52 - 55)
## map = CC22_core = 12 (11 - 12)
## map = CC30_core = 12 (12 - 12) 

