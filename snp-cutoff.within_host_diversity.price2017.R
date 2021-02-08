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

dataS1_file = "SupplementaryData3.xlsx"

require(gdata)
require(ggplot2)
require(svglite)

output_sufix = "st22_ref.whole";
output_sufix = "st22_ref.core";
output_sufix = "st30_ref.whole";
output_sufix = "st30_ref.core";

if(output_sufix == "st22_ref.whole"){ dataS1 = read.xls(dataS1_file, sheet = 1, header = T); }
if(output_sufix == "st22_ref.core"){ dataS1 = read.xls(dataS1_file, sheet = 2, header = T); }
if(output_sufix == "st30_ref.whole"){ dataS1 = read.xls(dataS1_file, sheet = 3, header = T); }
if(output_sufix == "st30_ref.core"){ dataS1 = read.xls(dataS1_file, sheet = 4, header = T); }

dim(dataS1)
# [1] 12058    11

# Total number of MRSA isolates used
length(unique(c(as.vector(dataS1$SequencingTag1),as.vector(dataS1$SequencingTag2))))
# [1] 1537 > price2017



##########################################################################################################
###                                     2. PERCENTAGE OF MIXED STRAINS                                ####
##########################################################################################################

# Number of individuals with more than one isolate
length(unique(dataS1$AnonymisedPatientId))
# [1] 280 > price2017

# Removing outliers
dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
dim(dataS1)
# [1] 11167    13 > price2017, st22_ref.whole


##########################################################################################################
###                           3. SNP DISTANCES AMONG ISOLATES COLLECTED ON THE SAME DAY               ####
##########################################################################################################

# Isolates from the same patient collected on the same day will be used to calculate the cloud of diversity

dataS1sd = dataS1[which(dataS1$TimeGap==0),];
dim(dataS1sd)
# [1] 3455   10

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
  # [1] 146  10
  percentiles_95 = c(percentiles_95, quantile(dataS1sd_max$SNPs, probs = 0.95))
}
quantile(percentiles_95)
# 0%   25%   50%   75%  100% 
# 18.75 19.75 21.50 22.25 22.75 > st22_ref.whole > all
# 18.00 20.50 21.25 22.00 24.00 > st30_ref.whole > all
#  9.75 10.00 10.75 10.75 12.00 > st22_ref.core > all
#  8.75  9.75 11.50 11.50 12.00  > st30_ref.core > all

# Keeping one isolate pair per patient (the one with the maximum SNP distance)
dataS1sd_max = dataS1sd[keepInd,]
dim(dataS1sd_max)
# [1] 146  10 > price2017 

##########################################################################################################
###                               4. EMPIRICAL DISTRIBUTION OF CLOUD OF DIVERSITY                     ####
##########################################################################################################

# The "cloud of diversity" follows an exponential distributio
quantile(dataS1sd_max$SNPs)
# 0%  25%  50%  75% 100% 
# 0    3    6   11   43  > price2017, st22_ref.whole 
#  0.00  3.00  6.00 10.75 44.00 > price2017, st30_ref.whole 
#  0    0    1    3   20 > price2017, st22_ref.core 
#  0.0  0.0  0.5  3.0 20.0 > price2017, st30_ref.core 

quantile(dataS1sd_max$SNPs, probs = 0.95)
# 95%
# 25.75 > price2017, st22_ref.whole 
#  26 >  price2017, st30_ref.whole 
# 12.75 > price2017, st22_ref.core 
# 12.75 > price2017, st30_ref.core 

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
plot_file = paste("empirical_clould_of_diversity.allCCs.",output_sufix,".price2017.pdf",sep="");
ggsave(plot_file, plot = g1, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")


##########################################################################################################
###                                   6. APPLYING LINEAR MIXED MODELS                                 ####
##########################################################################################################

# Linear mixed models are applied to calculate the SNP accumulation rate and to model the "cloud of diversity"
# The number of SNPs between MRSA isolates (SNPs) is modelled as a function of the time gap (TimeGap) between isolates
# The intercept (that is, number of SNPs at time 0) is interpreted as the "cloud of diversity" at time 0 and assumbed to 
# vary by patient (AnonPtNo, random variable)
# Since most of MRSA isolates belong to CC22, the linear mixed model was applied to CC22 isolates only
# as well as to all MRSA CCs too

library(lme4)
lmer_all  =  lmer(SNPs  ~ TimeGap   + (1|AnonPtNo), data=dataS1)
summary(lmer_all)


##########################################################################################################
###                                 7. CALCULATION OF THE MRSA SUBSTITUTION RATE                      ####
##########################################################################################################

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

### Using MRSA isolates from all CCs
timegap_coefficient_all = coef(summary(lmer_all))[2,1]
timegap_coefficient_all
# [1] 0.01166199  > price2017, st22_ref.whole 
# [1] 0.006505486 > price2017, st22_ref.core 
# [1] 0.01508212  > price2017, st30_ref.whole 
# [1] 0.006520493 > price2017, st30_ref.core 

## Converting units to SNPs per genome per year
substitution_rate_all = timegap_coefficient_all*365
substitution_rate_all
# [1] 4.256625 > price2017, st22_ref.whole 
# [1] 2.374502 > price2017, st22_ref.core 
# [1] 5.504972 > price2017, st30_ref.whole 
# [1] 2.37998 > price2017, st30_ref.core 

## 95% confidence interval
CI = confint(lmer_all, "TimeGap", level = 0.95)
CI_lower_bound = CI[1]*365
CI_lower_bound
# [1] 3.758648 > price2017, st22_ref.whole 
# [1] 2.107302 > price2017, st22_ref.core 
# [1] 4.977964 > price2017, st30_ref.whole 
# [1] 2.110005 > price2017, st30_ref.core 

CI_upper_bound = CI[2]*365
CI_upper_bound
# [1] 4.753981 > price2017, st22_ref.whole 
# [1] 2.641723 > price2017, st22_ref.core 
# [1] 6.031175 > price2017, st30_ref.whole 
# [1] 2.649974 > price2017, st30_ref.core 

## Converting units to substitutions per site per year
substitution_rate_all_ps = substitution_rate_all/(chromosome_length-mge_regions_length);
substitution_rate_all_ps
# [1] 1.59235e-06 > price2017, st22_ref.whole 
# [1] 1.351664e-06 > price2017, st22_ref.core 
# [1] 2.126999e-06 > price2017, st30_ref.whole 
# [1] 1.359707e-06 > price2017, st30_ref.core 

CI_lower_bound_ps = CI_lower_bound/(chromosome_length-mge_regions_length)
CI_lower_bound_ps
# [1] 1.406063e-06 > price2017, st22_ref.whole 
# [1] 1.199563e-06 > price2017, st22_ref.core 
# [1] 1.923375e-06 > price2017, st30_ref.whole 
# [1] 1.205468e-06 > price2017, st30_ref.core 

CI_upper_bound_ps = CI_upper_bound/(chromosome_length-mge_regions_length)
CI_upper_bound_ps
# [1] 1.778404e-06 > price2017, st22_ref.whole 
# [1] 1.503777e-06 > price2017, st22_ref.core 
# [1] 2.330312e-06 > price2017, st30_ref.whole 
# [1] 1.513958e-06 > price2017, st30_ref.core 


##########################################################################################################
###                      APPLYING LINEAR MIXED MODELS WITH SUB-SAMPLED DATA                          ####
##########################################################################################################

dim(dataS1)
# [1] 11167    17

library(lme4)
number_iterations = 100
mutation_rates_iterations = vector()
mutation_rates_lci_iterations = vector()
mutation_rates_uci_iterations = vector()
beta_zero_iterations = vector()

host_ids = unique(as.vector(dataS1$AnonymisedPatientId))
length(host_ids)
# [1] 255

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
      hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate1 == host_collection_dates[1] & dataS1$ClonalComplex1 == host_clonal_complex)
      if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(dataS1$SequencingTag1[hhhd1])); }
      hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate2 == host_collection_dates[1] & dataS1$ClonalComplex2 == host_clonal_complex)
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
        hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate1 == host_collection_dates[d] & dataS1$ClonalComplex1 == host_clonal_complex)
        if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(dataS1$SequencingTag1[hhhd1])); }
        hhhd1 = which(dataS1$AnonymisedPatientId == host_ids[h] & dataS1$CollectionDate2 == host_collection_dates[d] & dataS1$ClonalComplex2 == host_clonal_complex)
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
# 4.130730 4.660498 4.872606 5.079779 5.632804  > st22_ref.whole
# 5.074884 5.966094 6.237340 6.455723 7.072336  > st30_ref.whole
# 2.373727 2.669244 2.839473 2.939753 3.422637  >  st22_ref.core
# 2.327135 2.736639 2.828565 2.951420 3.280797  >  st30_ref.core

quantile(beta_zero_iterations)
# 0%      25%      50%      75%     100% 
#  9.167010  9.710830  9.948941 10.133470 10.858091 > st22_ref.whole
#  9.55004 10.92350 11.54625 11.85050 12.54120   > st30_ref.whole
#  2.177967 2.881564 3.143878 3.428547 4.946904   > st22_ref.core
# 2.094703 2.831938 3.010641 3.265376 4.450402  > st30_ref.core

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
# [1] "0.0134454107239625 (0.0116198721905772 - 0.0152666460326781)" > st22_ref.whole
# [1] "0.0173116453711501 (0.0149372750590745 - 0.019679732847411)" > st30_ref.whole
# [1] "0.00786090180552467 (0.00685597490722515 - 0.00886527411788905)" > st22_ref.core
# [1] "0.00764851295449261 (0.00679801395462071 - 0.00849821006350058)" > st30_ref.core

# Converting units to substitutions per site per year
print(paste(mutation_rate_median*365," (",mutation_rate_lci*365," - ",mutation_rate_uci*365,")",sep=""))
# [1] "4.90757491424633 (4.24125334956067 - 5.57232580192749)"  > st22_ref.whole
# [1] "6.31875056046978 (5.45210539656218 - 7.18310248930501)" > st30_ref.whole
# [1] "2.8692291590165 (2.50243084113718 - 3.2358250530295)" > st22_ref.core
# [1] "2.7917072283898 (2.48127509343656 - 3.10184667317771)" > st30_ref.core

print(paste(mutation_rate_median*365/(chromosome_length-mge_regions_length)," (",mutation_rate_lci*365/(chromosome_length-mge_regions_length)," - ",mutation_rate_uci*365/(chromosome_length-mge_regions_length),")",sep=""))
# [1] "1.83586200747514e-06 (1.58659949661326e-06 - 2.08453694783856e-06)" > st22_ref.whole
# [1] "2.44142533265966e-06 (2.10657282703493e-06 - 2.77539178301986e-06)" > st30_ref.whole
# [1] "1.6332830460183e-06 (1.42448638297809e-06 - 1.84196448108241e-06)" > st22_ref.core
# [1] "1.5949313504234e-06 (1.41757824577805e-06 - 1.77211723813572e-06)" > st30_ref.core

