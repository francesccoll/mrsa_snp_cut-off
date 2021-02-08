# If you find this code useful, please cite:
# Coll et al. Peacock SJ. Definition of a genetic relatedness cutoff to exclude recent transmission of meticillin-resistant Staphylococcus aureus: a genomic epidemiology analysis. Lancet Microbe. 2020 Dec;1(8):e328-e335. doi: 10.1016/S2666-5247(20)30149-X. PMID: 33313577; PMCID: PMC7721685. (https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30149-X/fulltext)

require(gdata)
require(ggplot2)

##########################################################################################################
###                                         1. INPUT FILES                                            ####
##########################################################################################################

working_dir =""
setwd(working_dir)

dataS2_file = paste(working_dir,"SupplementaryData2.xlsx",sep="")
dataS2 = read.xls(dataS2_file, sheet = 1, header = T)
dim(dataS2)
# [1] 779   7

dataS2 = dataS2[-which(dataS2$Epi_link_strength=="Same patient"),]
dim(dataS2)
# [1] 777   7

##########################################################################################################
###                           2. KEEPING SUBSET OF PATIENTS WITH EPI DATA                             ####
##########################################################################################################

#### Keeping entries with <= 50 SNPs
dataS2 = dataS2[which(dataS2$SNP_dis_st22 <= 50),];
dim(dataS2)
# [1] 429   7

#### Removing patients for which epi data could not be extracted (labelled as Unknown in Epi.link.strength)
dataS2 = dataS2[-which(dataS2$Epi_link_strength=="Unknown"),];
dim(dataS2)
# [1] 294   7



##########################################################################################################
###                                        3. PLOTTING THE DATA                                       ####
##########################################################################################################

##### whole-genome distances

# Plotting strength of epidemiological links with increasing SNP distances
size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
axis_text_size = 15; axis_title_size = 20; ann_text_size = 5;
g1 <- ggplot(dataS2, aes(x=SNP_dis_st22, fill = Epi_link_strength)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_bar() + 
  xlim(0, 50) + 
  ylim(0, 40) + 
  xlab("Number of SNPs") +
  ylab("Number of Patients") +
  ggtitle("Strength of Epidemiological Link") +
  theme(text = element_text(family = font)) +
  theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
plot_file = "Figure2A.pdf";
ggsave(plot_file, plot = g1, device = "pdf",width = 8, height = 5, dpi = 300, units = "in")


# Plotting strength of epidemiological links with increasing SNP distances -  core-genome distances
size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
axis_text_size = 15; axis_title_size = 20; ann_text_size = 5;
g1 <- ggplot(dataS2, aes(x=SNP_dis_st22_core, fill = Epi_link_strength)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_bar() + 
  ylim(0, 40) + 
  xlim(0, 50) + 
  xlab("Number of SNPs") +
  ylab("Number of Patients") +
  ggtitle("Strength of Epidemiological Link") +
  theme(text = element_text(family = font)) +
  theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
plot_file = "Figure2B.pdf";
ggsave(plot_file, plot = g1, device = "pdf",width = 8, height = 5, dpi = 300, units = "in")



# Plotting distribution of SNP distances among patients with Strong epi links

subset = "all"; # use this to study distribution of SNP distances among patients with strong epi links
subset = "hospital"; # use this to study distribution of SNP distances among patients with strong hospital epi links
subset = "community"; # use this to study distribution of SNP distances among patients with strong community epi links

plot_width = 6; plot_height = 5;
plot_between_host_diversity = function(data, text_x_offset, plot_title)
{
  size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
  axis_text_size = 15; axis_title_size = 20; ann_text_size = 5;
  co_x = round((nrow(data)/100)*95)
  co_y = as.numeric(data$SNP_dis_st22[co_x])
  g1 <- ggplot(data, aes(x=seq(1,nrow(data),1), y=SNP_dis_st22)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
    geom_point(shape = 21, colour = dot_color, fill = dot_color, size = size_dot) + 
    ylim(0, 50) + 
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
if(subset=="all"){ dataS2 = subset(dataS2, Epi_link_strength=='Strong'); }
if(subset=="hospital"){ dataS2 = subset(dataS2, Epi_link_strength=='Strong' & grepl('Same ward', Epi_link_description)); }
if(subset=="community"){ dataS2 = subset(dataS2, Epi_link_strength=='Strong' & (grepl('GP', Epi_link_description) | grepl('Postcode', Epi_link_description))); }
dataS2 = dataS2[order(as.numeric(dataS2$SNP_dis_st22)),];
plot_title = "Between-Hosts Diversity"
g1 = plot_between_host_diversity(dataS2,text_x_offset,plot_title)
if(subset=="all"){ plot_file = "Figure2C.pdf"; }
if(subset=="hospital"){ plot_file = "SupplementaryFigure1A.pdf"; }
if(subset=="community"){ plot_file = "SupplementaryFigure1B.pdf"; }
ggsave(plot_file, plot = g1, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")



# Plotting distribution of SNP distances among patients with Strong epi links - core-genome distances
plot_width = 6; plot_height = 5;
plot_between_host_diversity = function(data, text_x_offset, plot_title)
{
  size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
  axis_text_size = 15; axis_title_size = 20; ann_text_size = 5;
  co_x = round((nrow(data)/100)*95)
  co_y = as.numeric(data$SNP_dis_st22_core[co_x])
  g1 <- ggplot(data, aes(x=seq(1,nrow(data),1), y=SNP_dis_st22_core)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
    geom_point(shape = 21, colour = dot_color, fill = dot_color, size = size_dot) + 
    ylim(0, 50) + 
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
if(subset=="all"){ dataS2 = subset(dataS2, Epi_link_strength=='Strong'); }
if(subset=="hospital"){ dataS2 = subset(dataS2, Epi_link_strength=='Strong' & grepl('Same ward', Epi_link_description)); }
if(subset=="community"){ dataS2 = subset(dataS2, Epi_link_strength=='Strong' & (grepl('GP', Epi_link_description) | grepl('Postcode', Epi_link_description))); }
dataS2 = dataS2[order(as.numeric(dataS2$SNP_dis_st22_core)),];
plot_title = "Between-Hosts Diversity"
g1 = plot_between_host_diversity(dataS2,text_x_offset,plot_title)
if(subset=="all"){ plot_file = "Figure2D.pdf"; }
if(subset=="hospital"){ plot_file = "SupplementaryFigure1C.pdf"; }
if(subset=="community"){ plot_file = "SupplementaryFigure1D.pdf"; }
ggsave(plot_file, plot = g1, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")





