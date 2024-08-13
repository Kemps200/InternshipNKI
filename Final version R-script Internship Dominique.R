### FINAL VERSION R-SCRIPT INTERNSHIP DOMINIQUE KEMPS ###

#Set working directory
setwd()

#Load necessary packages
library(ggplot2)
library(openxlsx)
library(tidyr)
library(reshape2)
library(gganimate)
library(gifski)
library(av)
library(tidyverse)
library(svglite)
library(car)
library(dplyr)
library(ggpubr)

#Load Data (either as XLSX or CSV) I used xlsx files
Data_WT <- read.xlsx("Control_PGE1.xlsx")
Data_KO <- read.xlsx("AZGP1_KO_PGE1.xlsx")

Data_WT <- read.csv("Control_rep1.csv")
df_tidy1 <- gather(Data_WT, "Cell", "Lifetime_.ns.", -`time_.s.`)
df_tidy<- na.omit(df_tidy1)
df_summary<- data.frame(Time=Data_WT$`time_.s.`, n=tapply(df_tidy$`Lifetime_.ns.`, df_tidy$`time_.s.`, length), mean=tapply(df_tidy$`Lifetime_.ns.`, df_tidy$`time_.s.`, mean))


#Put the xlsx data in 'tidy' format
df_tidy <- gather(Data_WT, "Cell", "Lifetime_(ns)", -`time_(s)`)
df_tidy_KO <- gather(Data_KO, "Cell", "Lifetime_(ns)", -`time_(s)`)

#Make a summary dataframe from xlsx tidy data and calculate mean lifetime trace
df_summary<- data.frame(Time=Data_WT$`time_(s)`, n=tapply(df_tidy$`Lifetime_(ns)`, df_tidy$`time_(s)`, length), mean=tapply(df_tidy$`Lifetime_(ns)`, df_tidy$`time_(s)`, mean))
df_summary_KO<- data.frame(Time=Data_KO$`time_(s)`, n=tapply(df_tidy_KO$`Lifetime_(ns)`, df_tidy_KO$`time_(s)`, length), mean=tapply(df_tidy_KO$`Lifetime_(ns)`, df_tidy_KO$`time_(s)`, mean))


#Provide constant values (if necessary) #This needs to be adjusted according to the sample
stim_isop<-370.959
stim_forsIBMX<-2596.716

######Calculate CI and sd + sem values#####   didn't do this for the graphs in the report only for statistical testing
Conf_level= 0.95

df_summary$sd <- tapply(df_tidy$`Lifetime_(ns)`, df_tidy$`time_(s)`, sd)
df_summary_KO$sd <- tapply(df_tidy_KO$`Lifetime_(ns)`, df_tidy_KO$`time_(s)`, sd)
mean_sd<- mean(df_summary$sd)
mean_sd_KO<- mean(df_summary_KO$sd)

df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)

df_summary$CI_lower <- df_summary$mean + qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem

df_summary$CI_upper <- df_summary$mean - qt((1-Conf_level)/2, df=df_summary$n-1)*df_summary$sem



## Plot all observations in one plot for in the supplementals (we call this plot 'p')   Make sure to change plot names and data sets if necessary ##
p<- ggplot(df_summary, aes(x=Time, y=mean)) +
  geom_line(data=df_tidy, aes(x=`time_(s)`, y=`Lifetime_(ns)`, group=Cell), color="grey") +
  geom_line(data=df_summary, aes(x=Time, y=mean), linewidth=2, alpha=0.8, color="black") +
  labs(title = "HeLa H250 Wildtype Control HITS", x = "Time (s)", y = "Lifetime (ns)") +
  theme_classic() +
  #geom_vline(xintercept=stim_isop, linetype="dashed", color="black") +
  #geom_vline(xintercept=stim_forsIBMX, linetype="dashed", color="black") +
  #annotate("text", x=stim_isop, y=min(df_tidy$`Lifetime_(ns)`), label="10nM isoproterenol", vjust=-1.5, hjust=-0.05, color="black", size=2.9)+
  #annotate("text", x=stim_forsIBMX, y=min(df_tidy$`Lifetime_(ns)`), label="100µM IBMX\n 25µM and F", vjust=-0.45, hjust=-0.05, color="black", size=2.9)+ 
  scale_y_continuous(limits=c(1.5, 3.5)) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Customize title appearance
    axis.title = element_text(size = 12),  # Customize axis label appearance
    axis.text = element_text(size = 10))

#Save this graph
base_filename<- "HeLa H250 wildtype HITS plot of ALL traces"

ggsave(filename=paste0(base_filename, ".svg"), plot=p, width=7, height=5)
ggsave(filename=paste0(base_filename, ".png"), plot=p, width=7, height=5)
ggsave(filename=paste0(base_filename, ".pdf"), plot=p, width=7, height=5)

#Animate the plot
animated_line_plot <- p +
  transition_reveal(`time_(s)`) +
  ease_aes('linear')

# Display the animated line plot
animate(p, nframes = 100, fps = 10, width = 600, height = 400, renderer = gifski_renderer())

# Save the animation
anim_save("animated_control.gif", animation = last_animation())



## Select 20 random lines from the data that will be plotted for in the report ##
unique_traces<- unique(df_tidy$Cell)
selected_traces<- sample(unique_traces, 20)
selected_data<- df_tidy %>% filter(Cell %in% selected_traces)

#Plot the selected traces with different colors #For xlsx files
plot_20_selected<- ggplot(selected_data, aes(x=`time_(s)`, y=`Lifetime_(ns)`, group=`Cell`, color=`Cell`)) + geom_line() + 
  geom_line(data=df_summary, aes(x=Time, y=mean), linewidth=2, alpha=0.8, color="black", inherit.aes = FALSE)+
  labs(title = "HeLa H250 VPS26A KO Rep3", x = "Time (s)", y = "Lifetime (ns)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Customize title appearance
    axis.title = element_text(size = 13),  # Customize axis label appearance
    axis.text = element_text(size = 12)) +
  geom_vline(xintercept=stim_isop, linetype="dashed", color="black") +
  geom_vline(xintercept=stim_forsIBMX, linetype="dashed", color="black") +
  #annotate("text", x=stim_isop, y=min(df_tidy$`Lifetime_(ns)`), label="1µM IP", vjust = 0.5, hjust=0.5, color="black", size=2.9)+
#annotate("text", x=stim_forsIBMX, y=min(df_tidy$`Lifetime_(ns)`), label="100µM IBMX\n 25µM F", vjust= 1.0, hjust=-0.05, color="black", size=2.9)+ 
  scale_color_manual(values=scales::hue_pal()(20)) + theme(legend.position="none") + ylim(2.0,3.4) +
  theme(plot.margin=unit(c(1,1,2,2), "cm"))

#For csv files
plot_20_selected<- ggplot(selected_data, aes(x=`time_.s.`, y=`Lifetime_.ns.`, group=`Cell`, color=`Cell`)) + geom_line() +
  geom_line(data=df_summary, aes(x=Time, y=mean), linewidth=2, alpha=0.8, color="black", inherit.aes = FALSE)+
  labs(title = "HeLa H250 Wildtype", x = "Time (s)", y = "Lifetime (ns)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Customize title appearance
    axis.title = element_text(size = 13),  # Customize axis label appearance
    axis.text = element_text(size = 12)) +
  geom_vline(xintercept=stim_isop, linetype="dashed", color="black") +
  geom_vline(xintercept=stim_forsIBMX, linetype="dashed", color="black") +
  annotate("text", x=stim_isop, y=min(df_tidy$`Lifetime_(ns)`), label="1µM IP", vjust= 3, hjust=-0.05, color="black", size=2.9)+
  annotate("text", x=stim_forsIBMX, y=min(df_tidy$`Lifetime_(ns)`), label="100µM IBMX\n 25µM F", vjust= 1.0, hjust=-0.05, color="black", size=2.9)+ 
  scale_color_manual(values=scales::hue_pal()(200)) + theme(legend.position="none") + ylim(2.0,3.4) +
  theme(plot.margin=unit(c(1,1,2,2), "cm"))


#Save the plot
base_filename<- "HeLa H250 VPS26A rep3 of 20 selected traces bigger font"


ggsave(filename=paste0(base_filename, ".svg"), plot=plot_20_selected, width=7, height=5)
ggsave(filename=paste0(base_filename, ".png"), plot=plot_20_selected, width=7, height=5)
ggsave(filename=paste0(base_filename, ".pdf"), plot=plot_20_selected, width=7, height=5)


### Making a Box-plot ###
# Make a selection of the data from your df that you want to plot in the boxplot  In this case I want to plot two datapoints for each cell
df_boxplot <- df_tidy %>%
  group_by(Cell) %>%    # Group by 'Cell'
  filter(row_number() %in% c(7, 9)) %>%  # Selects the 7th and 9th row for each Cell
  ungroup()          

#df_boxplot <- tail(df, Frames_to_take)
df_boxplot <- data.frame(colMeans(df_boxplot[,-1]))
colnames(df_boxplot)[1] <- "X"

write.table(df_boxplot, "D:/Projects/LITE/REVISION/Data_exported_again_for_paper_revision/Boxplot_mean_data_extracted_for_every_plot/fdFLIM/Cos7H250_transient_fdFLIM_NDFIlterFigure3D_1000s.csv", sep =";" ,row.names=TRUE)   #new piece of code entered on 10-01-2024 to save the mean of the boxplot points

#Plotting the boxplot with all data points

boxplot = ggplot(df_boxplot, aes(x = "Lifetime_(ns)", y = X)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.3,              #width of errorbar
               alpha = 1,                 #Transparency of error bar
               linewidth = 0.6) +              #Size of error bar
  geom_boxplot(width = 0.6,                 #width of the boxplot    
               alpha = 1,                 #Transparency of boxplot
               color = "black",           #Border color of boxplot
               #outlier.color = 3
               outlier.shape = NA,         #Outlier color  
               lwd = 0.6, 
               fatten = 0.6)            +   #line width of boxplot(the box)
  coord_cartesian(ylim = c(2.25, 3.4)) +
  geom_point(position = position_jitter(seed = 0.5, width = 0.15), size =0.2, color ="#FF7F00") +
  theme_void(base_size = 16) + theme(aspect.ratio =5.5)

show(boxplot)

##### Statistical testing #####

###Test for normal distribution
#Histogram
p<- ggplot(subset_VPS26A, aes(x = `Lifetime_(ns)`)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = 'blue', alpha = 0.5) +
  geom_density(color = 'black') + 
  theme_minimal() +
  scale_x_continuous(limits=c(1.5,3.5)) +
  ggtitle("Histogram and Density Plot HeLa H250 gRNA population")

base_filename<- "HeLa H250 gRNA histogram"

ggsave(filename=paste0(base_filename, ".svg"), plot=p, width=7, height=5)
ggsave(filename=paste0(base_filename, ".png"), plot=p, width=7, height=5)
ggsave(filename=paste0(base_filename, ".pdf"), plot=p, width=7, height=5)

#QQ plot
p<- ggqqplot(subset_control$`Lifetime_(ns)`, title = "Q-Q Plot HeLa H250 gRNA population")

base_filename<- "HeLa H250 gRNA qqplot"

ggsave(filename=paste0(base_filename, ".svg"), plot=p, width=7, height=5)
ggsave(filename=paste0(base_filename, ".png"), plot=p, width=7, height=5)
ggsave(filename=paste0(base_filename, ".pdf"), plot=p, width=7, height=5)


### Test for equal variances
var.test(df_tidy$`Lifetime_(ns)`, df_tidy_KO$`Lifetime_(ns)`)

result1 <- var.test(df_FB$`Lifetime_(ns)`, df_FB_Hepes$`Lifetime_(ns)`)
result2 <- var.test(df_FB$`Lifetime_(ns)`, df_HCG$`Lifetime_(ns)`)
result3 <- var.test(df_FB_Hepes$`Lifetime_(ns)`, df_HCG$`Lifetime_(ns)`)
print(result3)


#Kruskal_wallis test (non-parametric alternative to one-way ANOVA)
#Non-normal distribution
lambda<- mean(df_tidy_low$`Lifetime_(ns)`)
df_tidy_low<- rpois(6970*3, lambda = lambda)
summary(df_tidy_low)

df_test_control<- rpois(1313,lambda=mean_sd)
df_test_KO<-rpois(798, lambda=mean_sd_KO)

#normal distribution
mean_WT<- mean(df_tidy_baseline$`Lifetime_(ns)`)
mean_High<- mean(df_tidy_high_baseline$`Lifetime_(ns)`)
df_tidy<- rnorm(5468*3, mean=mean_WT, sd=0.1661961)

#Perform Kruskal_wallis test
kruskal.test(list(df_test_control, df_test_KO))

# Welchs T-test
# Merge data frames by timepoints
df_summary$time_s<- c(0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575)
df_summary_AZGP1$time_s<- c(0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575)

print(df_summary)
merged_df<- merge(df_summary, df_summary_AZGP1, by="time_s", suffixes = c("control","AZGP1"))
print(merged_df)

merged_df$sdcontrol <- as.numeric(as.character(merged_df$sdcontrol))
merged_df$sdAZGP1 <- as.numeric(as.character(merged_df$sdAZGP1))

str(merged_df)

# Perform Levene's test
levene_test<- leveneTest(merged_df$sdcontrol, merged_df$sdAZGP1)

##Perform Welch test on mean sd values
t_test<- t.test(merged_df$sdcontrol, merged_df$sdAZGP1, var.equal=FALSE)
print(t_test)

##Perform Welch test on individual time points
df <- merged_df %>%
  select(time_s, sdAZGP1, sdcontrol)

# Group by time_s and perform Welch's t-test for sdAZGP1 vs sdcontrol
welch_results <- df %>%
  group_by(`Time`) %>%
  summarize(
    t_test_result = list(t.test(mean_sd, mean_sd_KO, var.equal = FALSE) %>% tidy())
  ) %>%
  ungroup()

#Wilcoxon Rank-Sum Test (Mann_whitney U test) Test is used to compare sd of two independent samples. It is non-parametric and not normally distributed
wilcox.test(mean_sd, mean_sd_KO)
