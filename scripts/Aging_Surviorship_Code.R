#Introduction####
#Title: Kaplan Meier Aging Survivorship

#Purpose: Estimate the survivorship of the flies. 

# USEFUL webpages
#http://reliawiki.org/index.php/Life_Data_Classification
#https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#what_is_censoring

#Created: 11/6/23, SJM
#Last Edited: 8/28/24, KMH

#Kaplan Meier Estimation####

#Directory Structure
#Master Directory/
#   Survivorship/
#     death.matrix.txt

# LOAD library
library(survival)

# READ in death matrix file
# (how many died/sampled per day)
death.matrix <- read.table(file = "Survivorship/death.matrix.txt", header = TRUE)

# GET a vector of the day in which every fly died
days.of.death <- as.numeric(unlist(apply(death.matrix,1,function(XXX) rep(XXX[1],XXX[2]))))
#We are repeating whatever the value in Day is by the number of dead. See use of rep
head(days.of.death)
#[1] 3 3 4 4 4 5

# GET a vector of the day in which every fly was sampled
days.of.censor <- as.numeric(unlist(apply(death.matrix,1,function(XXX) rep(XXX[1],XXX[3]))))
#Same as above except now repeating based on the value in Censored column
head(days.of.censor, n=51)
#[1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6

# MAKE 'days.of.death' into matrix (where 2nd column is all 1's)
matrix.of.death <- cbind(days.of.death,rep(1,length(days.of.death)))

# MAKE 'days.of.censor' into matrix (where 2nd column is all 0's)
matrix.of.censor <- cbind(days.of.censor,rep(0,length(days.of.censor)))

# PULL together the pair of matrices
all.matrix <- rbind(matrix.of.death,matrix.of.censor)

# CONVERT into a 'Surv' object
# COMPUTES survival curve using Kaplan-Meier estimate
my.surv.object <- Surv(all.matrix[,1],all.matrix[,2])
my.fit <- survfit(my.surv.object ~ 1)

# RETURN the Kaplan-Meier estimates per day (where min/max are those
#    in the original 'all.matrix' file)
# RETURN the day vector
KM.estimates <- summary(my.fit)$surv
KM.days <- summary(my.fit)$time
KM.se <- summary(my.fit)$std.err

# PLOT Kaplan-Meier curve
pdf(file="Survivorship/KaplanMeier_Survivorship.pdf",width=6,height=6)
plot(my.fit, main="Kaplan-Meier estimate with 95% confidence bounds",
	xlab="Time (days)", ylab="Survivorship")
graphics.off()

# SPIT out datapoints as a matrix
all.KM.data <- cbind(KM.days,KM.estimates,KM.se)

# MARK those rows conforming to days of sampling for RNA experiment
# ADD to matrix
# (file pasted below)
#Getting the surviorship for our RNAseq samples that were sent in for sequencing
sampling.days <- c(3,6,10,14,17,23,27,31,36,38,42,48,50,55,59)
sampling.days <- 1*(match(all.KM.data[,"KM.days"],sampling.days,nomatch=0)!=0)
all.KM.data <- cbind(all.KM.data,sampling.days)

# WRITE out
write.table(all.KM.data,file="Survivorship/output.KaplanMeier.survivorship.txt",quote=F,sep="\t",col.names=NA)

# SPLIT out the equivalent days/survivorship estimates
# WRITE out
# (file pasted below)
trunc.KM.data <- all.KM.data[all.KM.data[,"sampling.days"]==1,1:2]
colnames(trunc.KM.data) <- c("Day","Survivorship")
write.table(trunc.KM.data,file="Survivorship/A4.days.vs.survivorship.txt",quote=F,sep="\t",col.names=NA)

# PLOT Kaplan-Meier curve with sampling points
pdf(file="KaplanMeier_Survivorship_plusSamplingPoints.pdf",width=6,height=6)
plot(my.fit, main="Kaplan-Meier estimate with 95% confidence bounds",
	xlab="Time (days)", ylab="Survivorship")
points(trunc.KM.data[,1],trunc.KM.data[,2],pch=16,col="red")
graphics.off()

#Sup Figure 1####
all.KM.data <- read.table("Survivorship/output.KaplanMeier.survivorship.txt", header = TRUE)
colnames(all.KM.data) <- c("day", "Survivorship", "SE", "Sampling_Day")
all.KM.data$Add_Day_inf <- all.KM.data$Sampling_Day
all.KM.data$Sampling_Day <- as.factor(all.KM.data$Sampling_Day)
all.KM.data$Add_Day_inf[all.KM.data$day == 17 | all.KM.data$day == 31 | all.KM.data$day == 45] <- 3
all.KM.data$Add_Day_inf <- as.factor(all.KM.data$Add_Day_inf)

SurvivorvDay_plot <- 
  all.KM.data %>% ggplot(aes(x=day, y=Survivorship))+
  geom_line(linewidth = 0.5)+ #of 0.75
  geom_point(aes(color=Sampling_Day), show.legend = TRUE, size=1)+
  scale_color_manual(values = c("black", "red"), labels = c("Count Day","Sample Day"))+
  xlab("Day")+
  ggtitle("Survivorship Curve")+
  labs(color="Sampling Day")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=7), legend.position = c(0.2,0.25))

ggsave("Plots/SupFig1_Survival Curve", SurvivorvDay_plot, width = 4, height =2.5, units = "in")
