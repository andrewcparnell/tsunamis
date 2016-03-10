# A JAGS model for a constrained radiocarbon chronology for the tsunami data

# Clear the workspace if necessary
rm(list=ls())

# Set the wd and load in packages
setwd("~/GitHub/tsunamis")
library(rjags)
library(ggplot2)
library(Bchron)

### First sort data
cal = read.table('intcal13.14c',sep=',')
colnames(cal) = c('calbp','c14age','errorBP','Delta14C','Sigma_per_mil') 

# Give data
x_T12 = c(2822,2725) # T12 is first two dates - theta[1]
sig_T12 = c(20,20) 
x_T11 = c(2965,3065,3093,3210,3217) # T11 is 5 dates - theta[2]
sig_T11 = c(20,20,21,21,21)
x_T10 = c(3085,3069) # T10 has 2 dates - theta[3]
sig_T10 = c(21,20)
x_T9 = c(3078) # - theta[4]
sig_T9 = c(21)
x_T8 = c(3077) # - theta[5]
sig_T8 = c(20)
x_T7_max = c(3540) # max for theta[6]
sig_T7_max = c(30)
x_T6_min = c(4742)# - theta[7]
sig_T6_min = c(23)
# No dates for T5 - theta[8]
# No dates for T4 - theta[9]
x_T4_min_max = c(5090) # This is between T3 and T4 and corresponds to theta[10]
sig_T4_min_max = c(40)
# No dates for T3 - theta[11]
# No dates for T2 - theta[12]
x_T2_min_max = c(6060) # This is between T1 and T2 and corresponds to theta[13]
sig_T2_min_max = c(40)
# No dates for T1 - theta[14]
x_T1_max_1 = c(6788) # - still theta[14] - max age for T1 - calibrated this is 7585 to 7675
sig_T1_max_1 = c(26)
x_T1_max_2 = c(6560) # - still theta[14] - another max age for T1
sig_T1_max_2 = c(35)

o = order(cal$calbp) # - JAGS needs the calibration curve to be in order
data=list(c14=cal$c14age[o],calbp=cal$calbp[o],err=cal$errorBP[o],
          x_T12=x_T12,x_T11=x_T11,x_T10=x_T10,x_T9=x_T9,x_T8=x_T8,
          x_T7_max=x_T7_max,x_T6_min=x_T6_min,x_T4_min_max=x_T4_min_max,
          x_T2_min_max=x_T2_min_max,x_T1_max_1=x_T1_max_1,x_T1_max_2=x_T1_max_2,
          sig_T12=sig_T12,sig_T11=sig_T11,sig_T10=sig_T10,sig_T9=sig_T9,sig_T8=sig_T8,
          sig_T7_max=sig_T7_max,sig_T6_min=sig_T6_min,sig_T4_min_max=sig_T4_min_max,
          sig_T2_min_max=sig_T2_min_max,sig_T1_max_1=sig_T1_max_2,
          sig_T1_max_2=sig_T1_max_2,n_T12=length(x_T12),n_T11=length(x_T11),
          n_T10=length(x_T10),n_T9=length(x_T9),n_T8=length(x_T8))

# To get goo starting values, calibrate some radiocarbon dates and sort them
# Give some initial values
source('Initial_values_function_20150624.R')
  
# Model definition
modelstring ='
model{
  # Priors for the ages, make sure they are sorted
  for(i in 1:14) {
    # theta[1]=T12,theta[2]=T11,theta[3]=T10,theta[4]=T9,theta[5]=T8,theta[6]=T7,theta[7]=T6,theta[8]=T5,theta[9]=T4,theta[10]=date between T3 and T4,theta[11]=T3,theta[12]=T2,theta[13]=date between 1 and 2,theta[14]=T1
    theta_raw[i] ~ dunif(0,10000)
  }
  theta[1:14] <- sort(theta_raw)

  # First do T12 - top of core, youngest - Note this is theta[1]
  mu_cal_T12 <- interp.lin(theta[1],calbp,c14)
  sig_cal_T12 <- interp.lin(theta[1],calbp,err)
  for(i in 1:n_T12) {
      x_T12[i] ~ dnorm(mu_cal_T12,tau_all_T12[i])
      tau_all_T12[i] <- 1/sig_sq_all_T12[i]
      sig_sq_all_T12[i] <- pow(sig_T12[i],2) + pow(sig_cal_T12,2)
  }

  # Now T11 - must be older (bigger) than T12 but younger (smaller) than T10
  mu_cal_T11 <- interp.lin(theta[2],calbp,c14)
  sig_cal_T11 <- interp.lin(theta[2],calbp,err)
  for(i in 1:n_T11) {
      x_T11[i] ~ dnorm(mu_cal_T11,tau_all_T11[i])
      tau_all_T11[i] <- 1/sig_sq_all_T11[i]
      sig_sq_all_T11[i] <- pow(sig_T11[i],2) + pow(sig_cal_T11,2)
  }

  # Now T10 - must be older (bigger) than T11
  mu_cal_T10 <- interp.lin(theta[3],calbp,c14)
  sig_cal_T10 <- interp.lin(theta[3],calbp,err)
  for(i in 1:n_T10) {
      x_T10[i] ~ dnorm(mu_cal_T10,tau_all_T10[i])
      tau_all_T10[i] <- 1/sig_sq_all_T10[i]
      sig_sq_all_T10[i] <- pow(sig_T10[i],2) + pow(sig_cal_T10,2)
  }

  # Now T9
  mu_cal_T9 <- interp.lin(theta[4],calbp,c14)
  sig_cal_T9 <- interp.lin(theta[4],calbp,err)
  for(i in 1:n_T9) {
      x_T9[i] ~ dnorm(mu_cal_T9,tau_all_T9[i])
      tau_all_T9[i] <- 1/sig_sq_all_T9[i]
      sig_sq_all_T9[i] <- pow(sig_T9[i],2) + pow(sig_cal_T9,2)
  }

  # T8
  mu_cal_T8 <- interp.lin(theta[5],calbp,c14)
  sig_cal_T8 <- interp.lin(theta[5],calbp,err)
  for(i in 1:n_T8) {
      x_T8[i] ~ dnorm(mu_cal_T8,tau_all_T8[i])
      tau_all_T8[i] <- 1/sig_sq_all_T8[i]
      sig_sq_all_T8[i] <- pow(sig_T8[i],2) + pow(sig_cal_T8,2)
  }

  # T7 max - dates which are older than theta[6]
  theta_T7_max <- theta[6] + extra_T7
  extra_T7 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_T7_max <- interp.lin(theta_T7_max,calbp,c14)
  sig_cal_T7_max <- interp.lin(theta_T7_max,calbp,err)
  x_T7_max ~ dnorm(mu_cal_T7_max,tau_all_T7_max)
  tau_all_T7_max <- 1/sig_sq_all_T7_max
  sig_sq_all_T7_max <- pow(sig_T7_max,2) + pow(sig_cal_T7_max,2)

  # T6 min - date which is younger than T6
  theta_T6_min <- theta[7] - extra_T6
  extra_T6 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_T6_min <- interp.lin(theta_T6_min,calbp,c14)
  sig_cal_T6_min <- interp.lin(theta_T6_min,calbp,err)
  x_T6_min ~ dnorm(mu_cal_T6_min,tau_all_T6_min)
  tau_all_T6_min <- 1/sig_sq_all_T6_min
  sig_sq_all_T6_min <- pow(sig_T6_min,2) + pow(sig_cal_T6_min,2)

  # No dates at all for T5 - theta[8]

  # No dates for T4 - theta[9]

  # This is a date which is between T3 and T4
  mu_cal_T4_min_max <- interp.lin(theta[10],calbp,c14)
  sig_cal_T4_min_max <- interp.lin(theta[10],calbp,err)
  x_T4_min_max ~ dnorm(mu_cal_T4_min_max,tau_all_T4_min_max)
  tau_all_T4_min_max <- 1/sig_sq_all_T4_min_max
  sig_sq_all_T4_min_max <- pow(sig_T4_min_max,2) + pow(sig_cal_T4_min_max,2)

  # No dates for T3 - theta[11]
  
  # No dates for T2 - theta[12]

  # This is a date between T1 and T2 - theta[13]
  mu_cal_T2_min_max <- interp.lin(theta[13],calbp,c14)
  sig_cal_T2_min_max <- interp.lin(theta[13],calbp,err)
  x_T2_min_max ~ dnorm(mu_cal_T2_min_max,tau_all_T2_min_max)
  tau_all_T2_min_max <- 1/sig_sq_all_T2_min_max
  sig_sq_all_T2_min_max <- pow(sig_T2_min_max,2) + pow(sig_cal_T2_min_max,2)

  # T1 max - first date which is older than T1
  theta_T1_max_1 <- theta[14] + extra_T1_1
  extra_T1_1 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_T1_max_1 <- interp.lin(theta_T1_max_1,calbp,c14)
  sig_cal_T1_max_1 <- interp.lin(theta_T1_max_1,calbp,err)
  x_T1_max_1 ~ dnorm(mu_cal_T1_max_1,tau_all_T1_max_1)
  tau_all_T1_max_1 <- 1/sig_sq_all_T1_max_1
  sig_sq_all_T1_max_1 <- pow(sig_T1_max_1,2) + pow(sig_cal_T1_max_1,2)

  # T1 max_2 - first date which is older than T1
  theta_T1_max_2 <- theta[14] + extra_T1_2
  extra_T1_2 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_T1_max_2 <- interp.lin(theta_T1_max_2,calbp,c14)
  sig_cal_T1_max_2 <- interp.lin(theta_T1_max_2,calbp,err)
  x_T1_max_2 ~ dnorm(mu_cal_T1_max_2,tau_all_T1_max_2)
  tau_all_T1_max_2 <- 1/sig_sq_all_T1_max_2
  sig_sq_all_T1_max_2 <- pow(sig_T1_max_2,2) + pow(sig_cal_T1_max_2,2)

}'

# Run the model in jags
model=jags.model(textConnection(modelstring), data=data, inits=init,n.chains=4)
update(model,n.iter=10000) # Burnin period
output=coda.samples(model=model,
                    variable.names=c("theta",
                                     "extra_T7",
                                     "extra_T6",
                                     "extra_T1_1",
                                     "extra_T1_2"), 
                    n.iter=1e7, 
                    thin=5e3)
save(output,file='output_20150624.rda')
#load('output_20150624.rda')
par(mar=c(2,2,2,2))
plot(output)
summary(output)

##########################################################################################

## Extract the parameters - burnin needs to be after about 400,000
## Equivalent to the 800th iteration
output_2 = rbind(output[[1]][800:2000,],output[[2]][800:2000,],output[[3]][800:2000,],output[[4]][800:2000,])

# Now create upper and lower CIs
CIs = t(round(apply(output_2,2,'quantile',c(0.025,0.5,0.975)),0))

##########################################################################################

# First plot - tsunamis on y axis with CIs on x axis
Tsunami = factor(paste('T',12:1,sep=''),levels=paste('T',1:12,sep=''))

df_1 = data.frame(Tsunami,CIs=CIs[4+c(1:9,11:12,14),])
rownames(df_1)=NULL
colnames(df_1)=c('Tsunami','pc2.5','pc50','pc97.5')
  
ggplot(df_1,aes(x=pc50,y=Tsunami,colour=Tsunami))+ geom_point(size=2) + geom_errorbarh(aes(xmax = pc97.5, xmin =pc2.5),height=0.5)+theme_bw()+theme(legend.position='None',axis.title.y=element_text(angle=0,vjust=1,hjust=0))+ggtitle('Age of Tsunamis with 95% error bounds\n')+xlab('Age (thousands of years')+ scale_x_continuous(breaks=seq(2500,7600,by=500),limits=c(2500,7600))

##########################################################################################

# An alternative version - plot densities but facet by Tsunami
output_3 = c(output_2[,'theta[1]'],output_2[,'theta[2]'],output_2[,'theta[3]'],output_2[,'theta[4]'],output_2[,'theta[5]'],output_2[,'theta[6]'],output_2[,'theta[7]'],output_2[,'theta[8]'],output_2[,'theta[9]'],output_2[,'theta[11]'],output_2[,'theta[12]'],output_2[,'theta[14]'])
Tsunami = rep(factor(paste('T',12:1,sep=''),levels=paste('T',12:1,sep='')),each=4804)
df_2 = data.frame(Tsunami=Tsunami,Age=output_3)

ggplot(df_2,aes(x=Age,fill=Tsunami))+ geom_density(colour=NA) + facet_grid(Tsunami ~ .,scales='free')+ theme_bw()+theme(legend.position='None',axis.title.y=element_text(angle=0,vjust=1,hjust=0),axis.text.y = element_blank(),axis.ticks.y = element_blank())+ggtitle('Age of tsunamis\n')+xlab('Age (years BP)')+ scale_x_continuous(breaks=seq(7500,2500,by=-500),limits=c(2500,7600))+ylab("")+coord_trans(x="reverse", y="reverse")+theme(strip.text.y = element_text(size = 8, angle = 0),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
ggsave('Tsunami_fig1.pdf',height=8,width=8)

##########################################################################################

# Create a csv file of ages for Ben
CIs_2 = t(round(apply(output_2,2,'quantile',c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),0))[4+c(1:9,11:12,14),]
rownames(CIs_2) = paste('T',12:1,sep='')
write.csv(CIs_2,file='Tsunami_ages_20150624.csv',quote=FALSE)

##########################################################################################

# Finally calculate Age differences between Tsunamis - need to be careful here
output_4 = output_2[,4+c(1:9,11:12,14)]
colnames(output_4) = paste('T',12:1,sep='')
output_diff = t(apply(output_4,1,'diff'))
colnames(output_diff) = c('T11_12','T10_11','T9_10','T8_9','T7_8','T6_7','T5_6','T4_5','T3_4','T2_3','T1_2')
output_diff_2 = as.vector(output_diff)
Tsunami = rep(factor(c('T11_12','T10_11','T9_10','T8_9','T7_8','T6_7','T5_6','T4_5','T3_4','T2_3','T1_2'),c('T11_12','T10_11','T9_10','T8_9','T7_8','T6_7','T5_6','T4_5','T3_4','T2_3','T1_2')),each=4804)
df_3 = data.frame(Tsunami=Tsunami,Age=output_diff_2)
ggplot(df_3,aes(x=Age,fill=Tsunami))+ geom_density(colour=NA) + facet_grid(Tsunami ~ .,scales='free')+ theme_bw() +theme(legend.position='None',axis.title.y=element_text(angle=0,vjust=1,hjust=0),axis.text.y = element_blank(),axis.ticks.y = element_blank())+ggtitle('Age gaps between consecutive tsunamis\n')+xlab('Age gap') + scale_x_continuous(breaks=seq(0,2600,by=200),limits=c(0,2500))+ylab("")+theme(strip.text.y = element_text(size = 8, angle = 0),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
ggsave('Tsunami_fig2_20150624.pdf',height=8,width=8)

# Create a csv file of age gas for Ben
CIs_3 = t(round(apply(output_diff,2,'quantile',c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),0))
write.csv(CIs_3,file='Tsunami_age_gaps_20150624.csv',quote=FALSE)
