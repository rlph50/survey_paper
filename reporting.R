library(tidyverse)
library(scales)

df100 <- readRDS("sim_100/results_22_40_80.RDS")
be100 <- readRDS('sim_100/bayes_est.rds') %>%
  select(d_mean,f_mean,phi_mean,t) %>%
  rename(bayes_d=d_mean,bayes_f=f_mean,bayes_phi=phi_mean,bayes_t=t)
df100 <- cbind(df100,be100)
df100$len <- 'n = 100'
df100$trial <- 1:1000

df500 <- readRDS("sim_500/results_22_40_80.RDS")
be500 <- readRDS('sim_500/bayes_est.rds') %>%
  select(d_mean,f_mean,phi_mean,t) %>%
  rename(bayes_d=d_mean,bayes_f=f_mean,bayes_phi=phi_mean,bayes_t=t)
df500 <- cbind(df500,be500)
df500$len <- 'n = 500'
df500$trial <- 1:1000

df1000 <- readRDS("sim_1000/results_22_40_80.RDS")
be1000 <- readRDS('sim_1000/bayes_est.rds') %>%
  select(d_mean,f_mean,phi_mean,t) %>%
  rename(bayes_d=d_mean,bayes_f=f_mean,bayes_phi=phi_mean,bayes_t=t)
df1000 <- cbind(df1000,be1000)
df1000$len <- 'n = 1000'
#df1000$css_conv<-df1000$whittle_conv<-df1000$wll_conv<-df1000$qml_conv<-df1000$mb8_conv<-df1000$ar1_conv<-df1000$lw_conv<-0
df1000$trial <- 1:1000

res <- rbind(df100,df500,df1000)
res$len<-factor(res$len,levels=c('n = 100','n = 500','n = 1000'))

rm(df100,df500,df1000,be100,be500,be1000)

be100 <- readRDS('sim_100/bayes_est.rds') %>%
          mutate(len=100,trial=1:1000)
be500 <- readRDS('sim_500/bayes_est.rds') %>%
  mutate(len=500,trial=1:1000)
be1000 <- readRDS('sim_1000/bayes_est.rds') %>%
  mutate(len=1000,trial=1:1000)
bayes_df <- rbind(be100,be500,be1000)
rm(be100,be500,be1000)

resl<-res %>%
  pivot_longer(!c(len,trial),names_to='method_param',values_to='estimate') %>%
  separate(method_param,into=c('method','param'),sep='_') %>%
  mutate(true_value=ifelse(param=='d',0.22,ifelse(param=='f',0.4,ifelse(param=='phi',0.8,0))),
         bias = estimate-true_value) %>%
  filter(param %in% c('d','f','phi')) %>%
  mutate(method=case_when(method=='bayes' ~ 'Bayesian',
                          method=='ar1' ~ 'Arteche\nRobinson',
                          method=='lw'  ~ 'local\nWhittle',
                          method=='css' ~ 'CSS',
                          method=='qml' ~ 'QML',
                          method=='mb8' ~ 'Wavelet\nMB8',
                          method=='whittle' ~ 'Whittle',
                          method=='wll' ~ 'log-Whittle',
                          TRUE ~ 'Other'
                        )) %>%
  mutate(Method=factor(method,levels=c('Bayesian','Arteche\nRobinson','local\nWhittle','CSS','QML','Wavelet\nMB8','Whittle','log-Whittle')))


  
ggplot(resl,aes(x=Method,y=bias)) +
  geom_boxplot(outlier.size=0.3,size=0.2) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_grid(rows=vars(len),cols=vars(param),labeller=labeller(param=label_parsed)) + 
  ylab('Bias') + xlab('') + ggtitle('Comparing distribution of bias in parameter estimates.')+
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave("compare_sims.png",dpi=300,width=8,height=5)

# examine bayes stats
# only 1 with rhat too large
filter(bayes_df,phi_rhat>1.1|d_rhat>1.1|f_rhat>1.1)
# but 133 without sufficient effective samples
filter(bayes_df,phi_neff<400|d_neff<400|f_neff<400,len==1000)
filter(bayes_df,phi_neff<400|d_neff<400|f_neff<400) %>% count()
filter(bayes_df,phi_neff<400|d_neff<400|f_neff<400) %>%
  group_by(len) %>%
  count()

# get time and convergence stats
restc <- res %>%
  pivot_longer(!c(len,trial),names_to='method_param',values_to='estimate') %>%
  separate(method_param,into=c('method','param'),sep='_') %>%
  mutate(method=case_when(method=='bayes' ~ 'Bayesian',
                          method=='ar1' ~ 'Arteche\nRobinson',
                          method=='lw'  ~ 'local\nWhittle',
                          method=='css' ~ 'CSS',
                          method=='qml' ~ 'QML',
                          method=='mb8' ~ 'Wavelet\nMB8',
                          method=='whittle' ~ 'Whittle',
                          method=='wll' ~ 'log-Whittle',
                          TRUE ~ 'Other'
  )) %>%
  mutate(Method=factor(method,levels=c('Bayesian','Arteche\nRobinson','local\nWhittle','CSS','QML','Wavelet\nMB8','Whittle','log-Whittle'))) %>%
  select(len,trial,Method,param,estimate) %>%
  filter(param %in% c('t','conv')) %>%
  pivot_wider(id_cols=c(len,trial,Method),names_from='param',values_from='estimate') %>%
  mutate(med_t=t) %>%
  group_by(len,Method) %>%
  summarise(t=mean(t,na.rm=TRUE),med_t=median(med_t,na.rm=TRUE),conv=sum(conv,na.rm=TRUE))
  
head(restc)

resl %>% 
  group_by(len,Method,param) %>%
  summarise(avg_bias=mean(bias,na.rm=TRUE),
            mse = mean(bias^2,na.rm=TRUE)) %>%
  pivot_wider(names_from='param',values_from=c('avg_bias','mse')) %>%
  ungroup() %>%
  left_join(restc,by=c('len','Method')) %>%
  mutate(label=sprintf("%20.20s %7.4f (%6.4f) %7.4f (%6.4f) %7.4f (%6.4f) %6.2f %6.2f %d\n",
                       Method,avg_bias_d,mse_d,avg_bias_f,mse_f,avg_bias_phi,mse_phi,t,med_t,conv)) %>%
  select(len,label) %>% print(n=90)


res %>%
  select(len,ends_with('_t')) %>%
  pivot_longer(!len,names_to='method_param',values_to='estimate') %>%
  filter(method_param %in% c('whittle_t','css_t','ar1_t','lw_t','qml_t','mb8_t','bayes_t','wll_t'))

rest<-res %>%
  select(len,ends_with('_t')) %>%
  pivot_longer(!len,names_to='method_param',values_to='estimate') %>%
  separate(method_param,into=c('method','param'),sep='_') %>%
  mutate(method=case_when(method=='bayes' ~ 'Bayesian',
                          method=='ar1' ~ 'Arteche\nRobinson',
                          method=='lw'  ~ 'local\nWhittle',
                          method=='css' ~ 'CSS',
                          method=='qml' ~ 'QML',
                          method=='mb8' ~ 'Wavelet\nMB8',
                          method=='whittle' ~ 'Whittle',
                          method=='wll' ~ 'log-Whittle',
                          TRUE ~ 'Other'
  )) %>%
  mutate(Method=factor(method,levels=c('Bayesian','Arteche\nRobinson','local\nWhittle','CSS','QML','Wavelet\nMB8','Whittle','log-Whittle'))) %>%
  mutate(n=ifelse(len=='n = 100',100,ifelse(len=='n = 500',500,1000))) %>%
  group_by(n,Method) %>%
  summarise(avg_secs=mean(estimate,na.rm=TRUE))

arrange(rest,Method,n) %>%
print(n=90)

ggplot(rest,aes(x=n,y=avg_secs,linetype=Method)) +
  geom_line() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()
