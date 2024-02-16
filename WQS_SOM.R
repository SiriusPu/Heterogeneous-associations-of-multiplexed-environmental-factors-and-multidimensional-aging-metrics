# install packages gWQS and kohonen
install.packages('gWQS')
library(gWQS)
install.packages('kohonen')
library('kohonen')

#set path
#UKBdata is the dataset including all covariates, environmental exposures and aging metrcis
#Detailed definition of covariates and calculation of aging metrics have been already provided in the article
setwd("/Users/siriuspu/Desktop")
UKBdata <- read.csv('UKBdata.csv')

#toxic_chems is the list of nine environmental factors
toxic_chems <- names(UKBdata)[16:24]


#run WQS model for each aging metric
#Take PhenoAge as an example
results <- gwqs(PhenoAge~wqs+age+sex+bmi+ethnic+smoking+drinking+diet+exercise+CVD+cancer+nSES,mix_name = toxic_chems,data = UKBdata,q=4,validation = 0.6,b=500,b1_pos = FALSE,b_constr = FALSE,family = "gaussian",
                                          seed=100,plots=TRUE,tables=TRUE)


#SOM
set.seed(2021)
z1 <- UKBdata[,c("TRN_24h","TRN_16h","TRN_daytime","TRN_nighttime","TRN_evening","Nitrogen_oxides_2010","Nitrogen_dioxide_05_10_ave",
             "PM25_2010","PM10_07_10_ave","PM25_10_2010","Greenspace_1000m","Greenspace_300m","water_per_300","water_per_1000")]
z2 <- na.omit(z1)
X2 <- z2[,c('id')]
X <- as.matrix(scale(z2))
X <- na.omit(X)
model.som1 <- som(X=X, grid = somgrid(5,1,"hexagonal"),rlen=5000)

#Determine the best cluster number
#Choose the clusters with the smallest ch_index
#Take five clusters as an example, which represents the smallest ch_index
X2 <- as.data.frame(X2)
som_5_1 <- model.som1$unit.classif
som_5_1 <- data.frame(som_5_1)
som_5_1$number <- 1:435058
environmental_factor<- as.data.frame(X)
environmental_factor$number <- 1:435058
environmental_factor_som <- merge(environmental_factor,som_5_1,by = 'number',all = TRUE)##merge
t <- environmental_factor_som[,c(1:14)]
tt <- as.matrix(t)
quadratic_sum <- sum(tt^2)#Calculate the sum of squares
ssb <- 0 #Sum of Squares Between
ssw <- 0 #Sum of Squares Within
for(i in c(1:5)){
  g <- environmental_factor_som[environmental_factor_som$som_5_1==i,c(2:15)]
  ft <- c(mean(g$TRN_24h),mean(g$TRN_16h),mean(g$TRN_daytime),mean(g$TRN_nighttime),mean(g$TRN_evening),mean(g$Nitrogen_oxides_2010),mean(g$Nitrogen_dioxide_05_10_ave),mean(g$PM25_2010),mean(g$PM10_07_10_ave),mean(g$PM25_10_2010),mean(g$Greenspace_1000m),mean(g$Greenspace_300m),mean(g$water_per_300),mean(g$water_per_1000))
  ft <- as.matrix(ft)
  ft <- ft^2
  ssb <- ssb+length(g[,1])*sum(ft)
  g$TRN_24h <- g$TRN_24h-mean(g$TRN_24h)
  g$TRN_16h <- g$TRN_16h-mean(g$TRN_16h)
  g$TRN_daytime <- g$TRN_daytime-mean(g$TRN_daytime)
  g$TRN_nighttime <- g$TRN_nighttime-mean(g$TRN_nighttime)
  g$TRN_evening <- g$TRN_evening-mean(g$TRN_evening)
  g$Nitrogen_oxides_2010 <- g$Nitrogen_oxides_2010-mean(g$Nitrogen_oxides_2010)
  g$Nitrogen_dioxide_05_10_ave <- g$Nitrogen_dioxide_05_10_ave-mean(g$Nitrogen_dioxide_05_10_ave)
  g$PM25_2010 <- g$PM25_2010-mean(g$PM25_2010)
  g$PM10_07_10_ave <- g$PM10_07_10_ave-mean(g$PM10_07_10_ave)
  g$PM25_10_2010 <- g$PM25_10_2010-mean(g$PM25_10_2010)
  g$Greenspace_1000m <- g$Greenspace_1000m-mean(g$Greenspace_1000m)
  g$Greenspace_300m <- g$Greenspace_300m-mean(g$Greenspace_300m)
  g$water_per_300 <- g$water_per_300-mean(g$water_per_300)
  g$water_per_1000 <- g$water_per_1000-mean(g$water_per_1000)
  g <- as.matrix(g)
  g <- g^2
  ssw <- ssw+sum(g)
}
ssb <- quadratic_sum-ssw
#ch_index <- (ssb / (n - 1)) / (tatol / (435058 - n))
ch_index <- (ssb/4)/(ssw/(435058-5))