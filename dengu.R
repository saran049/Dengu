library(forecast) 
library(fpp) 
library(caret) 
library(neuralnet) 
library(randomForest) 
library(psych)
library(VIM) 
library(mice) 
library(ResourceSelection) 
library(corrplot) 
library(party)
library(magrittr)
library('dplyr')


par(mfrow=c(1,1))

#dengue_features_train_combined is one week shifted data:
            
train.data=read.csv("dengue_features_train_combined.csv")
test.data=read.csv("dengue_features_test.csv")

#check the missing data:
aggr(train.data, prop=F,numbers=F,plot = F)
aggr(test.data, prop=F,numbers=F,plot = F)


# Every column in the dataset
features = c('ndvi_ne','ndvi_nw','ndvi_se','ndvi_sw','precipitation_amt_mm',
             'reanalysis_air_temp_k','reanalysis_avg_temp_k','reanalysis_dew_point_temp_k',
             'reanalysis_max_air_temp_k','reanalysis_min_air_temp_k',
             'reanalysis_precip_amt_kg_per_m2','reanalysis_relative_humidity_percent',
             'reanalysis_sat_precip_amt_mm','reanalysis_specific_humidity_g_per_kg',
             'reanalysis_tdtr_k','station_avg_temp_c','station_diur_temp_rng_c',
             'station_max_temp_c','station_min_temp_c','station_precip_mm')


#Features that are highly correlated to Total cases:
f_related =c('reanalysis_min_air_temp_k',
             'reanalysis_specific_humidity_g_per_kg',
             'reanalysis_dew_point_temp_k',
             'station_avg_temp_c',
             'station_min_temp_c')


#adding missing values from the former values:
train.data[features] %<>% na.locf(fromLast = TRUE)
test.data[features] %<>% na.locf(fromLast = TRUE)


#dividing the two cities in train data:
a = nrow(train.data)
sj.train = data.frame()
iq.train = data.frame()
for(i in 1:a) {
  if (train.data$city[i] == "sj"){
    sj.train = rbind(sj.train,train.data[i,])
  }
  else {
    iq.train = rbind(iq.train,train.data[i,])
  }
}

#dividing the two cities in test data:
b = nrow(test.data)
sj.test = test.data[1:260,]
iq.test = test.data[261:416,]


#correlation between the column:
iq.cor=cor(iq.train[,5:25]) 
corrplot(iq.cor, method="circle", type = "lower", bg = "black")

sj.cor = cor(sj.train[,5:25])
corrplot(sj.cor, method = "circle", type = "lower",bg = "black")

cor(sj.train[,5:25])

mycor=cor(train.data[,5:25]) 
print(mycor)
highlyCorrelated <- findCorrelation(mycor, cutoff=0.1)

#formulating the time series:
sj.series=ts(sj.train$total_cases, frequency=52, start=c(1990,30,04)) 
iq.series=ts(iq.train$total_cases, frequency=52, start=c(2000,07,01)) 


# checking for trend and seasonality:
par(mfrow=c(1,1))
plot(sj.series,main = "City SJ",lwd = 2)
plot(iq.series,main = "City IQ")


plot(stl(sj.series,s.window = "periodic"))
acf(sj.series)
pacf(sj.series)

plot(stl(iq.series,s.window = "periodic"))
Acf(iq.series)
pacf(iq.series)


#removing trend and seasonality:
plot(stl(sj.series,s.window = "periodic"))
ndiffs(sj.series)
sj.diff1 <- diff(sj.series,lag = 1)
plot(sj.diff1,main = "City SQ")

Acf(sj.diff1)

plot(stl(iq.series,s.window = "periodic"))
ndiffs(iq.series)
iq.diff1 <- diff(iq.series,lag = 1)
plot(iq.diff1,main = "City IQ")
Acf(sj.diff1)


### auto arima.............

auto.sj <- auto.arima(sj.series, xreg =sj[f_related])
summary(auto.sj)

auto.iq <- auto.arima(iq.series, xreg = iq[f_related])
summary(auto.iq)



###.................... xreg seasonal arima...................


xreg.sj <- Arima(sj.series, order=c(1,1,1), xreg =sj[f_related] ) 
summary(xreg.sj)
#AIC = 7515.17, RMSE = 13.32
accuracy(xreg.sj)

xreg.sj.forecast <- (forecast(xreg.sj, h = 260, xreg = sj.test[f_related]))
tsdisplay(residuals(xreg.sj))
plot(xreg.sj.forecast,xlab = "Years",ylab = "range",
     main = "Forecast value for city SJ")
Box.test(residuals(xreg.sj.forecast),fitdf=7, lag=151, type="Ljung")  
xreg.sj.forecast$mean

xreg.iq <- Arima(iq.series, order=c(1,1,1), xreg =iq[f_related] ) 
summary(xreg.iq)
accuracy(xreg.iq)
#AIC = 3534.11, RMSE = 7.134
xreg.iq.forecast <- (forecast(xreg.iq, h = 156, xreg = iq.test[f_related]))
tsdisplay(residuals(xreg.iq))
plot(xreg.iq.forecast,xlab = "Years",ylab = "range",
     main = "Forecast value for city IQ")
Box.test(residuals(xreg.iq.forecast),fitdf=7, lag=151, type="Ljung")  

head(xreg.iq.forecast$mean)
head(round(xreg.iq.forecast$mean))
sj.final <- data.frame(sj.test[,1:3],total_cases = round(xreg.sj.forecast$mean))
iq.final <- data.frame(iq.test[,1:3],total_cases = round(xreg.iq.forecast$mean))  
final_data <- data.frame()
final_data <- bind_rows(sj.final,iq.final)

plot(sj.final$total_cases,sj.final$weekofyear)

aaa <- table(sj.final$weekofyear,sj.final$weekofyear)
barplot(aaa)
  
write.csv(final_data, file = 'xreg_arima_predicted_Solution.csv', row.names = F)
nrow(final_data)

#.................Seasonal Arima......................................


seasonal.sj <- Arima(sj.series, order=c(1,1,1), seasonal=c(0,0,0)) 
summary(seasonal.sj)
#AIC = 7517.73 , RMSE = 13.43
sj.forecast.seasonal=forecast(seasonal.sj, h = 260)
plot(sj.forecast.seasonal)
tsdisplay(residuals(sj.forecast.seasonal))
Box.test(residuals(sj.forecast.seasonal),fitdf=2, lag=151, type="Ljung")  

seasonal.iq <- Arima(iq.series, order=c(1,1,1), seasonal=c(0,0,0)) 
summary(seasonal.iq)
#AIC = 3525.37 , RMSE = 7.15
iq.forecast.seasonal=forecast(seasonal.iq, h = 156)
plot(iq.forecast.seasonal)
tsdisplay(residuals())
Box.test(residuals(iq.forecast.seasonal),fitdf=2, lag=151, type="Ljung")  

sj.final1 <- data.frame(sj.iq.forecast.seasonaltest[,1:3],total_cases = round(sj.forecast.seasonal$mean))
iq.final1 <- data.frame(iq.test[,1:3],total_cases = round(iq.forecast.seasonal$mean))  
final_data1 <- data.frame()
final_data1 <- bind_rows(sj.final1,iq.final1)
nrow(final_data1)
write.csv(final_data1, file = 'seasonal_arima_predicted_Solution.csv', row.names = F)

#............Random Forest model 1...................................

sj.rf.model <- randomForest(  total_cases~
                              ndvi_ne + ndvi_nw + ndvi_se + ndvi_sw + precipitation_amt_mm +	
                              reanalysis_air_temp_k + reanalysis_avg_temp_k + 
                              reanalysis_dew_point_temp_k + reanalysis_max_air_temp_k + 
                              reanalysis_min_air_temp_k + reanalysis_precip_amt_kg_per_m2 +
                              reanalysis_relative_humidity_percent + 
                              reanalysis_sat_precip_amt_mm + 
                              reanalysis_specific_humidity_g_per_kg + reanalysis_tdtr_k + 
                              station_avg_temp_c +
                              station_diur_temp_rng_c + station_max_temp_c + 
                              station_min_temp_c + station_precip_mm,
                              data = sj.train)

print(sj.rf.model)
plot(sj_rf_model)
summary(sj.rf.model)
sj.rf.prediction<- predict(object=sj.rf.model, sj.test)
head(sj.rf.prediction)

iq.rf.model <- randomForest(total_cases ~
                              ndvi_ne + ndvi_nw + ndvi_se + ndvi_sw + precipitation_amt_mm +	
                              reanalysis_air_temp_k + reanalysis_avg_temp_k + 
                              reanalysis_dew_point_temp_k + reanalysis_max_air_temp_k + 
                              reanalysis_min_air_temp_k + reanalysis_precip_amt_kg_per_m2 +
                              reanalysis_relative_humidity_percent + 
                              reanalysis_sat_precip_amt_mm +
                              reanalysis_specific_humidity_g_per_kg + reanalysis_tdtr_k +
                              station_avg_temp_c +
                              station_diur_temp_rng_c + station_max_temp_c +
                              station_min_temp_c + station_precip_mm
                            , data = iq.train)

print(iq.rf.model)
plot(iq.rf.model)
summary(iq.rf.model)
iq.rf.prediction<- predict(object=iq.rf.model, iq.test)

sj.final2 <- data.frame(sj.test[,1:3],total_cases = round(sj.rf.prediction))
iq.final2 <- data.frame(iq.test[,1:3],total_cases = round(iq.rf.prediction))  
final_data2 <- data.frame()
final_data2 <- bind_rows(sj.final2,iq.final2)


write.csv(final_data2, file = 'random_forest_model1.csv', row.names = F)


#..........Random Forest model 2....................#

sj.rf1.model <- randomForest(total_cases ~
                               reanalysis_specific_humidity_g_per_kg +
                               reanalysis_dew_point_temp_k + 
                               station_avg_temp_c +
                               station_min_temp_c
                             , data = sj.train)

print(sj.rf1.model)


sj.rf1.prediction<- predict(object=sj.rf1.model, sj.test)

iq.rf1.model <- randomForest(total_cases ~
                               reanalysis_specific_humidity_g_per_kg +
                               reanalysis_dew_point_temp_k + 
                               station_avg_temp_c +
                               station_min_temp_c
                             , data = iq.train)

print(iq.rf1.model)
plot(iq.rf1.model)

iq.rf1.prediction<- predict(object=iq.rf1.model, iq.test)


sj.final3 <- data.frame(sj.test[,1:3],total_cases = round(sj.rf1.prediction))
iq.final3 <- data.frame(iq.test[,1:3],total_cases = round(iq.rf1.prediction))  
final_data3 <- data.frame()
final_data3 <- bind_rows(sj.final3,iq.final3)


write.csv(final_data2, file = 'random_forest_model2.csv', row.names = F)

class(final_data)
str(final_data)


