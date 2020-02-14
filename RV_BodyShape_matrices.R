##############################
## RV coefficient of body shape matrices
##
## Matt Brachmann (PhDMattyB)
##
## 2020-02-13
##
##############################
dir.create('~/PhD/SNP Demographic modelling/RV_coeff')

setwd('~/PhD/SNP Demographic modelling/RV_coeff')

library(rsed)
library(tidyverse)
library(FactoMineR)

le_data = read_csv('GSTVMF_Morph_Eco_Geno.csv') %>% 
  filter(Lake %in% c('Galtabol', 
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn')) %>% 
  select(Specimen.ID, 
         Lake, 
         Morph, 
         BP, 
         Vector, 
         LaMorph, 
         carbon, 
         nitrogen, 
         q1)
fixed_id = sed_substitute(le_data$Specimen.ID, 
                          '-',
                          '') %>% 
  as_tibble()
le_data = bind_cols(le_data, 
                    fixed_id) %>% 
  na.omit()

Body_shape_noallo = read_csv('SIASamples_PWS_StandSexAllo_GSTV.csv') 
Body_shape_noallo = read_csv('AllLakes_AllometryExcluded_PWS_Combined.csv') %>% 
  filter(Lake %in% c('Galtabol', 
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn'))

id_fix = sed_substitute(Body_shape_noallo$Specimen.ID, 
                        '-',
                        '') %>% 
  as_tibble()
Body_shape_noallo = bind_cols(Body_shape_noallo, 
                              id_fix) %>% 
  select(value, 
         contains('PW'), 
         UNIX, 
         UNIY, 
         CS)


data = inner_join(le_data, 
                  Body_shape_noallo, 
                  by = 'value') %>% 
  group_by(Vector) #%>% 
  # filter(Vector == 'GSBPI')
  # nest()


# lapply alternative ------------------------------------------------------
## carbon and nitrogen regression

data_morphpair = split(data, data$Vector)

data_lapply = lapply(data_morphpair, function(x){
  lm(cbind(PW1X, PW1Y, PW2X, PW2Y,PW3X, PW3Y, PW4X, PW4Y, PW5X, PW5Y, PW6X, 
           PW6Y,PW7X, PW7Y, PW8X, PW8Y,PW9X, PW9Y, PW10X, PW10Y,PW11X, PW11Y, 
           PW12X, PW12Y,PW13X, PW13Y, PW14X, PW14Y,
           PW15X, PW15Y, PW16X, PW16Y,
           PW17X, PW17Y, PW18X, PW18Y, 
           PW19X, PW19Y, UNIX, UNIY) ~ carbon + nitrogen, data = x)
})

fitted = data_lapply %>% 
  lapply(function(x){
    fitted.values(x)
  })


## q1 regression on body shape

q1_lapply = lapply(data_morphpair, function(x){
  lm(cbind(PW1X, PW1Y, PW2X, PW2Y,PW3X, PW3Y, PW4X, PW4Y, PW5X, PW5Y, PW6X, 
           PW6Y,PW7X, PW7Y, PW8X, PW8Y,PW9X, PW9Y, PW10X, PW10Y,PW11X, PW11Y, 
           PW12X, PW12Y,PW13X, PW13Y, PW14X, PW14Y,
           PW15X, PW15Y, PW16X, PW16Y,
           PW17X, PW17Y, PW18X, PW18Y, 
           PW19X, PW19Y, UNIX, UNIY) ~ q1, data = x)
})

fitted_q1 = q1_lapply %>% 
  lapply(function(x){
    fitted.values(x)
  })


# RV coefficients of fitted matrices --------------------------------------

coeffRV(fitted$GSBPI, fitted_q1$GSBPI)
coeffRV(fitted$SLGBPI, fitted_q1$SLGBPI)
coeffRV(fitted$SLGBPL, fitted_q1$SLGBPL)
coeffRV(fitted$TLGBPL, fitted_q1$TLGBPL)
coeffRV(fitted$TSBPL, fitted_q1$TSBPL)
coeffRV(fitted$VBRSIL, fitted_q1$VBRSIL)
