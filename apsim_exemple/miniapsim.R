
library(tidyverse)
library(cowplot)
library(data.table)


meteo <- read_csv("meteo.csv") %>% 
  mutate(svpmax = 6.1078*exp(17.269*tmax/(237.3+tmax))*0.1) %>% 
  mutate(svpmin = 6.1078*exp(17.269*tmin/(237.3+tmin))*0.1) %>% 
  mutate(VPDcalc = 0.75*(svpmax-svpmin)*10)

plant <- read_csv("plants.csv")

soils <- read_csv("soils.csv")

    
temp <- plant %>% filter(plant == "maize")
tempsoil <- soils %>% filter(type == "dry")

input <- data.frame(RUE = temp$value[temp$param == "RUE"], 
                    RGR = temp$value[temp$param == "RGR"], 
                    kl1 = temp$value[temp$param == "kl1"], 
                    kl2 = temp$value[temp$param == "kl2"], 
                    kl3 = temp$value[temp$param == "kl3"], 
                    initASW1 = tempsoil$ASW1, 
                    initASW2 = tempsoil$ASW2, 
                    initASW3 = tempsoil$ASW3)
    

    
    # PARAMETERS
    RUE <- input$RUE
    RootGrowthRate <- input$RGR #20
    TEc <- 9
    PotentialDLAI <- 0.1
    k <- 0.45
    InitialLAI <- 1.5
    InitialBiomass <- 45
    VPDfrac <- 1
    
    tinit <- 30
    tmax <- 60
    
    
    meteo <- meteo %>% filter(type == "warm")
    
    soil <- data.frame(param = rep(c("depth", "ll", "ul", "sw", "kl"), 3),
                       layer = rep(c(1:3), each=5), 
                       value = rep(c(300, 50, 100, 100, 0.06), 3)
    )
    
    depth1 <- 300
    depth2 <- 300
    depth3 <- 300
    totDepth <- depth1 + depth2 + depth3
    
    kl1 <- input$kl1
    kl2 <- input$kl2
    kl3 <- input$kl3
    
    sd1 <- 0.5
    sd2 <- 1.5
    sd3 <- 4
    
    # Initialise dataframe
    sim <- data.frame(das = c(tinit)) %>% 
      mutate(root_depth = ifelse(das * RootGrowthRate < totDepth,
                                 das * RootGrowthRate , 
                                 totDepth)) %>% 
      mutate(ASW1 = input$initASW1, 
             ASW2 = input$initASW2, 
             ASW3 = input$initASW3) %>% 
      mutate(totASW = ASW1 + ASW2 + ASW3) %>% 
      
      mutate(ps1 = ifelse(root_depth >= depth1,
                          1 * ASW1 * kl1, 
                          (root_depth / depth1) * ASW1 * kl1)) %>% 
      mutate(ps2 = ifelse(root_depth <= depth1,
                          0, 
                          ifelse(root_depth > depth1+depth2,
                                 1 * ASW2 * kl2, 
                                 ((root_depth-depth1) / depth2) * ASW2 * kl2))) %>% 
      mutate(ps3 = ifelse(root_depth <= depth1+depth2, 
                          0, 
                          ((root_depth - depth1 - depth2)/depth3) * ASW3 * kl3)) %>% 
      mutate(potsupply = ps1 + ps2 + ps3) %>% 
      mutate(lai = InitialLAI) %>% 
      mutate(li = 1 - exp(-k*lai)) %>% 
      mutate(potdemand = meteo$Radn[meteo$DAS == das] * li * RUE / 
               (TEc / (meteo$VPDcalc[meteo$DAS == das]/10))) %>% 
      mutate(sd = potsupply / potdemand) %>% 
      mutate(leafexpeffect = ifelse( sd <= sd1, 
                                     0, 
                                     ifelse(sd > sd2,
                                            1, 
                                            (sd - sd1)/(sd2-sd1)))) %>% 
      mutate(Dlai = leafexpeffect * PotentialDLAI) %>% 
      mutate(transpiration = min(potdemand, potsupply)) %>% 
      mutate(SWaterUse = transpiration) %>% 
      mutate(BioWater = potsupply * TEc / (meteo$VPDcalc[meteo$DAS == das]/10)) %>% 
      mutate(BioLight = meteo$Radn[meteo$DAS == das] * li * RUE) %>% 
      mutate(DBiomass = ifelse(sd > 1, BioLight, BioWater)) %>% 
      mutate(biomass = DBiomass + InitialBiomass)
    
    
    # Update dataframe 
    
    rewater <- 0
    for(i in c((tinit+1):(tmax))){
      
      tempP <- sim[sim$das == i-1,]
      temp <- tempP %>% mutate(das = i)
      
      temp$ASW1 <- tempP$ASW1-(tempP$ps1/tempP$potsupply)*tempP$transpiration
      temp$ASW2 <- tempP$ASW2-(tempP$ps2/tempP$potsupply)*tempP$transpiration
      temp$ASW3 <- tempP$ASW3-(tempP$ps3/tempP$potsupply)*tempP$transpiration
      
      if(i == rewater){
        temp$ASW1 <- soils$ASW1[soils$type == "wet"]
        temp$ASW2 <- soils$ASW2[soils$type == "wet"]
        temp$ASW3 <- soils$ASW3[soils$type == "wet"]
        rewater = rewater + input$rewater
        print(paste0("rewater ",rewater))
      }
      
      temp <- temp %>%
        mutate(root_depth = ifelse(das * RootGrowthRate < totDepth,
                                   das * RootGrowthRate , 
                                   totDepth)) %>% 
        mutate(totASW = ASW1 + ASW2 + ASW3) %>% 
        mutate(ps1 = ifelse(root_depth >= depth1,
                            1 * ASW1 * kl1, 
                            (root_depth / depth1) * ASW1 * kl1)) %>% 
        mutate(ps2 = ifelse(root_depth <= depth1,
                            0, 
                            ifelse(root_depth > depth1+depth2,
                                   1 * ASW2 * kl2, 
                                   ((root_depth-depth1) / depth2) * ASW2 * kl2))) %>% 
        mutate(ps3 = ifelse(root_depth <= depth1+depth2, 
                            0, 
                            ((root_depth - depth1 - depth2)/depth3) * ASW3 * kl3)) %>% 
        mutate(potsupply = ps1 + ps2 + ps3) %>% 
        mutate(lai = tempP$lai + tempP$Dlai) %>% 
        mutate(li = 1 - exp(-k*lai)) %>% 
        mutate(potdemand = meteo$Radn[meteo$DAS == das] * li * RUE / 
                 (TEc / (meteo$VPDcalc[meteo$DAS == das]/10))) %>% 
        mutate(sd = potsupply / potdemand) %>% 
        mutate(leafexpeffect = ifelse( sd <= sd1, 
                                       0, 
                                       ifelse(sd > sd2,
                                              1, 
                                              (sd - sd1)/(sd2-sd1)))) %>% 
        mutate(Dlai = leafexpeffect * PotentialDLAI) %>% 
        mutate(transpiration = min(potdemand, potsupply)) %>% 
        mutate(SWaterUse = transpiration + tempP$SWaterUse) %>% 
        mutate(BioWater = potsupply * TEc / (meteo$VPDcalc[meteo$DAS == das]/10)) %>% 
        mutate(BioLight = meteo$Radn[meteo$DAS == das] * li * RUE) %>% 
        mutate(DBiomass = ifelse(sd > 1, BioLight, BioWater)) %>% 
        mutate(biomass = DBiomass + tempP$biomass)
      
      sim <- rbind(sim, temp)
    }
    
    sim$sim <- rs$counter
    rs$counter <- rs$counter+1
    
    
    rs$sim <- rbind(rs$sim, sim)
    
  })