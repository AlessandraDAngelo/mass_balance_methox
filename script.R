#R script for calculating the mass balance of experimental marine methane concentration as described in Uhlig et al. (2017).
#D'Angelo A. and Loose B. (2020) - University of Rhode Island, Graduate School of Oceanography

#load dataframe "data"

library(dplyr) 
library(ggplot2)

#create a column with values for dates
start <- as.Date(sort(data$Date)[1],format="%Y-%m-%d")
end   <- as.Date(sort(data$Date, TRUE)[1],format="%Y-%m-%d")
data$valueDate <- ifelse(data$Date != start -1,(data$Date - start),NA)
#to set values as 2-digit numbers
data$valueDate <- sprintf("%02d", as.numeric(as.character(data$valueDate)))
data$Date <- as.Date(data$Date, tryFormats = "%Y-%m-%d")
#order the df by sample
data <- data[order(data$bag.Num),]

#calculate the mass balance 
#subset df, based on samples
DF = data.frame()
Ubags = unique(data$bag.Num)

for (B in Ubags){
  
  BagID = data$bag.Num == B
  
  BagSub = data[BagID,]
  
  # Volumes 
  a=1027 #g/L
  Vw=(DF$bag.M/a) #L
  Vhs=0.1 #L
  
  # Ideal gas constant (mol/L)
  DF$IGC=(101.325/(8.314*272))   
  print(DF$IGC)
  
  # Pressure correction
  P_SSIM=BagSub$Pssim/(BagSub$airP/1.332) #ppm
  
  # Dilution factor correction
  BagSub$vol.L<- BagSub$vol.mL/1000 #L
  SSIM.dil.factor=(0.022/BagSub$vol.L) 
  
  # Calculate methane moles in headspace
  #pCH4 is the partial pressure of methane recorded from the Picarro analyzer
  pCH4.hs=SSIM.dil.factor*P_SSIM*BagSub$HP.12CH4.Mean #ppm
  chi.hs=pCH4.hs*10^-6 #ppp
  
  
  #Bunsen solubility coefficient: https://github.com/URIGSO/Gas-Solubility-Codes/tree/R_code
  
  #calculate the methane gas solubility according to the Bunsen solubility coefficient (Yamamoto et al., 1976)  bunsen <- function(T, S){
    
    A1 <- -67.1962
    A2 <- 99.1624
    A3 <- 27.9015
    B1 <- -0.072909
    B2 <- 0.041674
    B3 <- -0.0064603
    T<- 1
    T.K <- T + 273.15
    
    bunsen.coeff.L.L.x <- exp(A1 + A2*(100/T.K) + A3*(log(T.K/100)) + S*(B1 + B2*(T.K/100) + B3*((T.K/100)^2)))
    
    bunsen.coeff.L.L <- c(bunsen.coeff.L.L.x)
    
    return(bunsen.coeff.L.L)
}
#Calculate Bunsen solubility coefficients using water bath T and in-situ S
  beta_sol=bunsen(1,BagSub$Sal) 
  print(beta_sol)
  #beta solubility in moles
  beta_sol.moles=beta_sol*BagSub$IGC #mol/L
  print(beta_sol.moles)
  
  # Volumes: convert the volume of the standard injected from ml to L
  BagSub$Vsyr.spike.L<- BagSub$Vol.spike_ML/1000 #L
  BagSub$Vhs = Vhs+BagSub$Vsyr.spike.L[1] #Vhs initial = vol injected (100ml)+ Vspyke #L
  print(BagSub$Vhs)
  
  #create a column with the entries of easurements of bags
  BagSub$run_num<- 1:nrow(BagSub)
  #give an index for the run_num
  Run<- unique(BagSub$run_num)
  #set Run as a list
  r<- list(Run)
  
  #ifelse statement to account for Volume leftover at t=0 and Volume leftover after t=0
  BagSub$V.leftover<- NA
  for (i in r) {
    
    BagSub$V.leftover[i]<- BagSub$Vhs[i]- cumsum(BagSub$vol.L[i]) #L
  }
  BagSub$V.leftover[1] <- BagSub$Vhs[1] #L
  print(BagSub$V.leftover)#L
  
  for (i in r) {
    # Calculate CH4 moles in water
    BagSub$moles_w[i]=Vw[i]*beta_sol.moles[i]*chi.hs[i] #mol
    print(BagSub$moles_w[i])
    
    # Calculate CH4 moles in headspace
    BagSub$moles_hs= BagSub$IGC*chi.hs[i]*BagSub$V.leftover #mol
    print(BagSub$moles_hs)
    
    # Calculate CH4 moles sampled
    BagSub$sampled_moles<- BagSub$vol.L*chi.hs*BagSub$IGC
    print(BagSub$sampled_moles)
    
    # Calculate the cumulative sum of moles row by row
    cumsum.sampled_moles<- cumsum(BagSub$sampled_moles)
    print(cumsum.sampled_moles)
    
  }
  
  # Calculate the cumulative sum by days
  # Create a column with the entries of measurements by values of dates
  run_num_day<- BagSub%>%
    group_by(valueDate)%>%
    dplyr::summarise(run_num_day=  max(run_num))%>%
    as.data.frame()
  
  # Merge this with the dataframe BagSub
  BagSub<- left_join(BagSub, run_num_day)
  
  #create a column in run_num_day df with the cumulative sum of the sampled moles relative to the last run number of the day 
  #create an empty column
  BagSub$cumsum.sampled_moles_d<- NA
  
  for (i in r) {
    BagSub$cumsum.sampled_moles_d[i]<- cumsum.sampled_moles[BagSub$run_num_day[i]]
    print(BagSub$cumsum.sampled_moles_d)
  }
  
  #create an empty column for total moles
  BagSub$tot_mass = NA
  
  for (j in BagSub$run_num) {
    
    if (j==1){
      BagSub$tot_mass[1] <- BagSub$moles_hs[1]+BagSub$moles_w[1]
      print(BagSub$tot_mass)
    }
    else{
      BagSub$tot_mass[j] <- BagSub$moles_hs[j]+BagSub$moles_w[j]+BagSub$cumsum.sampled_moles_d[j-1]
    }
    print(BagSub$tot_mass)
    
  }
  
DF <- rbind(DF,BagSub)

#References:
#Yamamoto, S., J. B. Alcauskas, and T. E. Crozier. 1976. Solubility of methane in distilled water and seawater. J. Chem. Eng. Data 21: 78–80. doi:10.1021/je60068a029.
#Uhlig, C., Loose, B. 2017. Using stable isotopes and gas concentrations for independent constraints on microbial methane oxidation at Arctic Ocean temperatures. Limnol. Oceanogr. Methods, 15, 737-751. https://doi.org/10.1002/lom3.10199.
