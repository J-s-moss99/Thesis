
hist(Thorn$LngtClass[Thorn$LngtCode=="."])  # these are 0.1 mm
hist(Thorn$LngtClass[Thorn$LngtCode=="0"])  # these are mm
hist(Thorn$LngtClass[Thorn$LngtCode=="5"])  # these are cm
hist(Thorn$LngtClass[Thorn$LngtCode=="1"])  # these are cm
hist(Thorn$LngtClass[is.na(Thorn$LngtCode)])  # these do not exist

Thorn <- Thorn[!is.na(Thorn$LngtCode),]
Thorn$LngtCode <- as.character(Thorn$LngtCode)

# convert everything to cm 
Thorn$LngtClass[Thorn$LngtCode=="."] <- Thorn$LngtClass[Thorn$LngtCode=="."]/100   # I was wrong in my mail, mm to cm is not *10 but /10 ... so, now we have to divide by 100
Thorn$LngtClass[Thorn$LngtCode=="0"] <- Thorn$LngtClass[Thorn$LngtCode=="0"]/10
Thorn$LngtClass[Thorn$LngtCode=="5"] <- Thorn$LngtClass[Thorn$LngtCode=="5"]/1
Thorn$LngtClass[Thorn$LngtCode=="1"] <- Thorn$LngtClass[Thorn$LngtCode=="1"]/1

# check conversion
hist(Thorn$LngtClass[Thorn$LngtCode=="."])  # these are 0.1 mm
hist(Thorn$LngtClass[Thorn$LngtCode=="0"])  # these are mm
hist(Thorn$LngtClass[Thorn$LngtCode=="5"])  # these are cm
hist(Thorn$LngtClass[Thorn$LngtCode=="1"])  # these are cm
# looks good!

head(Thorn)

# in some cases, a subsample of was taken
table(Thorn$SubFactor)
# multiply the numbers with the subfactor (and round)
Thorn$HLNoAtLngt <- round(Thorn$HLNoAtLngt*Thorn$SubFactor)
hist(Thorn$HLNoAtLngt,breaks=50)
hist(Thorn$TotalNo,breaks=50)   # not at length


# MERGE HL WITH HH DATA
# stations (HH data) without RJC landings are not of interest, so, no problem if we lose them
# in case you would like to calculate an abundance index, 
# you should merge the other way around (to not lose you zero observations)

library(dplyr)

head(HH_BTS)
Thorn_w_date <- left_join(Thorn[,c("Survey", "Quarter", "Country", "Ship", "Gear", "StNo",  "HaulNo",  "Year", "LngtClass", "HLNoAtLngt")], 
                         HH_BTS[,c("Survey", "Quarter", "Country", "Ship", "Gear", "StNo", "HaulNo", "Year", "Month", "Day")])



# create an empty list
Elefan_rjc <- list()

# loop over all the data
for( i in 1:nrow(Thorn_w_date)){
  # create a dataframe per observation
  df.record <- data.frame(
  Length = rep(Thorn_w_date$LngtClass[i],Thorn_w_date$HLNoAtLngt[i]),
  Year =  rep(Thorn_w_date$Year[i],Thorn_w_date$HLNoAtLngt[i]),
  Month = rep(Thorn_w_date$Month[i],Thorn_w_date$HLNoAtLngt[i]),
  Day = rep(Thorn_w_date$Day[i],Thorn_w_date$HLNoAtLngt[i]))
  
  # add it in the list
  Elefan_rjc[[i]] <- df.record
}

# put all the list items together in a dataframe
Elefan_df <- do.call(rbind,Elefan_rjc)
head(Elefan_df)

Elefan_df$Date <- as.Date(with(Elefan_df, paste(Year, Month, Day,sep="-")), "%Y-%m-%d")
head(Elefan_df$Date)

library(TropFishR)

# create lfq data
lfq_dat <- lfqCreate(Elefan_df,Lname = "Length", Dname = "Date", aggregate_dates = TRUE,
                     length_unit = "cm", bin_size = 2, plot=TRUE, plus_group=c(TRUE,75))



