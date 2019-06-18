#Overall:   1) Breaks up image data on a per plate basis
#           2) MAKE SURE YOU HAVE SEPARATED YOUR GITTER OUTPUT FROM COL SIZE
#           3) Rearranges colony size data and combines with plate mev files
#           4) Creates new 'clean' dataframe that's easier to read


#Splitting MEV files into single plate files
#load file names
Mev_files <- dir("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/1_Data_FLEX_9by9_bothsignals/")
PathDirectory <- as.character("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/1_Data_FLEX_9by9_bothsignals/")
# set directory in which the colony size data is saved
setwd("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/1_Individual_plate_data/")

lapply(1:length(Mev_files), function(i){
  #1 create a character string containing the path to each image in the directory
  s = paste(as.character(PathDirectory), as.character(Mev_files[i]), sep="")
  #2 read each mev file
  j <- read.table(s, header = TRUE)
  
  #3 split into different grids
  w <- split (j, j$MR) [['1']]
  #x is 1:1, y is 1:2, z is 1:3that's DA 1 2 and 3 for 1-9 or 10,11,12 for 10-18
  x <- split (w, w$MC) [['1']]
  y <- split (w, w$MC) [['2']]
  z <- split (w, w$MC) [['3']]
  
  m <- split (j, j$MR) [['2']]
  #n is 2:1, o is 2:2, p is 2:3 that's DA 4 5 and 6 for 1-9 or 13,14,15 for 10-18
  n <- split (m, m$MC) [['1']]
  o <- split (m, m$MC) [['2']]
  p <- split (m, m$MC) [['3']]
  
  q <- split (j, j$MR) [['3']]
  #r is 3:1, s is 3:2, t is 3:3 that's DA 7  8 and 9 for 1-9 or 16,17,18 for 10-18
  r <- split (q, q$MC) [['1']]
  s <- split (q, q$MC) [['2']]
  t <- split (q, q$MC) [['3']]
  
  
  #now create new tables for each of these files
  #will keep raw data for each plate here
  write.table(x, paste(as.character(Mev_files[i]),"1.txt", sep = ""), sep = "\t", row.names=F)
  write.table(y, paste(as.character(Mev_files[i]),"2.txt", sep = ""), sep = "\t", row.names=F)
  write.table(z, paste(as.character(Mev_files[i]),"3.txt", sep = ""), sep = "\t", row.names=F)
  write.table(n, paste(as.character(Mev_files[i]),"4.txt", sep = ""), sep = "\t", row.names=F)
  write.table(o, paste(as.character(Mev_files[i]),"5.txt", sep = ""), sep = "\t", row.names=F)
  write.table(p, paste(as.character(Mev_files[i]),"6.txt", sep = ""), sep = "\t", row.names=F)
  write.table(r, paste(as.character(Mev_files[i]),"7.txt", sep = ""), sep = "\t", row.names=F)
  write.table(s, paste(as.character(Mev_files[i]),"8.txt", sep = ""), sep = "\t", row.names=F)
  write.table(t, paste(as.character(Mev_files[i]),"9.txt", sep = ""), sep = "\t", row.names=F)
  
  
})



#############################################
##########################

#Add strain name, and FLEX number based on 1536 array, add colony size from gitter package
#Add FLEX tables that have already been split into plates

#load file names
Single_files <- dir("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/1_Individual_plate_data/")
PathDirectory_single <- as.character("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/1_Individual_plate_data/")
setwd("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/3_Colony_size_and_Individual_plate_data/")


lapply(1:length(Single_files), function(a){
  #1/ create a character string containing the path to each image in the directory
  b = paste(as.character(PathDirectory_single), as.character(Single_files[a]), sep="")
  c = read.table(b, header = TRUE)
  

  Plate_Num_name <- (sub(".txt", "", Single_files[a]))
  Plate_Num_vec <- strsplit(Plate_Num_name, split = "_") [[1]]
  Plate_Num <- Plate_Num_vec[3]
  Strain <- Plate_Num_vec[1]
  Strain_ID <- rep(Strain, length(c$UID))
  FLEX_ID <- rep(Plate_Num, length(c$UID))
  
  e <- cbind(c,Strain_ID, FLEX_ID)
  
  #Invert all of the columns numbers for each row and sort by row number, followed by column number    
  invcolumn <- as.numeric(48:1)
  e <- cbind(data.frame(invcolumn,e))
  e <- e[order(e$invcolumn),]
  e <- e[order(e$R),]
  
  #INV COLUMN IS ACTUALLY THE CORRECT COLUMN NUMBER FOR THE FLEX ARRAY
  
  Col_directory <- dir("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/2_Colony_Size_output_gitter/")
  PathDirectory_col <- as.character("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/2_Colony_Size_output_gitter/")
  
  m = paste(as.character(PathDirectory_col), as.character(Col_directory[a]), sep="")
  o = read.table(m, sep="\t", row.names=NULL, skip=3)
  
#JUST KIDDING DONT DO THIS HERE
  #Invert all of the columns numbers for each row and sort by row number, followed by column number    
  #invcolumn <- as.numeric(48:1)
  #o <- cbind(data.frame(invcolumn,o))
  #o <- o[order(o$invcolumn),]
  #o <- o[order(o$V1),]
  
  
  names(o) <- c("row", "col", "size", "circ", "flags")
  
  #bind size table to mev data  
  y <- cbind(data.frame(o,e))
  
  write.table(y, paste("FLEXCol_Size_", as.character(Single_files[a]), sep = ""), sep = "\t", row.names = F)  
  
  
})


##################################
##################################
##################################

#Make new dataframe with cleaner data

#Make an array of 1-384 with every number doubled and then every 1-24 doubled 
#to line up with rows and columns that should be grouped
v <- (1:24)
v <- rep(v, each=2)
v1 <- (25:48)
v1 <- rep(v1, each=2)
v2 <- (49:72)
v2 <- rep(v2, each=2)
v3 <- (73:96)
v3 <- rep(v3, each=2)
v4 <- (97:120)
v4 <- rep(v4, each=2)
v5 <- (121:144)
v5 <- rep(v5, each=2)
v6 <- (145:168)
v6 <- rep(v6, each=2)
v7 <- (169:192)
v7 <- rep(v7, each=2)
v8 <- (193:216)
v8 <- rep(v8, each=2)
v9 <- (217:240)
v9 <- rep(v9, each=2)
v10 <- (241:264)
v10 <- rep(v10, each=2)
v11 <- (265:288)
v11 <- rep(v11, each=2)
v12 <- (289:312)
v12 <- rep(v12, each=2)
v13 <- (313:336)
v13 <- rep(v13, each=2)
v14 <- (337:360)
v14 <- rep(v14, each=2)
v15 <- (361:384)
v15 <- rep(v15, each=2)

# FLEX_quad is every 24 numbers (1 row in 384 plate), repeated 2 times (1 row in 1536 plate), 
#   then that row repeated again because quadruplets
FLEX_quad <- c(v, v, v1, v1, v2, v2, v3, v3, v4, v4, v5, v5, v6, v6, v7, v7, v8, v8, v9, v9, v10, v10, v11, v11, v12, v12, v13, v13, v14, v14, v15, v15)
Quad_vector <- (1:4) 
Quad_vec <- c(Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector, Quad_vector)
Quad_num <- c(Quad_vec, Quad_vec, Quad_vec, Quad_vec, Quad_vec, Quad_vec)
# Quad_number is  (1, 2, 3, 4) 384x
Quad_number <- c(Quad_num, Quad_num, Quad_num, Quad_num)


#Reshape data to combine quadruplets and extract columns that you want   

Col_files <- dir("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/3_Colony_size_and_Individual_plate_data/")
PathDirectory_quad <- as.character("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/3_Colony_size_and_Individual_plate_data/")
setwd("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/4_Reshaped_filtered_normalized_ind_plates_CORRECT_ALIGNMENT/")


FLEX_files <- dir("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/FLEX_split/")
PathDirectory_FLEX <- as.character("~/Documents/KD2018_Screen_data/20180227_42hours_FLEX_Cleancode/FLEX_split/")


lapply(1:length(Col_files), function(n){
  m = paste(as.character(PathDirectory_quad), as.character(Col_files[n]), sep="")
  o = read.table(m, header = TRUE)
  #add on which quadruplet number the colony is a part of  
  #FLEX_quad should correspond to the entire quadruplet and we will reshape based on this
  p <- cbind(data.frame(FLEX_quad, o))
  #order the quadruplets s.t. they're all together
  q = p[order(p$FLEX_quad),]
  #give each colony an arbitrary number to ID it from the rest of the quadruplet in the row
  q = cbind(data.frame(Quad_number, q))
  q = q[order(q$Quad_number),]
  q = q[order(q$FLEX_quad),]
  
  FLEXname = paste(as.character(PathDirectory_FLEX), as.character(FLEX_files[n]), sep="") 
  FLEX = read.csv(FLEXname, header = TRUE, sep="\t") 
  
  q = cbind(data.frame(q, FLEX))
  #Clean up data and put it in a format that's easy to read. 
  #Can always go back to original "individual plate" files.
  
  GFP_Mean <- as.numeric(q$MNA)
  RFP_Mean <- as.numeric(q$MNB)
  GFP_Med <- as.numeric(q$MedA)
  RFP_Med <- as.numeric(q$MedB)
  size_Ty <- as.numeric(q$SA)
  size_Scan <- as.numeric(q$size)
  circ_Scan <- as.numeric(q$circ)
  Strain <- as.character(q$Strain_ID)
  Plate <- as.character(q$FLEX_ID)
  Row <- as.numeric(q$R)
  Col <- as.numeric(q$invcolumn)
  #FLEX_quad is going to give the quadruplet number, but on a per plate basis
  #So, you should have 1-384 for every plate
  FLEX_quad <- as.numeric(q$FLEX_quad)
  Quad_number <- as.numeric(q$Quad_number)
  
 
  dat_final <- cbind(data.frame(Strain, Plate, Row, Col, q$X384_Plates, q$X384_Rows, q$X384_Column, q$SGD, q$Symbol, q$comment, q$Description, q$Q, FLEX_quad, Quad_number, size_Ty, size_Scan, circ_Scan, GFP_Mean, RFP_Mean, GFP_Med, RFP_Med))

########FILTER##############
  Array_filter1_border <- subset(dat_final, q.SGD != "empty FLEX vector") 
  Array_filter2_blank <- subset(Array_filter1_border, q.SGD != "blank") 
  Array_filter3_blank <- subset(Array_filter2_blank, q.SGD != "probably YOR202W") 
  
  #Replace all zero values to NA
  Array_filter3_blank[Array_filter3_blank == 0] <- NA
  
  ####  intensityFBorder <- subset(intensityFBorder, ORF != "probably YOR202W")
  
  # Filter Size
  SizeFilterSmall <- 1500 
  SizeFilterLarge <- 6000
  Array_filter3_blank$GFP_Mean[Array_filter3_blank$size_Scan < SizeFilterSmall] <- NA
  Array_filter3_blank$GFP_Mean[Array_filter3_blank$size_Scan > SizeFilterLarge] <- NA
  Array_filter3_blank$GFP_Med[Array_filter3_blank$size_Scan < SizeFilterSmall] <- NA
  Array_filter3_blank$GFP_Med[Array_filter3_blank$size_Scan > SizeFilterLarge] <- NA
  
  
  Filtd_Array = Array_filter3_blank
  
  
  ############# Do LOESS analysis now ############# 
  Filtd_Array <- within(Filtd_Array, { 
    meanlogRatio <- log2(GFP_Mean/RFP_Mean)
    meanlogBrightness <- log2(GFP_Mean*RFP_Mean)
    medlogRatio <- log2(GFP_Med/RFP_Med)
    medlogBrightness <- log2(GFP_Med*RFP_Med)  
    meanRatio <- (GFP_Mean/RFP_Mean)
    meanBrightness <- (GFP_Mean*RFP_Mean)
    medRatio <- (GFP_Med/RFP_Med)
    medBrightness <- (GFP_Med*RFP_Med)  
  })
  
  # Define LOESS span
  # When you set the span to 1, you make the distribution normal. 
  # It really gets rid of the background colonies that gave extremely low brightness
  # but you might inherently be getting rid of your low GFP high RFP things, so maybe 0.75
  span = 1
  
  # LOESS normalization
  
  # - apply normalization
  Filtd_Array <- within(Filtd_Array, {
    meanlogLOESS <- predict(loess(meanlogRatio ~ meanlogBrightness, span = span, na.action = na.exclude))
    medlogLOESS <- predict(loess(medlogRatio ~ medlogBrightness, span = span, na.action = na.exclude))
    meanlogRatioNorm <- meanlogRatio - meanlogLOESS
    medlogRatioNorm <- medlogRatio - medlogLOESS
  })
  

###########Now we have all of the data for normalized mean values, need to reshape data to add Z score for each individual plate###########

  
  y <- reshape(Filtd_Array, idvar= c("FLEX_quad", "Plate"), timevar = "Quad_number", direction="wide")
  
  
  
   New_form <- cbind(data.frame(y$Plate, y$FLEX_quad, y$Strain.1, y$q.X384_Plates.1, y$q.X384_Rows.1, y$q.X384_Column.1, y$q.SGD.1, y$q.Symbol.1, y$Row.1, y$Col.1, y$size_Ty.1, y$size_Scan.1, y$circ_Scan.1, y$GFP_Mean.1, y$RFP_Mean.1, y$GFP_Med.1, y$RFP_Med.1, y$Row.2, y$Col.2, y$size_Ty.2, y$size_Scan.2, y$circ_Scan.2, y$GFP_Mean.2, y$RFP_Mean.2, y$GFP_Med.2, y$RFP_Med.2,  y$Row.3, y$Col.3, y$size_Ty.3, y$size_Scan.3, y$circ_Scan.3, y$GFP_Mean.3, y$RFP_Mean.3, y$GFP_Med.3, y$RFP_Med.3, y$Row.4, y$Col.4, y$size_Ty.4, y$size_Scan.4, y$circ_Scan.4, y$GFP_Mean.4, y$RFP_Mean.4, y$GFP_Med.4, y$RFP_Med.4, y$medlogRatioNorm.1, y$medlogRatioNorm.2, y$medlogRatioNorm.3, y$medlogRatioNorm.4, y$meanlogRatioNorm.1, y$meanlogRatioNorm.2, y$meanlogRatioNorm.3, y$meanlogRatioNorm.4))
  names(New_form) <- c("Plate", "FLEX_quad", "Strain", "X384_Plates", "X384_Rows", "X384_Column", "SGD", "Symbol", "Row.1", "Col.1", "size_Ty.1", "size_Scan.1", "circ_Scan.1", "GFP_Mean.1", "RFP_Mean.1", "GFP_Med.1", "RFP_Med.1", "Row.2", "Col.2", "size_Ty.2", "size_Scan.2", "circ_Scan.2", "GFP_Mean.2", "RFP_Mean.2", "GFP_Med.2", "RFP_Med.2", "Row.3", "Col.3", "size_Ty.3", "size_Scan.3", "circ_Scan.3", "GFP_Mean.3", "RFP_Mean.3", "GFP_Med.3", "RFP_Med.3", "Row.4", "Col.4", "size_Ty.4", "size_Scan.4", "circ_Scan.4", "GFP_Mean.4", "RFP_Mean.4", "GFP_Med.4", "RFP_Med.4", "medlogRatioNorm.1", "medlogRatioNorm.2", "medlogRatioNorm.3", "medlogRatioNorm.4", "meanlogRatioNorm.1", "meanlogRatioNorm.2", "meanlogRatioNorm.3", "meanlogRatioNorm.4")
  
#####Now Add Avg mean logGFP/RFP avg median logGFP/RFP and Z score############
  Avg_medlogratioNorm <- rowMeans(New_form[,45:48], na.rm = TRUE)
  Avg_meanlogratioNorm <- rowMeans(New_form[,49:52], na.rm = TRUE)
  
  New_form_avgs <- cbind(data.frame(New_form, Avg_medlogratioNorm, Avg_meanlogratioNorm))  

  
  Sort_NA_New_form <- New_form_avgs[with(New_form_avgs, order(Avg_meanlogratioNorm)), ]
  Sort_NA_New_form <- subset(Sort_NA_New_form, Avg_meanlogratioNorm != "NA") 
  Sort_NA_New_form <- subset(Sort_NA_New_form, Avg_meanlogratioNorm != "NaN") 
  
  GFPRFPlog <- Sort_NA_New_form$Avg_meanlogratioNorm
  
  POPstddevGFPRFP_logmeanratio <- sd(GFPRFPlog)*sqrt((length(GFPRFPlog)-1)/(length(GFPRFPlog)))
  
  POP_logmeanGFPRFP <- mean(GFPRFPlog)
  
  
  Sort_NA_New_form$Z_LOG <- (Sort_NA_New_form$Avg_meanlogratioNorm - POP_logmeanGFPRFP)/POPstddevGFPRFP_logmeanratio

  
  write.table(Sort_NA_New_form, paste("Clean_", as.character(Col_files[n]), sep = ""), sep = "\t", row.names=F, quote=FALSE)
  
  
})

