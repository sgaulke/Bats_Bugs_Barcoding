#' MIST NET SITES LAND COVER ANALYSES
#' NHD Data set Set-Up- Downloaded from here-
#' (https://gdg.sc.egov.usda.gov/GDGOrder.aspx?order=QuickState)
#' Land Cover Data set downloaded from here- 
#' (https://www.nass.usda.gov/Research_and_Science/Cropland/SARS1a.php)

#download and install needed packages
p <- c("sp", "raster", "rgdal", "rgeos", "foreign")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(p, require, character.only = TRUE)

##Set Working Directory to Personal Computer
#You will need to download the Spatial Data folder from Box under IBCP- Shared_Subfolders- IBCP_Field_Team- Guano Metabarcoding and set your working directory to where ever you downloaded and unzipped this folder 
wd <- "~/UIUC/Publications/Diet Metabarcoding/Spatial_Data/"

setwd(wd) 

##Load mist netting sites and reproject it
mn <- shapefile("Mist_Netting/MN_Diet.shp")
crs <- CRS("+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs ")
mn <- spTransform(mn, crs)

##Load hydrography data- str is streams, wb is water bodies, and ar is area of water (larger rivers are polygons etc.) and re-projec them
str <- shapefile("hydrography_NHD24K_il_3863606_01/hydrography/nhd24kst_l_il.shp")
str <- spTransform(str, crs)
wb <- shapefile("hydrography_NHD24K_il_3863606_01/hydrography/nhd24kwb_a_il.shp")
wb <- spTransform(wb, crs)
ar <- shapefile("hydrography_NHD24K_il_3863606_01/hydrography/nhd24kar_a_il.shp")
ar <- spTransform(ar, crs)


#Subset water sources to perennial sources that are biologically relevant. 
strsub <- str[str$FCODE == 46006,]
wbsub <- wb[(wb$FCODE == 39004 | 39009 | 39010 | 39011| 39012| 43615) & (wb$AREASQKM > 0),] 
arsub <- ar[(ar$FCODE == 31200 | 46006 |44500) & (ar$AREASQKM > 0), ]

##Land Cover Set-up for Cropscape
cdl <- raster(paste0(wd, "CDL_LandCover_il_2019_20201202151100/CDL_2019_clip_20201203161054_2030686002.tif"))
#read in and edit key to raster values
rat <- read.dbf(paste0(wd, "CDL_LandCover_il_2019_20201202151100/CDL_2019_clip_20201203161054_2030686002.tif.vat.dbf"))
rat <- na.omit(rat)
row.names(rat) <- 1:nrow(rat)
#'reclassify raster values to broader categories- 
#'1- Agriculture
#'2- Forest 
#'3- Urban 
#'4- Water 
#'5- Other
rc <- cbind(rat, c(0, rep(1,51),2,5,5, rep(1,11),5, 3, 4, 4, 5,rep(4,3), rep(3,4), 5, rep(2,3), 5, 1, 2, 4, rep(1,48)))
cat <- as.matrix(rc[,c(1,7)], nrow = 134, ncol = 2, byrow = FALSE)
lc <- reclassify(cdl, cat)

##This is a function to do the buffer analyses. The input for the function is the size of the buffer in meters. 
Buffer_Analysis <- function(size) {
  
  #add buffer around each mist net site
  mn_buff <- buffer(mn, size, dissolve = FALSE)
  
  #create data frames to write data into 
  per <- as.data.frame(c("Ag", "For", "Urb", "Water", "Oth"))
  dist <- as.data.frame(c("DistWat"))
  colnames(per) <- "Legend"
  
  #this for loop does the analysis for each site 
  for (i in 1:nrow(mn_buff)) {
    
    #specify site and crop land cover data to each buffer
    site <- mn_buff[i,]
    temp <- crop(lc, site)
    lr <- mask(temp, site)
    
    ##extract each cell value from each cell in raster 
    ext <- raster::extract(lr, site, method = "simple")
    
    #tabulate frequency of each value
    class <- lapply(ext, table)
    class2 <- as.data.frame(class)
    
    #calculate percentage of each value 
    for (j in 1:5) {
      class2[j,3] <- (class2[j,2])/(sum(class2[,2]))*100
    }
    #add data to dataframe and name each column to mist net site
    per <- cbind(per, class2[,3])
    colnames(per)[i+1] <- paste(site@data[5])
    
    ##Nearest Water- calculate closest distance to each water source
    str_temp <- gDistance(mn[i,], strsub)
    wb_temp <- gDistance(mn[i,], wbsub)
    ar_temp <- gDistance(mn[i,], arsub)
    
    #pick minimum distance to water source, add to data frame, and rename to mist net site
    dist <- cbind(dist, c(min(c(str_temp, wb_temp, ar_temp))))
    colnames(dist)[i+1] <- paste(site@data[5])
    colnames(dist)[1] <- "Legend"
    
  }
  #combine water and land cover values
  alldat2 <- rbind(per, dist)
  
  #create dataframe outside of function
  assign(paste0("dat", size/1000), alldat2, envir = .GlobalEnv)
}

##Check Land Cover Data for Accuracy
#create 10km buffer around mn
mn_buff <- buffer(mn, 10000, dissolve = FALSE)

#randomly add 100 points total spread out inside buffers 
sample <- spsample(mn_buff, n = 150, type = "stratified")

#extract cell value at each random point  
temp_lc <- extract(lc, sample, df = TRUE)

#write cell value and coordinates into dataframe, then export to Excel
temp2 <- cbind(as.data.frame(sample@coords), temp_lc)
write.table(temp2, "LandCoAccuracy_cdl.csv", sep = ",")
