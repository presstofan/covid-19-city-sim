library(raster)
library(rgdal)
library(tidyverse)
library(RANN)

# Read and Clean Data -----------------------------------------------------

# japan_pop <- raster("japan_map/jpn_ppp_2020.tif")
# # japan_pop@crs
# 
# # # crop the raster, then plot the new cropped raster
# cropbox1 <- extent(138.92915, 140.77516,
#                    34.97162, 36.24657)
# cropbox2 <- extent(139.55015933490043, 139.9607733485723,
#                    35.54205126667411, 35.84597272385362)
# tokyo <- crop(japan_pop, cropbox1)
# 
# # # lower the resolution to 1km
# #
# tokyo_1km <- aggregate(tokyo, fact=10, fun=sum)
# res(tokyo_1km)
# 
# tokyo_2km <- aggregate(tokyo, fact=20, fun=sum)
# res(tokyo_2km)
# 
# tokyo_3km <- aggregate(tokyo, fact=30, fun=sum)
# res(tokyo_3km)
# 
# tokyo_4km <- aggregate(tokyo, fact=40, fun=sum)
# res(tokyo_4km)
# 
# tokyo_5km <- aggregate(tokyo, fact=50, fun=sum)
# res(tokyo_5km)
# 
# writeRaster(tokyo_1km, filename = "japan_map/greater_tokyo_1km.tif", format="GTiff", overwrite=TRUE)
# writeRaster(tokyo_2km, filename = "japan_map/greater_tokyo_2km.tif", format="GTiff", overwrite=TRUE)
# writeRaster(tokyo_3km, filename = "japan_map/greater_tokyo_3km.tif", format="GTiff", overwrite=TRUE)
# writeRaster(tokyo_4km, filename = "japan_map/greater_tokyo_4km.tif", format="GTiff", overwrite=TRUE)
# writeRaster(tokyo_5km, filename = "japan_map/greater_tokyo_5km.tif", format="GTiff", overwrite=TRUE)

tokyo <- raster("japan_map/greater_tokyo_3km.tif")
res(tokyo)
# Pop <- values(tokyo)
# dfPop <- rasterToPoints(tokyo) %>% as.data.frame()
dfPop <- cbind(xyFromCell(tokyo, 1:ncell(tokyo)), values(tokyo)) %>% as.data.frame()
names(dfPop) <- c("x", "y", "pop")
dfPop$id <- 1:nrow(dfPop)
dfPop$pop <- replace_na(dfPop$pop, 0)
dfPop$pop <- round(dfPop$pop)

# clean OD flow data

# covert_sp <- function(path, filename){
# 
#   path <- "japan_map/S05-b-10_SYUTO_GML/"
#   filename <- "S05-b-10_SYUTO-1-g_PersonTripODAmount.shp"
# 
#   shp <- readOGR(dsn = paste0(path, filename))
#   dfFlow <- shp@data
#   temp <- raster::geom(shp) %>% as.data.frame()
#   temp_orig <- temp[seq(2,nrow(temp),2), c("x", "y")]
#   names(temp_orig) <- c("orig_x", "orig_y")
#   temp_dest <- temp[seq(1,nrow(temp),2),  c("x", "y")]
#   names(temp_dest) <- c("dest_x", "dest_y")
#   temp_comb <- cbind(temp_orig, temp_dest)
#   dfFlow <- cbind(dfFlow, temp_comb)
#   return(dfFlow)
# }
# 
# path <- "japan_map/S05-b-10_SYUTO_GML/"
# filename1 <- "S05-b-10_SYUTO-1-g_PersonTripODAmount.shp"
# filename2 <- "S05-b-10_SYUTO-2-g_PersonTripODAmount.shp"
# filename3 <- "S05-b-10_SYUTO-3-g_PersonTripODAmount.shp"
# 
# dfFlow1 <- covert_sp(path, filename1)
# dfFlow2 <- covert_sp(path, filename2)
# dfFlow3 <- covert_sp(path, filename3)
# 
# dfFlow <- rbind(dfFlow1, dfFlow2, dfFlow3)
# rm(dfFlow1, dfFlow2, dfFlow3)
# 
# saveRDS(dfFlow, "dfFlow.rds")

dfFlow <- readRDS("dfFlow.rds")

# get the total number of trip by transportation means

cols <- names(dfFlow)[5:35]
dfFlow[cols] <- sapply(dfFlow[cols],as.character)
dfFlow[cols] <- sapply(dfFlow[cols],as.numeric)

dfFlow_sub <- dfFlow %>%
  mutate(orig_code = as.character(S05b_003),
         dest_code = as.character(S05b_004),
         ntrip = S05b_035,
         ntrip_bus = S05b_016,
         ntrip_rail = S05b_010,
         ntrip_car = S05b_017,
         ntrip_public = ntrip_bus + ntrip_rail,
         ntrip_school = (S05b_006 + S05b_012 + S05b_018 + S05b_024 + S05b_030),
         ntrip_work = (S05b_005 + S05b_011 + S05b_017 + S05b_023 + S05b_029) +  # work
           S05b_008 + S05b_014 + S05b_020 + S05b_026 + S05b_032, # business trips
         ntrip_no_public = ntrip - ntrip_public,
         ntrip_no_rail = ntrip - ntrip_rail,
         ntrip_no_bus = ntrip - ntrip_bus,
         ntrip_no_car = ntrip - ntrip_car,
         ntrip_no_public_car = ntrip - ntrip_public - ntrip_car,
         ntrip_no_work = ntrip - ntrip_work,
         ntrip_no_work_school = ntrip - ntrip_work - ntrip_school) %>%
  select(orig_code, dest_code, 
         ntrip, ntrip_no_public, ntrip_no_rail, ntrip_no_bus, ntrip_no_car,
         ntrip_no_public_car, ntrip_no_work, ntrip_no_work_school,
         orig_x, orig_y, dest_x, dest_y)

# match the cells in the population raster file with the area code in
# the mobility OD flow file
area_coord <- dfFlow_sub %>% 
  select(orig_code, orig_x, orig_y) %>%
  filter(!duplicated(orig_code))

nearest <- nn2(data = area_coord[, c("orig_x", "orig_y")],
               query = dfPop[, c("x", "y")],
               searchtype = "standard", 
               k = 1)
sum(nearest[["nn.idx"]]==0)
dfPop$code <- as.character(area_coord[nearest$nn.idx, "orig_code"])

# creating the OD matrix
create_OD_matrix <- function(pop_df, flow_df, 
                             trip_var = "ntrip", filter_trip_var = NULL,
                             normalised = TRUE, reduce_OD = NULL) {
  # Debug
  # pop_df <- dfPop
  # flow_df <- dfFlow_sub
  # trip_var <- "ntrip"
  # filter_trip_var <- "ntrip_no_public"
  
  # use cartesian product to find all the combination of cell id
  # but here we use area code
  dfOD <- expand.grid(pop_df$code, pop_df$code) 
  names(dfOD) <- c("orig_code", "dest_code")
  dfOD$orig_code <- as.character(dfOD$orig_code)
  dfOD$dest_code <- as.character(dfOD$dest_code)
  dfOD <- left_join(dfOD, 
                    flow_df[, c("orig_code", "dest_code", trip_var, filter_trip_var)],
                    by = c("orig_code", "dest_code"))
  
  od_vector <- dfOD[, trip_var]
  od_vector[is.na(od_vector)] <- 0
  
  # now each cell in od_vector has the same flow number as the larger area, which is not right
  # we need to divide the flow number by number of unique flow between two areas
  # this is the average flow between a pair of cells from two areas
  # note that we also need to acount for the cells without population
  
  id_vector <- expand.grid(pop_df$id, pop_df$id)
  names(id_vector) <- c("orig_id", "dest_id")
  id_vector <- left_join(id_vector, pop_df[, c("id", "code")], by = c("orig_id" = "id"))
  names(id_vector) <- c("orig_id", "dest_id", "orig_code")
  id_vector <- left_join(id_vector, pop_df[, c("id", "code")], by = c("dest_id" = "id"))
  names(id_vector) <- c("orig_id", "dest_id", "orig_code", "dest_code")
  ## filter out the cells without population (e.g. water)
  id_vector <- left_join(id_vector, pop_df[, c("id", "pop")], by = c("orig_id" = "id"))
  id_vector <- left_join(id_vector, pop_df[, c("id", "pop")], by = c("dest_id" = "id"))
  ## get the number of duplicated orig_code dest_code match
  id_vector_by_orig_code <- id_vector %>%
    filter(pop.x != 0 & pop.y != 0)  %>%
    group_by(orig_code, dest_code) %>%
    summarise(n = n())
  scale_vector <- left_join(id_vector, id_vector_by_orig_code, by = c("orig_code", "dest_code")) %>%
    .[, "n"] # this will generate some NAs because the water-to-water flow
  
  if (is.null(filter_trip_var)){
    # calculated the average flow per cell and then scale it by row (origin)
    od_vector <- od_vector / scale_vector
    OD <- matrix(od_vector, nrow = nrow(pop_df), ncol = nrow(pop_df), byrow = FALSE)
    
    if (normalised == TRUE) {
      OD <- OD / rowSums(OD, na.rm = T) * pop_df$pop
      
      if (!is.null(reduce_OD)) {
        OD <- OD * reduce_OD
      }
    } 
    
    OD[is.na(OD)] <- 0
    return(OD)
    
  } else {
    # we still need to calculate the all trip OD matrix so that we can normalise the OD_filter
    # by the rowSum(OD), this act as a transportation alpha vector
    od_vector <- od_vector / scale_vector
    od_vector[is.na(od_vector)] <- 0
    OD <- matrix(od_vector, nrow = nrow(pop_df), ncol = nrow(pop_df), byrow = FALSE)
    
    od_vector_filter <- dfOD[, filter_trip_var]
    od_vector_filter[is.na(od_vector_filter)] <- 0
    od_vector_filter <- od_vector_filter / scale_vector
    OD_filter <- matrix(od_vector_filter, nrow = nrow(pop_df), ncol = nrow(pop_df), byrow = FALSE)
    
    if (normalised == TRUE) {
      OD_filter <- OD_filter / rowSums(OD, na.rm = T) * pop_df$pop # scale it using the all trip OD
      
      if (!is.null(reduce_OD)) {
        OD_filter <- OD_filter * reduce_OD
      }
      
      OD_filter[is.na(OD_filter)] <- 0
      return(OD_filter)
    }
  }
}

OD_alltrip <- create_OD_matrix(dfPop, dfFlow_sub, "ntrip")
OD_no_public <- create_OD_matrix(dfPop, dfFlow_sub, "ntrip", "ntrip_no_public")
OD_essential <- create_OD_matrix(dfPop, dfFlow_sub, "ntrip", reduce_OD = 0.1)

saveRDS(OD_alltrip, "OD_alltrip.rds")
saveRDS(OD_no_public, "OD_no_public.rds")
saveRDS(OD_essential, "OD_essential.rds")

