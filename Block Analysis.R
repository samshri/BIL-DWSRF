library(tidyverse)
library(sf)
library(censusapi)
library(tigris)
options(tigris_use_cache = TRUE)
library(mapview)
Sys.setenv(CENSUS_KEY="02f0380ed2127343cb56f24b45b2737a454f0dbf")
#----POPULATION----------
pws_all <- st_read("/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/Final data/PWS_shapefile.gpkg")%>% 
  st_transform(26910)

pws <- pws_all %>% filter(state_code == 'WA')

# Retrieve census blocks for Washington state and transform to appropriate projection (EPSG:26910)
blocks <- blocks("WA", year = 2020) %>% 
  select(GEOID20) %>% 
  st_transform(26910)

# Retrieve county information for Washington state
wa_counties <- counties("WA")

# Retrieve population data for counties within Washington state
pop <- wa_counties$COUNTYFP %>% 
  map_dfr(function(county){
    tryCatch({
      getCensus(
        name = "dec/pl",
        vintage = 2020,
        region = "block:*", 
        regionin = paste0("state:53+county:", county),
        vars = c(
          "P1_001N",  # total population
          "P2_002N",  # hispanic or latino
          "P2_005N",  # white alone
          "P2_006N",  # black or african american alone
          "P2_007N",  # american indian and alaska native alone
          "P2_008N",  # asian alone
          "P2_009N"   # native hawaiian and other pacific islander alone
        )
      )
    }, error = function(e) NULL)  # Return NULL on error
  }) %>% 
  transmute(
    GEOID20 = paste0("53", county, tract, block),
    pop = P1_001N,
    hisp = P2_002N,
    white = P2_005N,
    black = P2_006N,
    aian = P2_007N,
    asian = P2_008N,
    nhpi = P2_009N,
    other = pop - hisp - white - black - aian - asian - nhpi
  )

# Join blocks with population data, filtering blocks with population > 0
blocks_pop <- blocks %>% 
  left_join(pop) %>% 
  filter(pop > 0)

# Extract unique pwsids for Washington state
ids <- unique(pws$pwsid)

# Initialize dataframe to store processed pws population data
wa_pws_pop <- data.frame()

# Loop through each pwsid and calculate adjusted population statistics
for (row in 1:length(ids)) {
  
  print(ids[row])
  
  temp_pws <- pws %>% 
    filter(pwsid == ids[row]) %>% 
    st_union()
  
  temp_blocks <- blocks_pop[temp_pws,] %>% 
    mutate(
      original_area = as.numeric(st_area(.))
    ) %>% 
    st_intersection(temp_pws) %>% 
    mutate(
      new_area = as.numeric(st_area(.)),
      perc_area = new_area / original_area
    ) %>% 
    mutate(across(
      pop:other,
      ~(.*perc_area),
      .names = "{.col}_adjust"
    )) %>% 
    st_drop_geometry() %>% 
    select(-GEOID20, -original_area, -new_area, -perc_area) %>% 
    summarize_all(sum) %>% 
    mutate(
      pwsid = ids[row]
    )
  
  wa_pws_pop <- wa_pws_pop %>% 
    rbind(temp_blocks)
  
}

# Define the full path to save the CSV file
save_path <- "/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/"

# Write the processed data to CSV at the specified location
write.csv(wa_pws_pop, file.path(save_path, "wa_pws_pop.csv"), row.names = FALSE)


#----INCOME----------

# Define the state abbreviation for California
state <- "MA"

# Fetch the boundary for California using tigris
state_boundary <- states(cb = TRUE) %>% 
  filter(STUSPS == state) %>% 
  st_transform(26910) # Transform to a common CRS (NAD83 / UTM zone 10N)

# Get the FIPS code for California
state_fp <- state_boundary %>% 
  pull(STATEFP)

state_counties <- counties(state = state, cb = TRUE)

# Fetch block group income data for California
income <- map_dfr(state_counties$COUNTYFP, function(county) {
  tryCatch(
    getCensus(
      name = "acs/acs5",
      vintage = 2021,
      region = "block group:*", 
      regionin = paste0("state:", state_fp, "+county:", county),
      vars = c(
        "B01001_001E",  # Total population
        "B19013_001E"   # Median Household Income
      )
    ),
    error = function(e) NULL
  )
}) %>% 
  transmute(
    GEOID = paste0(state_fp, county, tract, block_group),
    pop = B01001_001E,
    mhi = ifelse(B19013_001E > 0, B19013_001E, NA) # Validating MHI values
  )

# Get block group geometries for California using tigris
cbgs <- block_groups(state = state, year = 2020, cb = TRUE) %>% 
  select(GEOID) %>% 
  st_transform(26910) # Transform to a common CRS (NAD83 / UTM zone 10N)

# Join the income data with block group geometries
cbgs_mhi <- cbgs %>% 
  left_join(income, by = "GEOID") %>% 
  filter(pop > 0, !is.na(mhi)) # Filter out block groups with no population or MHI


# Filter PWS data for California
pws <- st_intersection(pws_all, state_boundary)

# Extract unique PWS IDs
ids <- unique(pws$pwsid)

# Initialize a dataframe to store PWS MHI calculations
pws_mhi <- data.frame()

# Iterate over each PWS and calculate weighted and unweighted MHI
for(row in seq_along(ids)){
  
  print(paste("Processing PWS:", row, "of", length(ids)))
  
  temp_pws <- pws %>% 
    filter(pwsid == ids[row]) %>% 
    st_union()
  
  temp_cbgs <- cbgs_mhi[temp_pws, ] %>% 
    st_drop_geometry() %>% 
    summarize(
      mhi_all_weighted = weighted.mean(mhi, pop, na.rm = TRUE),
      mhi_all_unweighted = mean(mhi, na.rm = TRUE)
    ) %>% 
    mutate(
      pwsid = ids[row]
    ) %>% 
    left_join(
      st_centroid(cbgs_mhi)[temp_pws, ] %>% 
        st_drop_geometry() %>% 
        summarize(
          mhi_centroid_weighted = weighted.mean(mhi, pop, na.rm = TRUE),
          mhi_centroid_unweighted = mean(mhi, na.rm = TRUE)
        ) %>% 
        mutate(
          pwsid = ids[row]
        ),
      by = "pwsid"
    )
  
  pws_mhi <- bind_rows(pws_mhi, temp_cbgs)
}

write.csv(pws_mhi, file.path(save_path, "ma_mhi.csv"), row.names = FALSE)

# Print the time taken for the operation
print(Sys.time() - start)

#----VOTE----------
vote  = read.csv('/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/2021blockgroupvoting_1.csv')
vote <- vote %>%rename("GEOID"="BLOCKGROUP_GEOID")
vote$GEOID <- as.character(vote$GEOID)

block_groups <- block_groups("WA", year = 2020) %>% 
  select(GEOID) %>% 
  st_transform(26910)

pws_all <- st_read("/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/Final data/PWS_shapefile.gpkg")%>% 
  st_transform(26910)
wa_pws = pws_all%>%filter(state_code=='WA')
blocks_vote <- block_groups %>% 
  left_join(vote) 

ids <- unique(wa_pws$pwsid)


wa_votes=data.frame()
for (row in 1:length(ids)) {
  print(row)
  temp_pws <- wa_pws %>% filter(pwsid == ids[row]) %>% st_union()
  temp_blocks <- blocks_vote[temp_pws,] %>% 
    mutate(
      original_area = as.numeric(st_area(.))
    ) %>% 
    st_intersection(temp_pws)%>% 
    mutate(
      new_area = as.numeric(st_area(.)),
      perc_area = new_area/original_area
    ) %>% 
    mutate(across(
      REP:OTH,
      ~(.*perc_area),
      .names = "{.col}_adjust"
    )) %>% 
    st_drop_geometry() %>% 
    select(-GEOID,-original_area,-new_area,-perc_area, -AREA, -GAP, -PRECINCTS, -STATE) %>% 
    summarize_all(sum) %>% 
    mutate(
      pwsid = ids[row]
    )
  wa_votes <- wa_votes %>% rbind(temp_blocks)
}
write.csv(wa_votes, "/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/wa_votes.csv")


#---------Merge some WA stuff---------
pws_all <- st_read("/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/Final data/PWS_shapefile.gpkg")%>% 
  st_transform(26910)

pws <- pws_all %>% filter(state_code == 'WA')

# Retrieve census blocks for Washington state and transform to appropriate projection (EPSG:26910)
tracts <- tracts("WA", year = 2010) %>%   st_transform(26910)
intersection <- st_intersection(pws, tracts)

# Calculate the area of intersection for each polygon
intersection <- intersection %>%
  mutate(intersection_area = st_area(.))

# Group by PWS and find the tract with the maximum intersection area
max_overlap <- intersection %>%
  group_by(pwsid) %>%
  filter(intersection_area == max(intersection_area)) %>%
  ungroup()

# Join the maximum overlap information back to the original PWS data

#write.csv(max_overlap, "/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/wa_pws_with_tract.csv")
wa_pop = read.csv("/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/wa_pws_pop.csv")
wa_votes = read.csv("/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/wa_votes.csv")
wa_mhi = read.csv("/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/wa_mhi.csv")
wa_svi = read.csv('/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/WA_SVI.csv')%>%rename("SVI"="Rank")
wa_ehd = read.csv('/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/WA_EHD.csv')%>%rename("EHD"="Rank")
wa_SVI_EHD = merge(wa_svi,wa_ehd,by="State.FIPS.Code" )%>%rename("FIPS"="State.FIPS.Code")
wa_data <- merge(max_overlap, wa_SVI_EHD, by.x = "GEOID10", by.y = "FIPS", all.x=TRUE)
wa_final <- wa_data %>% mutate(DAC=ifelse(SVI>=7|EHD>=7, 1, 0))

merged_data <- wa_pop %>%
  left_join(wa_votes, by = "pwsid")

# Merge the result with 'df3'
merged_data <- merged_data %>%
  left_join(wa_mhi, by = "pwsid")

# Finally, merge with the 'merged_pre_post' dataframe
wa_analysis <- merged_data %>%
  left_join(wa_final, by = "pwsid")
install.packages("writexl")  # Install if not already installed
library(writexl)             # Load the package
write_xlsx(wa_analysis,'/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/wa_analysis.xlsx')


#----Education - CA----
library(tidycensus)
library(dplyr)
library(tigris)
library(sf)

# Set state abbreviation
state <- "CA"

# Fetch the boundary for California using tigris
state_boundary <- states(cb = TRUE) %>%
  filter(STUSPS == state) %>%
  st_transform(26910) # Transform to a common CRS (NAD83 / UTM zone 10N)

# Get the FIPS code for California
state_fp <- state_boundary %>%
  pull(STATEFP)

# Get county information for the state
state_counties <- counties(state = state, cb = TRUE)

# Fetch block group educational attainment data for California
education <- map_dfr(state_counties$COUNTYFP, function(county) {
  tryCatch(
    get_acs(
      geography = "block group",
      state = state,
      county = county,
      variables = paste0("B15003_0", 17:25),  # High school diploma and above
      summary_var = "B15003_001"              # Population 25 years and older (denominator)
    ),
    error = function(e) NULL
  )
}) %>%
  group_by(GEOID) %>%
  summarize(
    n_hs_above = sum(estimate, na.rm = TRUE),
    n_hs_above_moe = moe_sum(moe, estimate),
    n_pop_over_25 = first(summary_est),
    n_pop_over_25_moe = first(summary_moe)
  ) %>%
  mutate(
    pct_hs_above = n_hs_above / n_pop_over_25,
    pct_hs_above_moe = moe_prop(n_hs_above, n_pop_over_25, n_hs_above_moe, n_pop_over_25_moe)
  ) %>%
  ungroup()

# Get block group geometries for California using tigris
cbgs <- block_groups(state = state, year = 2020, cb = TRUE) %>%
  select(GEOID) %>%
  st_transform(26910) # Transform to a common CRS (NAD83 / UTM zone 10N)

# Join the educational data with block group geometries
cbgs_education <- cbgs %>%
  left_join(education, by = "GEOID") %>%
  filter(!is.na(n_pop_over_25)) # Filter out groups with missing population data

# Filter PWS data for California
pws <- st_intersection(pws_all, state_boundary)

# Extract unique PWS IDs
ids <- unique(pws$pwsid)

# Initialize a dataframe to store PWS education calculations
pws_education <- data.frame()

# Iterate over each PWS and calculate weighted and unweighted high school attainment
for(row in seq_along(ids)){
  
  print(paste("Processing PWS:", row, "of", length(ids)))
  
  temp_pws <- pws %>%
    filter(pwsid == ids[row]) %>%
    st_union()
  
  temp_cbgs <- cbgs_education[temp_pws, ] %>%
    st_drop_geometry() %>%
    summarize(
      pct_hs_above_weighted = weighted.mean(pct_hs_above, n_pop_over_25, na.rm = TRUE),
      pct_hs_above_unweighted = mean(pct_hs_above, na.rm = TRUE)
    ) %>%
    mutate(
      pwsid = ids[row]
    )
  
  pws_education <- bind_rows(pws_education, temp_cbgs)
}
pws_education <- pws_education %>% mutate(pct_hs_below_weighted = 100-(pct_hs_above_weighted*100))%>%select(pwsid, pct_hs_below_weighted)
save_path <- "/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/Final Data/Age_and_education"



# Save the results to CSV
write.csv(pws_education, file.path(save_path, "ca_education.csv"), row.names = FALSE)


