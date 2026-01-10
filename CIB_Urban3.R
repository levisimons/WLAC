rm(list=ls())
require(data.table)
require(dplyr)
require(phyloseq)
require(sf)
require(ggplot2)
require(vegan)
require(emmeans)
require(tidyr)

wd <- "~/Desktop/Archive/CIB/"
setwd(wd)

set.seed(1)

#Set sf settings to reduce risk of errors
sf_use_s2(FALSE)

#Read in Metropolitan Statistical Area (MSA) shapefiles for 2023
#https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2023&layergroup=Urban+Areas
MSA_Boundaries <- sf::st_read(paste(wd,"/tl_2023_us_uac20/tl_2023_us_uac20.shp",sep=""))
#Reproject
MSA_Boundaries <- sf::st_transform(MSA_Boundaries,crs=3857)
#Read in Metropolitan Statistical Area (MSA) shapefiles for 2000
#https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2000&layergroup=Urban+Areas
#MSA_Boundaries <- sf::st_read(paste(wd,"/tl_2008_us_uac00/tl_2008_us_uac00.shp",sep=""))
#Reproject
#MSA_Boundaries <- sf::st_transform(MSA_Boundaries,crs=3857)

#Get boundaries of California
require(rnaturalearth)
# Get US states
usa <- ne_states(country = "United States of America", returnclass = "sf")
# Subset California
california <- usa[usa$name == "California", ]

#GBIF.org (5 December 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.pnhtw9
gbif_input <- fread(input="0030378-251120083545085.zip")
gbif_spatial <- st_as_sf(gbif_input, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
#Reproject
gbif_spatial <- sf::st_transform(gbif_spatial,crs=3857)

#Intersect urban and non-urban areas with species occurrences.
Sampled_Taxa <- st_join(gbif_spatial,MSA_Boundaries)
Sampled_Taxa$Area <- ifelse(is.na(Sampled_Taxa$NAME20),"Non-Urban","Urban")
tmp <- as.data.frame(st_coordinates(Sampled_Taxa))
colnames(tmp) <- c("decimalLongitude","decimalLatitude")
Sampled_Taxa <- st_drop_geometry(Sampled_Taxa)
Sampled_Taxa <- cbind(Sampled_Taxa,tmp)
setDT(Sampled_Taxa)
Sampled_Taxa[, sampleid := .GRP, by = .(eventDate,decimalLongitude,decimalLatitude)]

gbif_filtered <- Sampled_Taxa %>%
  group_by(sampleid) %>%
  filter(n() >= 1) %>%
  ungroup()

gbif_filtered <- gbif_filtered[gbif_filtered$species!="",]

taxa_count <- gbif_filtered[,c("sampleid","taxonKey")]
taxa_count <- taxa_count %>%
  group_by(sampleid, taxonKey) %>%
  mutate(taxonCount = n()) %>%
  ungroup()
taxa_count <- taxa_count[!duplicated(taxa_count),]

#Create a presence/absence data table.
gbif_pa <- taxa_count %>%
  pivot_wider(
    names_from = sampleid,
    values_from = taxonCount,
    values_fill = 0   # fill missing combinations with 0 (or NA if you prefer)
  )

#Create an otu matrix
otu_mat <- as.data.frame(gbif_pa)
rownames(otu_mat) <- otu_mat$taxonKey
otu_mat$taxonKey <- NULL
otu_mat <- as.matrix(otu_mat)
OTU <- otu_table(otu_mat,taxa_are_rows=T)
#Create a taxonomy matrix
tax_mat <- st_drop_geometry(gbif_filtered[,c("kingdom","phylum","class","order","family","genus","species","taxonKey")])
tax_mat <- tax_mat[!duplicated(tax_mat),]
TAX <- tax_table(as.matrix(tax_mat))
taxa_names(TAX) <- tax_mat$taxonKey
#Create a metadata matrix
sample_mat <- st_drop_geometry(gbif_filtered[,c("Area","sampleid","year")])
sample_mat <- sample_mat[!duplicated(sample_mat),]
rownames(sample_mat) <- sample_mat$sampleid
Sample <- sample_data(sample_mat)
sample_names(Sample) <- sample_mat$sampleid

#Create full phyloseq object
physeq_initial <- phyloseq(OTU,TAX,Sample)

fungal_phyla <- c("Ascomycota","Basidiomycota","Entorrhizomycetes","Blastocladiomycota","Chytridiomycota","Cryptomycota","Microsporidia","Mucoromycota","Nephridiophaga","Olpidiomycota","Sanchytriomycota","Zoopagomycota")
invertebrate_phyla <- c("Platyhelminthes","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Annelida","Mollusca","Arthropoda")
fungal_orders <- unique(na.omit(gbif_filtered[gbif_filtered$phylum %in% fungal_phyla,"order"]))
fungal_orders <- fungal_orders[fungal_orders$order!="",]$order
invertebrate_orders <- unique(na.omit(gbif_filtered[gbif_filtered$phylum %in% invertebrate_phyla,"order"]))
invertebrate_orders <- invertebrate_orders[invertebrate_orders$order!="",]$order

#Summarize alpha and beta diversity statistics between urban and non-urban environments by order.
i <- 1
j <- 1
alpha_summary <- c()
beta_summary <- c()
for(order in c(fungal_orders,invertebrate_orders)){
  #Filter phyloseq object by taxonomy
  ps_filtered <- prune_taxa(
    tax_table(physeq_initial)[, "order"] %in% c(order),
    physeq_initial
  )
  ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
  #Compare urban versus non-urban richness
  richness_df <- estimate_richness(ps_filtered, measures = "Chao1")
  richness_df$sampleid <- sample_names(ps_filtered)
  richness_df <- merge(richness_df, sample_mat, by = "sampleid")
  #Create a glm with a gamma distribution of chao-1 as a function of the interaction between substrate and urbanization (family is gamma and link equals log). Run an anova on this glm.
  if(length(unique(richness_df$Area))>1){
    richness_model <- glm(Chao1 ~ factor(Area),data=richness_df,family="Gamma"(link='log'))
    em_richness <- emmeans(richness_model, pairwise ~ Area)
    em_contrasts <- as.data.frame(summary(em_richness$contrasts))
    em_contrasts$order <- order
    if(order %in% fungal_orders){em_contrasts$type <- "fungi"}
    if(order %in% invertebrate_orders){em_contrasts$type <- "invertebrates"}
    alpha_summary[[i]] <- em_contrasts
    print(paste("alpha",order,i,length(c(fungal_orders,invertebrate_orders))))
    i=i+1
  }
  if(min(table(sample_data(ps_filtered)$Area))>1 & length(unique(sample_data(ps_filtered)$Area))>1){
    #Check if beta diversity distributions are significantly different between urban and non-urban areas
    beta_dist <- phyloseq::distance(ps_filtered, "chao")
    tmp <- as.data.frame(TukeyHSD(betadisper(beta_dist,group=sample_data(ps_filtered)$Area,bias.adjust=T))$group)
    tmp$order <- order
    if(order %in% fungal_orders){tmp$type <- "fungi"}
    if(order %in% invertebrate_orders){tmp$type <- "invertebrates"}
    urbanization <- data.frame(t(data.frame(table(sample_data(ps_filtered)$Area))))
    colnames(urbanization) <- urbanization[1,]
    urbanization <- urbanization[-1,]
    tmp <- cbind(tmp,urbanization)
    beta_summary[[j]] <- tmp
    print(paste("beta",order,j,length(c(fungal_orders,invertebrate_orders))))
    j=j+1
  }
  
}
alpha_summary <- rbindlist(alpha_summary)
beta_summary <- rbindlist(beta_summary)
write.table(x=alpha_summary,file="alpha_summary.txt",sep="\t",row.names=F)
write.table(x=beta_summary,file="beta_summary.txt",sep="\t",row.names=F)
alpha_summary <- fread(input="alpha_summary.txt",sep="\t")
beta_summary <- fread(input="beta_summary.txt",sep="\t")
alpha_summary_filtered <- alpha_summary[alpha_summary$p.value<=0.05,]
beta_summary_filtered <- beta_summary[beta_summary$`p adj`<=0.05,]

#Compare the importance of sample year versus urban/non-urban status on alpha diversity split by order
alpha_summary <- c()
i=1
for(order in c(fungal_orders,invertebrate_orders)){
  ps_filtered <- prune_taxa(
    tax_table(physeq_initial)[, "order"] %in% c(order),
    physeq_initial
  )
  ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
  richness_df <- estimate_richness(ps_filtered, measures = "Chao1")
  richness_df$sampleid <- sample_names(ps_filtered)
  richness_df <- merge(richness_df, sample_mat, by = "sampleid")
  if(length(unique(richness_df$Area))>1){
    richness_model <- glm(Chao1 ~ factor(Area)*year,data=richness_df,family="Gamma"(link='log'))
    if(sum(is.na(richness_model$coefficients))<=0){
      tmp <- anova(richness_model,test="F")
      tmp$factor <- rownames(tmp)
      tmp$order <- order
      if(order %in% fungal_orders){tmp$community <- "Fungi"}
      if(order %in% invertebrate_orders){tmp$community <- "Invertebrates"}
      tmp <- tmp[2:4,]
      alpha_summary[[i]] <- tmp
      print(paste(order,i,nrow(tmp),length(c(fungal_orders,invertebrate_orders))))
      i=i+1
    }
  }
}
alpha_summary <- rbindlist(alpha_summary)
write.table(x=alpha_summary,file="alpha_order_anova.txt",sep="\t",row.names=F)
alpha_summary_filtered <- alpha_summary[alpha_summary$`Pr(>F)`<=0.05,]
mean(alpha_summary_filtered[alpha_summary_filtered$community=="Fungi" & alpha_summary_filtered$factor=="factor(Area)","F"]$F)
mean(alpha_summary_filtered[alpha_summary_filtered$community=="Fungi" & alpha_summary_filtered$factor=="year","F"]$F)

#Check temporal trends in alpha/beta diversity
i <- 1
j <- 1
alpha_summary <- c()
beta_summary <- c()
for(year_selected in unique(sample_data(physeq_initial)$year)){
  #Filter phyloseq object by year and taxonomic group
  ps_filtered <- prune_taxa(
    tax_table(physeq_initial)[, "phylum"] %in% fungal_phyla,
    physeq_initial
  )
  ps_filtered <- subset_samples(ps_filtered, as.numeric(as.character(year)) == as.numeric(as.character(year_selected)))
  ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
  #Compare urban versus non-urban richness
  richness_df <- estimate_richness(ps_filtered, measures = "Chao1")
  richness_df$sampleid <- sample_names(ps_filtered)
  richness_df <- merge(richness_df, sample_mat, by = "sampleid")
  #Create a glm with a gamma distribution of chao-1 as a function of the interaction between substrate and urbanization (family is gamma and link equals log). Run an anova on this glm.
  if(length(unique(richness_df$Area))>1){
    richness_model <- glm(Chao1 ~ factor(Area),data=richness_df,family="Gamma"(link='log'))
    em_richness <- emmeans(richness_model, pairwise ~ Area)
    em_contrasts <- as.data.frame(summary(em_richness$contrasts))
    em_contrasts$year <- year_selected
    if(sum(tax_table(ps_filtered)[,"phylum"] %in% fungal_phyla)>1){em_contrasts$type <- "fungi"}
    if(sum(tax_table(ps_filtered)[,"phylum"] %in% invertebrate_phyla)>1){em_contrasts$type <- "invertebrates"}
    alpha_summary[[i]] <- em_contrasts
    print(paste("alpha",year_selected,i,length(unique(sample_data(physeq_initial)$year))))
    i=i+1
  }
  if(min(table(sample_data(ps_filtered)$Area))>1 & length(unique(sample_data(ps_filtered)$Area))>1){
    #Check if beta diversity distributions are significantly different between urban and non-urban areas
    beta_dist <- phyloseq::distance(ps_filtered, "chao")
    tmp <- as.data.frame(TukeyHSD(betadisper(beta_dist,group=sample_data(ps_filtered)$Area,bias.adjust=T))$group)
    tmp$year <- year_selected
    if(sum(tax_table(ps_filtered)[,"phylum"] %in% fungal_phyla)>1){tmp$type <- "fungi"}
    if(sum(tax_table(ps_filtered)[,"phylum"] %in% invertebrate_phyla)>1){tmp$type <- "invertebrates"}
    urbanization <- data.frame(t(data.frame(table(sample_data(ps_filtered)$Area))))
    colnames(urbanization) <- urbanization[1,]
    urbanization <- urbanization[-1,]
    tmp <- cbind(tmp,urbanization)
    beta_summary[[j]] <- tmp
    print(paste("beta",year_selected,j,length(unique(sample_data(physeq_initial)$year))))
    j=j+1
  }
  
}
alpha_summary <- rbindlist(alpha_summary)
beta_summary <- rbindlist(beta_summary)
write.table(x=alpha_summary,file="alpha_fungi_time.txt",sep="\t",row.names=F)
write.table(x=beta_summary,file="beta_fungi_time.txt",sep="\t",row.names=F)
alpha_summary <- fread(input="alpha_fungi_time.txt",sep="\t")
beta_summary <- fread(input="beta_fungi_time.txt",sep="\t")

#Variations in beta diversity as a function of urbanization and year per order
beta_summary <- c()
i=1
for(order in c(fungal_orders,invertebrate_orders)){
  ps_filtered <- prune_taxa(
    tax_table(physeq_initial)[, "order"] %in% c(order),
    physeq_initial
  )
  ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
  if(min(table(sample_data(ps_filtered)$Area))>1 & length(unique(sample_data(ps_filtered)$Area))>1){
    beta_dist <- phyloseq::distance(ps_filtered, "chao")
    tmp <- adonis2(beta_dist ~ Area*year,data=data.frame(sample_data(ps_filtered)))
    tmp$Factor <- rownames(tmp)
    tmp$order <- order
    if(order %in% fungal_orders){tmp$community <- "Fungi"}
    if(order %in% invertebrate_orders){tmp$community <- "Invertebrates"}
    tmp <- tmp[1:3,]
    beta_summary[[i]] <- tmp
    print(paste(order,i,length(c(fungal_orders,invertebrate_orders))))
    i <- i+1
  }
}
beta_summary <- rbindlist(beta_summary)
write.table(x=beta_summary,file="beta_orders.txt",sep="\t",row.names=F)

require(SpadeR)
#Calculate Observed richness
ps_filtered <- prune_taxa(
  tax_table(physeq_initial)[, "phylum"] %in% invertebrate_phyla,
  physeq_initial
)
ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
richness_data <- plot_richness(ps_filtered,measures=c("Observed"))$data
colnames(richness_data)[which(names(richness_data) == "value")] <- "Observed"
#Calculate Chao-1 richness
tmp <- plot_richness(ps_filtered,measures=c("Chao1"))$data
colnames(tmp)[which(names(tmp) == "value")] <- "Chao1"
richness_data <- dplyr::left_join(richness_data,tmp[,c("sampleid","Chao1")])
#Calculate iChao1 and Chao1-bc richness
otu_table_retain <- data.frame(otu_table(ps_filtered))
iChao_richness_values <- list()
i=1
for(sample in colnames(otu_table_retain)){
  iChao_richness_value <- data.frame(matrix(nrow=1,ncol=3))
  colnames(iChao_richness_value) <- c("Sample.ID","iChao1","Chao1-bc")
  iChao_richness_value$Sample.ID <- sample
  if(sum(otu_table_retain[,sample]>1)>4){
    tmp <- as.data.frame(ChaoSpecies(otu_table_retain[,sample],datatype="abundance")$Species_table)
    iChao_richness_value$iChao1 <- tmp[5,"Estimate"]
    iChao_richness_value$`Chao1-bc` <- tmp[4,"Estimate"]
  } else{
    iChao_richness_value$iChao1 <- NA
    iChao_richness_value$`Chao1-bc` <- NA
  }
  
  iChao_richness_values[[i]] <- iChao_richness_value
  i=i+1
  print(i)
}
iChao_richness_values <- rbindlist(iChao_richness_values)
iChao_richness_values$sampleid <- as.numeric(sub("^X", "", iChao_richness_values$Sample.ID))
richness_data <- dplyr::left_join(richness_data,iChao_richness_values)
#Export richness data
write.table(richness_data,"Invertebrate_Richness.txt", row.names=FALSE, sep="\t",quote = FALSE)
richness_data <- fread(input="Invertebrate_Richness.txt",sep="\t")
#Visualize correlations between richness values
require(PerformanceAnalytics)
chart.Correlation(richness_data[,c("Chao1","iChao1","Chao1-bc")])

#Compare LCBD values between urban and non-urban samples.
require(adespatial)
LCBD_summary <- c()
i=1
for(order in c(fungal_orders,invertebrate_orders)){
  ps_filtered <- prune_taxa(
    tax_table(physeq_initial)[, "order"] %in% order,
    physeq_initial
  )
  ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
  if(nrow(tax_table(ps_filtered))>1){
    beta_dist <- phyloseq::distance(ps_filtered, "chao")
    LCDB_samples <- LCBD.comp(beta_dist,save.D = TRUE)
    metadata_order <- data.frame(sample_data(ps_filtered))
    metadata_order$LCBD <- LCDB_samples$LCBD
    if(length(unique(metadata_order$Area))>1){
      LCBD_model <- glm(LCBD ~ factor(Area),data=metadata_order,family="Gamma"(link='log'))
      em_LCBD <- emmeans(LCBD_model, pairwise ~ Area)
      em_LCBD_contrasts <- as.data.frame(summary(em_LCBD$contrasts))
      em_LCBD_contrasts$order <- order
      if(order %in% fungal_orders){em_LCBD_contrasts$Community <- "Fungi"}
      if(order %in% invertebrate_orders){em_LCBD_contrasts$Community <- "Invertebrates"}
      LCBD_summary[[i]] <- em_LCBD_contrasts
      print(paste(i,length(c(fungal_orders,invertebrate_orders))))
      i=i+1
    }
  }
}
LCBD_summary <- rbindlist(LCBD_summary)
write.table(LCBD_summary,"Orders_LCBD.txt", row.names=FALSE, sep="\t",quote = FALSE)
LCBD_summary_filtered <- LCBD_summary[LCBD_summary$p.value<=0.05,]

#Get the number of exclusively urban or non-urban taxa.
urban_fungi <- unique(gbif_filtered[gbif_filtered$Area=="Urban" & gbif_filtered$phylum %in% fungal_phyla,"species"]$species)[!(unique(gbif_filtered[gbif_filtered$Area=="Urban" & gbif_filtered$phylum %in% fungal_phyla,"species"]$species) %in% unique(gbif_filtered[gbif_filtered$Area=="Non-Urban" & gbif_filtered$phylum %in% fungal_phyla,"species"]$species))]
non_urban_fungi <- unique(gbif_filtered[gbif_filtered$Area=="Non-Urban" & gbif_filtered$phylum %in% fungal_phyla,"species"]$species)[!(unique(gbif_filtered[gbif_filtered$Area=="Non-Urban" & gbif_filtered$phylum %in% fungal_phyla,"species"]$species) %in% unique(gbif_filtered[gbif_filtered$Area=="Urban" & gbif_filtered$phylum %in% fungal_phyla,"species"]$species))]
urban_invertebrates <- unique(gbif_filtered[gbif_filtered$Area=="Urban" & gbif_filtered$phylum %in% invertebrate_phyla,"species"]$species)[!(unique(gbif_filtered[gbif_filtered$Area=="Urban" & gbif_filtered$phylum %in% invertebrate_phyla,"species"]$species) %in% unique(gbif_filtered[gbif_filtered$Area=="Non-Urban" & gbif_filtered$phylum %in% invertebrate_phyla,"species"]$species))]
non_urban_invertebrates <- unique(gbif_filtered[gbif_filtered$Area=="Non-Urban" & gbif_filtered$phylum %in% invertebrate_phyla,"species"]$species)[!(unique(gbif_filtered[gbif_filtered$Area=="Non-Urban" & gbif_filtered$phylum %in% invertebrate_phyla,"species"]$species) %in% unique(gbif_filtered[gbif_filtered$Area=="Urban" & gbif_filtered$phylum %in% invertebrate_phyla,"species"]$species))]

#Map species richness
require(leaflet)
require(leaflegend)
ps_filtered <- prune_taxa(
  tax_table(physeq_initial)[, "order"] %in% fungal_orders,
  physeq_initial
)
ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
richness_df <- estimate_richness(ps_filtered, measures = "Chao1")
richness_df$sampleid <- sample_names(ps_filtered)
richness_df <- merge(richness_df, gbif_filtered[,c("sampleid","decimalLongitude","decimalLatitude")], by = "sampleid")
richness_df <- richness_df[,c("Chao1","decimalLongitude", "decimalLatitude")]
richness_df <- richness_df[!duplicated(richness_df),]
richness_points <- st_as_sf(richness_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 3857)
richness_points <- st_transform(richness_points,crs=4326)
#
MSA_Boundaries <- sf::st_transform(MSA_Boundaries,crs=4326)
MSA_California <- sf::st_filter(MSA_Boundaries,california)
# Define a color palette based on Chao1
pal <- colorNumeric("plasma", log10(richness_points$Chao1))
#Map richness
leaflet() %>%
  addProviderTiles(providers$OpenStreetMap) %>%     # OSM basemap
  setView(lng = -119.5, lat = 37.25, zoom = 6) %>%  # California center
  addPolygons(
    data = MSA_California,
    color = "black",         # border color
    weight = 2,             # border width
    fillColor = "magenta",# fill color
    fillOpacity = 0.8            # optional, a column from shapefile
  ) %>%
  addCircleMarkers(
    lng = st_coordinates(richness_points)[,1],
    lat = st_coordinates(richness_points)[,2],
    radius = 3*(0.5+log10(richness_points$Chao1)),                # adjust or scale by Chao1 if desired
    color = pal(log10(richness_points$Chao1)),
    fillOpacity = 0.8,
    stroke = FALSE
  ) %>%
  addLegendNumeric(
    position="bottomright",
    pal = pal,
    values = log10(richness_points$Chao1),
    title = "Log(Chao1)",
    decreasing=F)

#Map LCBD for fungi and invertebrates
ps_filtered <- prune_taxa(
  tax_table(physeq_initial)[, "order"] %in% invertebrate_orders,
  physeq_initial
)
ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
beta_dist <- phyloseq::distance(ps_filtered, "chao")
LCDB_samples <- LCBD.comp(beta_dist,save.D = TRUE)
metadata_order <- data.frame(sample_data(ps_filtered))
metadata_order$LCBD <- LCDB_samples$LCBD
richness_df <- merge(metadata_order, gbif_filtered[,c("sampleid","decimalLongitude","decimalLatitude")], by = "sampleid")
richness_df <- richness_df[,c("LCBD","decimalLongitude", "decimalLatitude")]
richness_df <- richness_df[!duplicated(richness_df),]
richness_points <- st_as_sf(richness_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 3857)
richness_points <- st_transform(richness_points,crs=4326)
# Define a color palette based on Chao1
pal <- colorNumeric("plasma", log10(richness_points$LCBD))
#Map richness
leaflet() %>%
  addProviderTiles(providers$OpenStreetMap) %>%     # OSM basemap
  setView(lng = -119.5, lat = 37.25, zoom = 6) %>%  # California center
  addPolygons(
    data = MSA_California,
    color = "black",         # border color
    weight = 2,             # border width
    fillColor = "magenta",# fill color
    fillOpacity = 0.8            # optional, a column from shapefile
  ) %>%
  addCircleMarkers(
    lng = st_coordinates(richness_points)[,1],
    lat = st_coordinates(richness_points)[,2],
    radius = 2,                # adjust or scale by Chao1 if desired
    color = pal(log10(richness_points$LCBD)),
    fillOpacity = 0.8,
    stroke = FALSE
  ) %>%
  addLegendNumeric(
    position="bottomright",
    pal = pal,
    values = log10(richness_points$LCBD),
    title = "Log(LCBD)",
    decreasing=F)
