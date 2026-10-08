# Code 
# "Understanding coexistence mechanisms through functional niche partitioning 
# between Quercus and Lobatae oak sections" by 
# Goretty Mendoza-Juárez, Felipe García-Oliva, Fernando Pineda-García, Ricardo Gaytán-Legaria,
# Antonio González-Rodríguez
# Date: september 2026

# Libraries
library(pacman)
p_load(ggplot2, ggfortify, dplyr, corrplot, dplyr, TPD, hypervolume, gridExtra, 
       GGally, grid, vegan, ggvegan, purrr, tidyr, ecodist)

# Input/Output dir
indir  <- "data/"
outdir <- "results/"
DatosFull <- paste0(indir,"RawTraits.csv")
DatosFull <- read.csv(DatosFull, stringsAsFactors =TRUE)

# t-test Q. castanea and Q. obtusata
CastObt <- DatosFull %>% filter(SPECIES %in% c("Q. castanea", "Q. obtusata"))

FullZ <- CastObt %>% select(1,3, 9:22) %>%
  mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))

hist(FullZ$CHLOROPHYLL) 
shapiro.test(FullZ$CHLOROPHYLL)
t.test(CHLOROPHYLL ~ SPECIES, data= FullZ) 

hist(FullZ$PH) 
shapiro.test(FullZ$PH)
t.test(PH ~ SPECIES, data= FullZ) 
# ... rest of the traits 

#
# Boxplots (functional foliar traits)
aa <- ggplot(CastObt, aes(x = SPECIES, y = CHLOROPHYLL, fill = SPECIES)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#AA4B8F", "#FFC500")) +
  labs(x = "SPECIES", y = "CC (spad values)") + ggtitle("Chlorophyll concentration") + annotate("text", x = 0.75, y = 50.2, label = "**", size = 7) + theme(axis.text.x = element_text(face = "italic")) +
  theme(
    axis.title.x = element_text(size = 15),  
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 13),   
    axis.text.y = element_text(size = 13), 
    plot.title = element_text(size = 14.8, face = "bold"),
    legend.position = "none" 
  )
# ... rest of the traits 
#

#
# permanova by section and species
FullZ <- DatosFull %>% select(1, 2, 3, 9:22) %>%  
         mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))
colMeans(FullZ[, 4:ncol(FullZ)]); print(apply(FullZ, 2, sd)) 
# Distance matrix
FullZ_dist <- FullZ[, 4:ncol(FullZ)] 
FullZ_dist <- vegdist(FullZ_dist, method = "euclidean")
# Permanova
perman <- adonis2(FullZ_dist ~ SECTION, data = FullZ, permutations = 9999)

# permanova Q. cast and Q. obt
CastObt <- DatosFull %>% filter(SPECIES %in% c("Q. castanea", "Q. obtusata"))
FullZ <- CastObt %>% select(1, 2, 3, 9:22) %>%
         mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))
FullZ_dist <- FullZ[, 4:ncol(FullZ)] 
FullZ_dist <- vegdist(FullZ_dist, method = "euclidean")
perman <- adonis2(FullZ_dist ~ SPECIES, data = FullZ, permutations = 9999)

# permanova excluding Q. cast and Q.. obtusata
NoCastObt <- DatosFull %>% filter(!SPECIES %in% c("Q. castanea", "Q. obtusata"))
FullZ <- NoCastObt %>% select(1, 2, 3, 9:22) %>%
  mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))
FullZ_dist <- FullZ[, 4:ncol(FullZ)] 
FullZ_dist <- vegdist(FullZ_dist, method = "euclidean")
perman <- adonis2(FullZ_dist ~ SECTION, data = FullZ, permutations = 9999)

#
# Fric and functional dissimilarity by Section
# (FullZ data.- Col 4:CHLOROPHYLL - Col 17: N:P)
FullZ <- DatosFull %>% select(1, 2, 3, 9:22) %>%  
  mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))

TPDs <- TPDs(species = FullZ$SECTION, traits = FullZ[4]) # clorofila
RicEveDiv <- REND(TPDs = TPDs)
dissim_Secc <- dissim(TPDs) 

TPDs <- TPDs(species = FullZ$SECTION, traits = FullZ[5]) # pH
RicEveDiv <- REND(TPDs = TPDs) 
dissim_Secc <- dissim(TPDs)
# ... rest of the traits 

# Global functional dissimilarity by species pairs
FullZ <- DatosFull %>% select(1, 2, 3, 9:22) %>%  
  mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))

TPDs <- TPDs(species = FullZ$SPECIES, traits = FullZ[4]) # clorofila
dissim_Spp <- dissim(TPDs) 
# ... rest of the traits to construct the full matrix

# Fric and functional dissimilarity by Q. castanea and Q. obtusata
# (FullZ data.- Col 4:CLOROFILA - Col 17: N:P)
CastObt <- DatosFull %>% 
  filter(SPECIES %in% c("Q. castanea", "Q. obtusata")) %>% droplevels()

FullZ <- CastObt %>% select(1, 2, 3, 9:22) %>% droplevels() %>%
  mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))

TPDs <- TPDs(species = FullZ$SPECIES, traits = FullZ[4]) # Clorofila
RicEveDiv <- REND(TPDs = TPDs) 
dissim <- dissim(TPDs)
# ... rest of the traits 

#
# PCA using the full data
FullZ <- DatosFull %>% select(1, 2, 3, 9:22) %>%  
         mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))
dat_traits <- FullZ[c(4:14)] 
colnames(dat_traits) <- c("CC", "PH", "SLA", "LT", "WC", "SD", "SL", "TanC","TC", "TN", "TP")
pca <- prcomp(dat_traits, scale = TRUE, center =TRUE, retx=TRUE)
summary(pca)
pca$rotation[,1:3]
plotPca <- autoplot(pca, x = 1, y = 2, data = DatosFull, 
                    colour = "SECTION", shape = "SPECIES", size = 4,
                    loadings = TRUE, loadings.colour = '8', 
                    loadings.label.colour = 'black', loadings.label = TRUE,
                    loadings.label.repel = TRUE, loadings.label.size = 6) +
  scale_color_manual(values = c("#AA4B8F", "#FFC500")) +  
  labs(colour = "SECTION") +
  scale_shape_manual(values = c(15, 17, 7, 18, 19, 23, 11)) +
  geom_point(aes(shape = SPECIES, colour = SECTION),
             stroke = 2, size = 4) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold")
  )

# PCA excluding Q. cast y Q. obt
NoCastObt <- DatosFull %>% filter(!SPECIES %in% c("Q. castanea", "Q. obtusata"))
FullZ <- NoCastObt %>% select(1, 2, 3, 9:22) %>%
  mutate(across(where(is.numeric), ~ scale(., center = TRUE, scale = TRUE)))
# continue the PCA using the lines above

#
# functional hypervolumes 
# dataframe using the scores of the PCA 
ID <- DatosFull[1:3] %>% mutate(ID = 1:35) %>% dplyr::select(ID, SECTION, POP, SPECIES)
scores <- as.data.frame(pca$x[,1:3]) %>% mutate(ID = 1:35) %>% dplyr::select(ID, PC1, PC2, PC3)
scorespca <- left_join(ID, scores, by = "ID")

# Calculating by Section
SeccL_hv <- hypervolume(data = subset(scorespca, SECTION == "Lobatae")[,5:7], method = "gaussian", name=("Lobatae"))
SeccQ_hv <- hypervolume(data = subset(scorespca, SECTION == "Quercus")[,5:7], method = "gaussian", name =("Quercus"))
ListhvSecc <- hypervolume_join(SeccL_hv, SeccQ_hv)    
summary(ListhvSecc)
get_volume(ListhvSecc)
#plot
colors <- c("Lobatae" = "#AA4B8F", "Quercus" = "#FFC500")
plot(ListhvSecc, 
     limits = c(-5.5, 5.5), 
     show.contour = TRUE, contour.lwd = 4, 
     contour.type = 'kde', contour.kde.level = 1e-02, 
     show.centroid = TRUE, cex.centroid = 3.5, 
     colors = colors, 
     point.alpha.min = 0.5, point.dark.factor = 0.5, 
     cex.data = 1.2, 
     cex.axis = 1.8,  
     cex.names = 3,   
     cex.legend = 3.7) 
# Sorensen Overlap
hv_set <- hypervolume_set(SeccQ_hv, SeccL_hv, check.memory=FALSE)
hypervolume_overlap_statistics(hv_set)

# Calculating by Q. castanea and Q. obtusata (apply for the rest of the spp)
Q.cast_hv <- hypervolume(data = subset(scorespca, SPECIES == "Q. castanea")[,5:7], method = "gaussian", name=("Q. castanea"))
Q.obt_hv <- hypervolume(data = subset(scorespca, SPECIES == "Q. obtusata")[,5:7], method = "gaussian", name =("Q. obtusata"))
ListhvSpp <- hypervolume_join(Q.cast_hv, Q.obt_hv)    
summary(ListhvSpp)
get_volume(ListhvSpp)
# Sorensen Overlap
hv_set <- hypervolume_set(Q.cast_hv, Q.obt_hv, check.memory=FALSE)
hypervolume_overlap_statistics(hv_set)

# Rest of the species to build the paired matrix


# Climatic niche
t1=Sys.time()
# libraries
library(ecospat)
library(dplyr)
library(ade4)
library(raster)
library(rworldmap)
library(sp)
setwd(choose.dir())

# get climate data for the two area
clim_sp1<-stack(choose.files())
clim_sp1
clim_puntos_sp1<- rasterToPoints(clim_sp1[[1]], fun=NULL, spatial=TRUE)

clim_sp2<-stack(choose.files())
clim_sp2
clim_puntos_sp2<- rasterToPoints(clim_sp2[[1]], fun=NULL, spatial=TRUE)

# load occurrence sites for the species (column names should be x,y)
sp1<-as.data.frame(read.csv(choose.files()))
sp1_1 <- sp1%>%rename(x=longitude,y=latitude)
#sp1_1 <- sp1

sp2<-as.data.frame(read.csv(choose.files()))
sp2_1 <- sp2%>%rename(x=longitude,y=latitude)
#sp2_1 <- sp2
# remove duplicated occurrences
sp1_dup<-sp1_1[!duplicated(sp1_1),]
sp2_dup<-sp2_1[!duplicated(sp2_1),]

# remove occurrences closer than a minimum distance to each other (remove aggregation). 
#-Setting min.dist=0 will remove no occurrence. 
#-If your dataset is in degree, 0.008333 correspond to 1km at the equator.
#-Si busco eliminar registros a una diatancia de 4km, sería igual a 4X0.008333 = 0.033332

sp1_occ <- ecospat.occ.desaggregation(sp1_dup[,2:3],min.dist = 0.033332)
sp2_occ <- ecospat.occ.desaggregation(sp2_dup[,2:3],min.dist = 0.033332)

# create sp occurrence dataset by extracting climate variables from the rasters
env.sp1 <-  na.exclude(data.frame(extract(clim_sp1,clim_puntos_sp1)))
env.sp2 <-  na.exclude(data.frame(extract(clim_sp2,clim_puntos_sp2)))
env.occ.sp1 <- na.exclude(extract(clim_sp1,sp1_occ[,1:2]))
env.occ.sp2 <- na.exclude(extract(clim_sp2,sp2_occ[,1:2]))

#### niche quantifications #####
#calibration of PCA-env 
pca.env <-dudi.pca(rbind(env.sp1,env.sp2), center = T, scale = T, scannf = F, nf = 2)
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)

# predict the scores on the PCA axes
scores.sp<- pca.env$li
scores.sp1<- suprow(pca.env,env.sp1)$lisup
scores.sp2<- suprow(pca.env,env.sp2)$lisup
scores.occ.sp1<- suprow(pca.env,env.occ.sp1)$lisup
scores.occ.sp2<- suprow(pca.env,env.occ.sp2)$lisup

#write.csv(scores.occ.sp1
"C:/Users/Ricardo/Documents/Filogeografia y variacion genetica Q. glaucoides/Articulo/Revision_JBI/Minor revision/Comparacion_quercus/PCA_Quercus2.csv")

# calculation of occurence density
z1<- ecospat.grid.clim.dyn(scores.sp,scores.sp1,scores.occ.sp1,R=100)
z2<- ecospat.grid.clim.dyn(scores.sp,scores.sp2,scores.occ.sp2,R=100)
equ<-ecospat.niche.equivalency.test(z1,z2,rep=10, overlap.alternative = "lower")

# test of niche equivalency and similarity
sim1<-ecospat.niche.similarity.test(z1,z2,rep=100,overlap.alternative = "higher",rand.type = 1) #niches randomly shifted in both area
sim2<-ecospat.niche.similarity.test(z1,z2,rep=100,overlap.alternative = "higher",rand.type = 2) #niche randomly shifted only in invaded area

# overlap corrected by availabilty of background conditions
ecospat.niche.overlap(z1,z2,cor=T) 
# uncorrected overlap
ecospat.niche.overlap(z1,z2,cor=F) 

#### niche visualizations  #####
# occurrence density plots
ecospat.plot.niche(z1,title="PCA-env - Q_deserticola niche",name.axis1="PC1",name.axis2="PC2")
ecospat.plot.niche(z2,title="PCA-env - Q_glaucoides niche",name.axis1="PC1",name.axis2="PC2")

# contribution of original variables
ecospat.plot.contrib(pca.env$co,pca.env$eig)

# niche tests plots
ecospat.plot.overlap.test(equ,"D","Equivalency")
ecospat.plot.overlap.test(sim1,"D","Similarity 1<->2") #niches randomly shifted in both areas
ecospat.plot.overlap.test(sim2,"D","Similarity 1->2") #niche randomly shifted only in invaded area

ecospat.plot.niche.dyn(z1,z2,quant=0.5,name.axis1="PC1",name.axis2="PC2",interest=1)
ecospat.plot.niche.dyn(z1,z2,quant=0.5,name.axis1="PC1",name.axis2="PC2",interest=2)

# occurrence density in geography
plot(ecospat.niche.zProjGeo(z1,clim_sp1))
points(sp1_occ,cex=0.2,pch=19)
plot(ecospat.niche.zProjGeo(z2,clim_sp2))
points(sp2_occ,cex=0.2,pch=19)

#graficar las elipsoides 
# Installing and loading packages
if(!require(devtools)){
  install.packages("devtools")
}
if(!require(ellipsenm)){
  devtools::install_github("marlonecobos/ellipsenm")
}
library(ellipsenm)

####Niche overlap using ellipsenm
occurrences1<- read.csv(choose.files())
occurrences2<- read.csv(choose.files())
vars1 <- raster::stack(choose.files()) #seleccionar set de capas climaticas para la especie 1 en su M (seleccionar solo las capas ".asc")
vars2 <- raster::stack(choose.files()) #seleccionar set de capas climaticas para la especie 2 en su M (seleccionar solo las capas ".asc")

# preparing overlap objects to perform analyses
niche1 <- overlap_object(occurrences1, species =  "species", longitude = "longitude", 
                         latitude = "latitude", method = "covmat", level = 95, 
                         variables = vars1)
niche2 <- overlap_object(occurrences2, species =  "species", longitude = "longitude", 
                         latitude = "latitude", method = "covmat", level = 95, 
                         variables = vars2)
# niche overlap analysis
overlap <- ellipsoid_overlap(niche1, niche2, overlap_type = "full")
View(overlap)

# niche overlap analysis with test of significance
overlap_st <- ellipsoid_overlap(niche1, niche2, overlap_type = "back_union",
                                significance_test = TRUE, replicates = 100)
# plotting only ellipsoids
plot_overlap(overlap,  niche_col = c("blue","green"), 
             niches = c(1,2),
             data_col= c("blue","green") )
# plotting ellispodis and background for full overlap
plot_overlap(overlap, background = TRUE, proportion = 0.6, background_type = "full")
?plot_overlap
# plotting ellispodis and background for overlap based on accessible environments
plot_overlap(overlap, background = TRUE,  proportion = 1, background_type = "back_union")

t2=Sys.time()
t2-t1 


# Simpatry within my sites (following the methodology of Cannon (2024)
indir  <- "data/"
outdir <- "results/"
matrix <- paste0(indir,"SimpatryMatrix.csv")
matrix <- read.csv(matrix, stringsAsFactors =TRUE)

# getting the species names
species <- colnames(matrix)[-1]
# metrics
calculate_metrics <- function(sp1, sp2) {
  N_A <- sum(matrix[[sp1]], na.rm = TRUE)
  N_B <- sum(matrix[[sp2]], na.rm = TRUE)
  N_AB <- sum(matrix[[sp1]] == 1 & matrix[[sp2]] == 1, na.rm = TRUE)
  
  if (N_A > 0 & N_B > 0) {
    PS_A_to_B <- N_AB / N_A
    PS_B_to_A <- N_AB / N_B
    RS <- 1 - (N_AB / (N_A + N_B - N_AB))
  } else {
    PS_A_to_B <- NA
    PS_B_to_A <- NA
    RS <- NA}
  
  return(tibble(
    POPS_A = N_A,
    POPS_B = N_B,
    POPS_coocurren = N_AB,
    PS_A_to_B = round(PS_A_to_B, 4),
    PS_B_to_A = round(PS_B_to_A, 4),
    RS = round(RS, 4)))}

results <- combn(species, 2, simplify = FALSE) %>%
  map_df(~{
    sp1 <- .x[1]
    sp2 <- .x[2]
    metrics <- calculate_metrics(sp1, sp2)
    tibble(Sp1 = sp1, Sp2 = sp2) %>%
      bind_cols(metrics)
  }) %>%
  mutate(
    Interpretacion_RS = case_when(
      is.na(RS) ~ "Insufficient data",
      RS > 0.7 ~ "Strong segregation",
      RS < 0.3 ~ "High co-occurrence",
      TRUE ~ "Intermediate pattern"
    )
  ) %>%
  arrange(Sp1, Sp2)

# selecting certain columns 
results <- results %>%
  select(1, 2, 6, 7, 8) %>%   
  rename(PSAB = PS_A_to_B) %>%
  rename(PSBA = PS_B_to_A) %>%
  mutate(sympatry_maxv = pmax(PSAB, PSBA, na.rm = TRUE))


# Mantel test
indir  <- "data/"
outdir <- "results/"
data <- paste0(indir,"MantelTest.csv")
data <- read.csv(data, stringsAsFactors =TRUE)
# pairwise matrix 
data <- data %>%
  separate(Species, into = c("sp1", "sp2"), sep = "_", remove = FALSE)
especies <- sort(unique(c(data$sp1, data$sp2)))
n <- length(especies)
crear_matriz_simetrica <- function(data, col_valor, sp1_col = "sp1", sp2_col = "sp2") {
  especies <- sort(unique(c(data[[sp1_col]], data[[sp2_col]])))
  n <- length(especies)
  # empty matr
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(especies, especies))
  # data on it
  for(i in 1:nrow(data)) {
    esp1 <- data[[sp1_col]][i]
    esp2 <- data[[sp2_col]][i]
    valor <- data[[col_valor]][i]
    mat[esp1, esp2] <- valor
    mat[esp2, esp1] <- valor}
  diag(mat) <- 0
  return(mat)}
# choosing cols
cols <- c("functdiss", "overl_sorens", "D.index", "I.index", "sympatry_maxv",
                        "C.N","C.P","N_P","WC","SD","SL","TanC","Ct","Nt","Pt",
                        "CC","pH","SLA","LT")
# List
matrices <- lapply(cols, function(col) {
  crear_matriz_simetrica(data, col_valor = col)
})
names(matrices) <- cols
dists <- lapply(matrices, as.dist)

# functional global dissimilarity and climatic overlap
mantel(functdiss ~ I.index, data = dists, nperm = 9999)
mantel(functdiss ~ D.index, data = dists, nperm = 9999)
# sympatry and functional dissimilarity, functional and climatic overlap
mantel(sympatry_maxv ~ functdiss, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ overl_sorens, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ D.index, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ I.index, data = dists, nperm = 9999)
# sympatry and functional traits
mantel(sympatry_maxv ~ CC, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ pH, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ SLA, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ LT, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ WC, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ SD, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ SL, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ TanC, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ Ct, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ Nt, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ Pt, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ C.N, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ C.P, data = dists, nperm = 9999)
mantel(sympatry_maxv ~ N_P, data = dists, nperm = 9999)

a <- ggplot(data, aes(x = I.index, y = functdiss)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "CLIMATIC NICHE OVERLAP\n Hellinger's I index", y = "FUNCTIONAL DISSIMILITUDE") +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.5\n P = 0.06", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() + theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 11),  
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11))

b <- ggplot(data, aes(x = I.index, y = functdiss)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "CLIMATIC NICHE OVERLAP\n Schoener's D index", y = "FUNCTIONAL DISSIMILITUDE") +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.6\n P = 0.03", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() + theme(
             panel.background = element_rect(fill = "white"), 
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),  
             panel.border = element_rect(color = "black", fill = NA), 
             axis.text.x = element_text(size = 11),  # Etiquetas eje X más grandes
             axis.text.y = element_text(size = 11),  # Etiquetas eje Y más grandes
             axis.title.x = element_text(size = 11),
             axis.title.y = element_text(size = 11))
grid.arrange(a,b, nrow = 2, ncol = 1)

# ecological metrics
mantel(sympatry_maxv ~ functdiss, data = data, nperm = 9999)
mantel(sympatry_maxv ~ overl_sorens, data = data, nperm = 9999)
mantel(sympatry_maxv ~ I.index, data = data, nperm = 9999)
mantel(sympatry_maxv ~ D.index, data = data, nperm = 9999)

aa <- ggplot(data, aes(x = sympatry_maxv, y = functdiss)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = "GLOBAL FUNCTIONAL DISSIMILITUDE") +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.7\n P = 0.006", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 11),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 11),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

bb <- ggplot(data, aes(x = sympatry_maxv, y = overl_sorens)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = "HYPERVOLUME OVERLAP") +
  annotate("text", x = 0.98, y = 0.02, label = "r = 0.3\n P = 0.13", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 11),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 11),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

cc <- ggplot(data, aes(x = sympatry_maxv, y = D.index)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("CLIMATIC OVERLAP", "(Schoener'" * s * " D)"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = 0.5\n P = 0.06", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 11),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 11),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

dd <- ggplot(data, aes(x = sympatry_maxv, y = I.index)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("CLIMATIC OVERLAP", "(Hellinger'" * s * " I)"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = 0.5\n P = 0.03", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 11),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 11),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

grid.arrange(aa, bb, cc, dd, nrow = 2, ncol = 2)


# Mantel test foliar traits
a <- ggplot(data, aes(x = sympatry_maxv, y = CC)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "CHLOROPHYLL CONCENTRATION"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.6\nP = 0.03", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

b <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "pH"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.6\n P = 0.03", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

c <- ggplot(data, aes(x = sympatry_maxv, y = SLA)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "SPECIFIC LEAF AREA"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.7\n P = 0.006", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

d <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "LEAF THICKNESS"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.6\n P = 0.03", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

e <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "WATER CONTENT"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.5\n P = 0.06", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

f <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "STOMATAL DENSITY"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.2\n P = 0.18", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

g <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "STOMATAL PORE LENGTH"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.3\n P = 0.15", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

h <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "TANNIN CONTENT"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.4\n P = 0.09", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

i <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "TOTAL CARBON"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.3\n P = 0.18", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

j <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "TOTAL NITROGEN"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.07\n P = 0.37", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

k <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "TOTAL PHOSPHORUS"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.5\n P = 0.03", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

l <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "C.N RATIO"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.0003\n P = 0.49", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

m <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "C.P RATIO"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.5\n P = 0.05", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

n <- ggplot(data, aes(x = sympatry_maxv, y = pH)) + 
  geom_point(shape = 16, size = 3, color = "black") +
  labs(x = "SYMPATRY", y = expression(atop("FUNCTIONAL DISSIMILITUDE OF", "N.P RATIO"))) +
  annotate("text", x = 0.98, y = 0.02, label = "r = -0.4\n P = 0.06", vjust = 0,      
           hjust = 1, size = 4) + theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA), 
    axis.text.x = element_text(size = 8),  # Etiquetas eje X más grandes
    axis.text.y = element_text(size = 8),  # Etiquetas eje Y más grandes
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )
# plots for the rest of the traits
grid.arrange(a, b, c, d, e, f, g, h, i, j, k, nullGrob(), l, m, n, nullGrob(), nrow = 4, ncol = 4)




###