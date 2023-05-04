## Script for making soil texture ternary plot

## Tutorial from: https://saryace.github.io/flipbook_soiltexture_en/#36, https://stackoverflow.com/questions/29136168/add-texture-classes-to-soil-classification-triangle-via-ggplot2 
## Confirmed with: https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/?cid=nrcs142p2_054167

setwd("/Users/juliefowler/OneDrive - Colostate/PhD Work/Studies/Burn Piles/Montana Burn Pile Work/Montana-Burn-Piles-Project")


library(ggtern)
data(USDA)
head(USDA, 10)


library(dplyr)
USDA_text <- USDA  %>% group_by(Label) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
USDA_text


soil_data <- read.delim("SoilTexture_forR_TernaryPlot.txt",header = T,check.names=FALSE)
soil_data <- soil_data[1:9,]
summary(soil_data)



Soil_Textural_Triangle <- ggplot(data = USDA, aes(
  y = Clay,
  x = Sand,
  z = Silt
)) +
  coord_tern(L = "x", T = "y", R = "z") +
  geom_polygon(
    aes(fill = Label),
    alpha = 0.0,
    size = 0.5,
    color = "black"
  ) +
  geom_text(data = USDA_text,
            aes(label = Label),
            color = 'black',
            size = 3) +
  theme_showarrows() +
  theme_clockwise() +
  guides(fill=FALSE, color=FALSE) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(size= 10)) +
  theme(axis.title = element_text(size = 15)) +
  theme_nomask() +
  theme(legend.position="right") 


Soil_Textural_Triangle <- Soil_Textural_Triangle + 
  geom_point(data = soil_data, size=3,
  aes(x = Sand,y = Clay,z = Silt, shape = Site)) 


Soil_Textural_Triangle


