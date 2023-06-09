## Map for Publication

## Using guide from: https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggspatial)

setwd("/Users/juliefowler/OneDrive - Colostate/PhD Work/Studies/Burn Piles/Montana Burn Pile Work/Spatial Data") 

## Google API key needed here - will have to create one to utilize this code


## Getting the background satellite map 

colors <- c("#fee0d2", "#fc9272", "#de2d26")

sites <- read.delim("Burn_Pile_Locations.txt", header=T)
sites$Site <- sites$Pile.ID
sites$Site <- c("Lonesome Wood Lake", "Lonesome Wood Lake","Lonesome Wood Lake","Lonesome Wood Lake","Lonesome Wood Lake",
                "Lonesome Wood Ridge", "Lonesome Wood Ridge", "Lonesome Wood Ridge", "Lonesome Wood Ridge", "Lonesome Wood Ridge",
                "Rendezvous", "Rendezvous", "Rendezvous", "Rendezvous", "Rendezvous")
sites <- as.data.frame(sites)

sbbox <- make_bbox(lon = sites$Long, lat = sites$Lat, f = .1)
sbbox

# First get the map. By default it gets it from Google. I want it to be a satellite map
sm <- with(sites,c(Long=median(Long),Lat=median(Lat)))
sq_map <- get_map(location = sm, zoom=11,
                   source = "google", maptype = "satellite")
#sq_map <- get_map(location = sbbox, maptype = "satelite", source = "google")
#> Warning: bounding box given to google - spatial extent only approximate.
#> converting bounding box to center/zoom specification. (experimental)
ggmap(sq_map) + geom_point(data = sites, mapping = aes(x = Long, y = Lat, color = Site), size = 4, shape = 17) +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.4)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  theme(legend.title.align=0.5, legend.title = element_text(face="bold", size = 12),
        legend.position = c(0.835, 0.895),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(vjust=-1)) +
  xlab("Longitude") + ylab("Latitude") 



## Sites - Getting points to overlay on satellite map

library("rnaturalearth")
library("rnaturalearthdata")
library("maps")
library("sf")
library("lwgeom")

sf_use_s2(FALSE)
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
counties <- st_as_sf(map("county", plot = FALSE, fill = TRUE))
counties <- subset(counties, grepl("montana", counties$ID))
counties$area <- as.numeric(st_area(counties))

sites <- read.delim("Burn_Pile_Locations.txt", header=T)
sites <- as.data.frame(sites)
sites <- st_as_sf(sites, coords = c("Long", "Lat"), 
                  crs = 4326, agr = "constant")
sites$Site <- sites$Pile.ID
sites$Site <- c("Lonesome Wood Lake", "Lonesome Wood Lake","Lonesome Wood Lake","Lonesome Wood Lake","Lonesome Wood Lake",
                "Lonesome Wood Ridge", "Lonesome Wood Ridge", "Lonesome Wood Ridge", "Lonesome Wood Ridge", "Lonesome Wood Ridge",
                "Rendezvous", "Rendezvous", "Rendezvous", "Rendezvous", "Rendezvous")

ggplot(data = world) +
  geom_sf(fill = "white") +
  geom_sf(data = counties, color = "black", fill = "darkgrey", show.legend = F) +
  geom_sf(data = states, fill = NA) + 
  geom_sf(data = sites, size = 4, shape = 17, aes(color = Site)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-111.48, -111.02), ylim = c(44.55, 44.85), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.4)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(vjust=-1))



## USA - inset map in top right corner
## Guide: https://remiller1450.github.io/s230s19/Intro_maps.html 

MainStates <- map_data("state")
ggplot() + 
  geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),
               color="black", fill="white") +
  xlab("Longitude") + ylab("Latitude") +
  geom_sf(data = sites, size = 4, shape = 17, aes(color = Site)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.4)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  scale_color_manual(values = colors) +
  theme(legend.position="none") 



## Will have to combine these maps in Illustrator or another software








