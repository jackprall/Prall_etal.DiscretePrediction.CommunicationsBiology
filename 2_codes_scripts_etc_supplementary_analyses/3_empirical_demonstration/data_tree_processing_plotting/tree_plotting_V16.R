# plotting Dodo data
# Prall et al. 2025
# Data Storchova & Horak 2017
# topology from timetree.org Feb 2025
# labeling major bird orders

library(ggplot2)
library(ggtree)
library(geiger)
library(ggnewscale)
library(export)

tree<- read.nexus("./ChrisBirds/BirdTreeCleaned.tre")
tree$tip.label<- gsub("_", " ", tree$tip.label)

dat<- read.table("./ChrisBirds/DataCleaned.txt")
dat$V1<- gsub("_", " ", dat$V1)

dat[,2]<- as.factor(dat[,2])
rownames(dat)<- dat[,1]
d<- ggimage::phylopic_uid(c("Hieraaetus", "Raphus cucullatus"))
dat$V5_tippoint<- dat$V3
dat$V5_tippoint[dat$V1 %in% c("Apus pallidus", 
                              "Buteo rufinus",
                              "Lanius minor", 
                              "Passer domesticus", 
                              "Sylvia mystacea")]<- "Prediction Accuracy Tests"

selects<- which(tree$tip.label %in% c("Apus pallidus", 
                                      "Buteo rufinus",
                                      "Lanius minor", 
                                      "Passer domesticus", 
                                      "Sylvia mystacea"))

dat$V5_tippoint[dat$V1 == "Raphus cucullatus"]<- "Raphus cucullatus"
dat$V5_tippoint[dat$V1 == "Hieraaetus moorei"]<- "Hieraaetus moorei"

desired_order <- c("0", "1", "Prediction Accuracy Tests", "Hieraaetus moorei", "Raphus cucullatus")
dat$V5_tippoint <- factor(dat$V5_tippoint, levels = desired_order)

d$uid[1]<- "e038a748-825f-4cdc-8909-afe6c29c26fc"
d$uid[2]<- "01bbffec-312c-493a-bcd3-e9361ed3435d"


theme_set(theme(text = element_text(family = "sans")))
treeplot <- ggtree(tree, color = "grey60",
                   layout = "circular",
                   barsize = 2) %<+% dat

treeplot_2 <- treeplot +
  geom_tiplab2(size = 1,
               offset = 7,
               aes(color = V5_tippoint),
               show.legend = FALSE) +
  scale_color_manual(values = c("grey30","grey30","dodgerblue4",
                                "deeppink4","chartreuse4")) +
  new_scale_color() +
  geom_tippoint(aes(color = V2), size = 0.75) +
  scale_color_manual(values = c("pink","orchid4","dodgerblue2",
                                "chartreuse3","deeppink2")) +
  geom_tippoint(aes(color = V5_tippoint),
                position = position_nudge(x = 4),
                size = 0.75) +
  
  geom_cladelabel(node=498, label="Columbiformes",
                  color="#C45C8A", offset=30,
                  angle=30, offset.text=3,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=430, label="Gruiformes",
                  color="#CE4441FF", offset=30,
                  offset.text=3, angle=30,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=410, label="Pelicaniformes",
                  color="#EE8577FF", offset=30,
                  angle=30, offset.text=3,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=436, label="Charadriiformes",
                  color="#EB7930FF", offset=30,
                  offset.text=30,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=516, label="Afroaves",
                  color="orange1", offset=30,
                  offset.text=30,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=530, label="Accipitriformes",
                  color="#6D8255", offset=30,
                  offset.text=49,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=555, label="Passeriformes",
                  color="#62929AFF", offset=30,
                  offset.text=4,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=751, label="Galliformes",
                  color="#575AAB", offset=30,
                  offset.text=5,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=720, label="Anseriformes",
                  color="#7D57AB", offset=30,
                  offset.text=3,
                  barsize=1.25, fontface="bold") +
  
  geom_cladelabel(node=402, label="Apodiformes",
                  color="#B550AB", offset=30,
                  offset.text=3,
                  barsize=1.25, fontface="bold") +
  
  geom_tippoint(aes(subset=node==104),
                position=position_nudge(x=4),
                color="deeppink2", size=2) +
  
  geom_tippoint(aes(subset=node==153),
                position=position_nudge(x=4),
                color="chartreuse3", size=2) +
  
  geom_tippoint(aes(subset=node %in% selects),
                position=position_nudge(x=4),
                color="dodgerblue2", size=2) +
  
  labs(color="Trait Values:\nNest Type (inner)\nNesting Behavior (outer)") +
  
  theme(legend.position=c(0.2,0.9),
        legend.text=element_text(size=10)) +
  
  guides(colour=guide_legend(override.aes=list(size=3)))
  
  

treeplot_2

ggsave(filename="~/Desktop/birdtree_plot3.pdf", width=10.5, height=10.5, units= "in", dpi=500)
