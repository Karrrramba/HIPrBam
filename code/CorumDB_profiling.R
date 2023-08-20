#load libraries
library(tidyverse)
library(reshape2)
library(readxl)


#1 load data 
my_dat <- read_excel(file.choose())
my_dat 

#1.1 load human complex annotation file from CORUM
corum_hsa <- read.delim("humanComplexes.txt")
corum_hsa 

#1.2 OPTIONAL load Functional complex groups annotation file from CORUM
# uncomment to use

corum_fcg <- read.delim("fcg.txt") %>%
  subset(., Organism == "Human")

#load a pre-filtered list of functional complexes of interest
my_fcg  <- subset(corum_fcg, Functional.Complex.Group == "Toll-like receptor signaling")  
my_fcg

#2 split annotated Corum entries in your data file into individual complexes 
#and group by Corum name and condition
cor_set <- my_dat %>%
  mutate(Corum = strsplit(as.character(Corum), ";")) %>%
  unnest(Corum) %>%
  group_by(Corum, Condition) %>%
  summarise_each(funs(mean)) %>%
  .[-3] 

#3 match CORUM terms from own data set with database 
#and retrieve IDs and annotations
data_annot <- merge(x = cor_set, y = corum_hsa[ , c("ComplexName", 
                                                 "ComplexID", 
                                                 "GO.ID", 
                                                 "GO.description",
                                                 "FunCat.ID",
                                                 "FunCat.description",
                                                 "Complex.comment")],
                 by.x = "Corum", by.y = "ComplexName") %>%
  .[,c(38,1:37,39:43)]

#save annotated data set
write.csv(data_annot, "corum_matched_complexes.csv") 

#3.1 OPTIONAL match and annotate user data with functional complex groups
# uncomment to use
# mydat_fcg <- merge(x = cor_set, y = my_fcg[ , c("ComplexName", 
#                                                     "ComplexID", 
#                                                     "GO.ID", 
#                                                     "GO.description",
#                                                     "FunCat.ID",
#                                                     "FunCat.description",
#                                                     "Complex.comment",
#                                                     "category_name")],
#                     by.x = "Corum", by.y = "ComplexName") %>%
#   .[,c(38,1:37,39:43)]

#saved annotated data set
write.csv(data_fcg_annot, "data_fcg_matched.csv") 

#4 tidy up data frame from annotations and 
#transform into long format for plotting
cor_plot <- data_annot %>% 
  .[, -c(2, 39:43)] %>% 
  group_by(ComplexID, Condition) %>%
  pivot_longer(cols=3:37, names_to = "Fraction",  values_to = "RelativeIntensity")

#4.1 tidy up data frame from annotations and 
#transform into long format for plotting
# uncomment to use
#fcg_plot <- mydat_fcg %>% 
  .[, -c(2, 39:44)] %>% 
  group_by(ComplexID, Condition) %>%
  pivot_longer(cols=3:37, names_to = "Fraction",  values_to = "RelativeIntensity")
# in #6 "cor_plot must be changed to "fcg_plot
  
#5 create list of CorumIDs for plotting loop
uniq_comp = unique(cor_plot$ComplexID)
uniq_comp

#6 loop for profile plots of both conditions for each CorumID
for (i in uniq_comp) {
  temp_plot = ggplot(data= subset(cor_plot, ComplexID == i), aes(fct_inorder(Fraction), RelativeIntensity))  +
    geom_point(aes(group=Condition, color=Condition)) +
    geom_line(size=2.0, aes(group=Condition, color = Condition), alpha=0.5) +
    ggtitle(i) +
    scale_x_discrete(breaks = seq(1, 35, 2)) +
    scale_color_manual(values = c("#00FFFF", "#CC0066", "#003333")) +
    labs(x="Fraction",
         y="Relative Intensity [%]") +
    theme(text=element_text(size = 14, family = "Arial"),
          plot.title = element_text(face = "bold",
                                    color = "black",
                                    size = "18",
                                    hjust = .5,
                                    vjust = 2))
  
  ggsave(plot = temp_plot, file=paste("CorumID_", i, ".png", dev = "png",
                                      dpi = 300, width = 14, height = 10, 
                                      units = "cm"))
}



  temp_plot2 = ggplot(data= subset(cor_plot, ComplexID == "1743"), aes(fct_inorder(Fraction), RelativeIntensity))  +
    geom_point(aes(group=Condition, color=Condition)) +
    geom_line(size=1.0, aes(group=Condition, color = Condition), alpha=0.7) +
    ggtitle(i) +
    scale_x_discrete(breaks = seq(1, 35, 2)) +
    scale_color_manual(values = c("#00FFFF", "#CC0066", "#003333")) +
    labs(x="Fraction",
         y="Relative Intensity [%]") +
    theme(text=element_text(size = 14, family = "Arial"),
          plot.title = element_text(face = "bold",
                                    color = "black",
                                    size = "18",
                                    hjust = .5,
                                    vjust = 2))

  ggsave(plot = temp_plot2, file=paste("CorumID_1743.jpeg", dev = "jpeg",
                                      dpi = 300, width = 5, height = 5, 
                                      units = "cm"))
