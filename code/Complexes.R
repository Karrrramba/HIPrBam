# delete "#" to install packages 
#install.packages("tidyverse")
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("readxl")

#load libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(readxl)


#load data 
data_comp <- read_excel(file.choose())
data_comp 

#split annotated Corum entries into individual complexes and group elution profiles of each condition by gene name
cor_profile <- complex_dat %>%
  mutate(Corum = strsplit(as.character(Corum), ";")) %>%
  unnest(Corum) %>%
  group_by(Corum, Condition) %>%
  summarise_each(funs(mean)) %>%
  .[-3] %>%
  pivot_longer(cols=3:37, names_to = "Fraction",  values_to = "RelativeIntensity")

#create list of profiles grouped by complex name
uniq_comp = unique(cor_profile$Corum)
uniq_comp

#export list of Corum terms
capture.output(uniq_comp, file = "CorumIDs.txt")

#read Corum ID list and use it as template for plot name generation
IDlist <- read.table("CorumIDs.txt")


#loop for profile plots of both conditions for each protein
for (i in uniq_comp) {
  temp_plot = ggplot(data= subset(cor_profile, Corum == i), aes(fct_inorder(Fraction), RelativeIntensity))  +
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
  
  
  ggsave(temp_plot, next_file, width = 14, height = 10, units = "cm")
}

#function for sequential file numbering 

next_file = function(basename = 'CorumID_', fileext = 'png', filepath = '~/plots'){
  old.fnames = grep(paste0(basename,' \\d+\\.', fileext,'$'), 
                    list.files(filepath), value = T)
  lastnum = gsub(paste0(basename,' (\\d+)\\.', fileext,'$'), '\\1', old.fnames)
  if (!length(lastnum)) { 
    lastnum = 1 
  } else {
    lastnum = sort(as.integer(lastnum),T)[1] + 1L 
  }
  return(paste0(basename, ' ', sprintf('%03i', lastnum), '.', fileext))
}

ggsave(next_file())




match
row.names(match("(E.F.G) complex", IDlist$V2))
IDlist
row.names(IDlist$V2 == "(E.F.G) complex")

#meanprofile plot
temp_plotmean = data_gr %>% filter(Condition == "Ctrl" | Condition == "IBR") %>%
  ggplot(.) +
  stat_summary(aes(fct_inorder(Fraction), RelativeIntensity, group=Condition, color=Condition), 
               fun = mean, geom = "line", size = 1, alpha) +
  ggtitle("Cytokine") +
  labs(x="Fraction",
       y="Relative Intensity [%]") +
  scale_x_discrete(breaks = seq(1, 35, 2)) +
  scale_color_manual(values = c("#00FFFF", "#CC0066", "#003333")) +
  theme(text = element_text(size = 14, family = "Arial"),
        plot.title = element_text(face = "bold",
                                  color = "black",
                                  size = "18",
                                  hjust = .5,
                                  vjust = 2))

temp_plotmean



#loop for profile plots of both conditions for each protein
for (i in uniq_comp) {
  temp <- ggplot(data= subset(cor_plot, ComplexID == i), aes(fct_inorder(Fraction), RelativeIntensity))  +
    geom_point(aes(group=Condition, color=Condition)) +
    geom_line(size=2, aes(group=Condition, color = Condition), alpha=0.4) +
    ggtitle(i) +
    scale_x_discrete(breaks = seq(1, 35, 2)) +
    scale_color_manual(values = c("#00FFFF", "#CC0066", "#003333")) +
    labs(x="Fraction",
         y="Relative Intensity [%]") +
    theme_bw(base_size = 12, base_family = "Arial", ) +
    theme(plot.title = element_text(face = "bold",
                                    color = "black",
                                    size = rel(1),
                                    hjust = .5,
                                    vjust = 2))
  print(temp)
}


#loop for saving plots

filename <- paste0('CorumID_',cor_plot$ComplexID, '.png') 
path <- '~/Corum'
for(i in temp){
  ggsave(filename = filename, plot = temp[[i]], dpi = 300, device = "png",
         height = 14, width = 10, units = "cm")
}



