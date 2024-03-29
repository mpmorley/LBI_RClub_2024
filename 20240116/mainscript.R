###################################################################################################
# 1/16/2024 Created by MPM fro R Club
#
#################################################################################################






# install non-CRAN packages -----------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBImage")

#From github, we need a special library/package to do so. 
#Can use Rstudio tools->Install packages or this command
install.packages("devtools") 


#Load the library devtools
library(devtools)

install_github("BlakeRMills/MetBrewer")

# or we do need to need to load the library into memmory. 

devtools::install_github("BlakeRMills/MetBrewer")



# Let's some code!! -------------------------------------------------------

#Let's setup the structure of a project directory using some R code, 
#yes we code this manually but what fun is that. 

# Code creates directoriies for project organization
#  bin => Keep R scripts, Can rename to scripts, code etc
#  data => raw data files, meta data etc. Typically data not generated by us
#  results => Analysis results, csv, excel, etc. 
#  plots => Where plots. 
#  misc =>

# Brute force method ------------------------------------------------------
dir.create('bin', showWarnings = )
dir.create('data',showWarnings = T)
dir.create('plots',showWarnings = T)
dir.create('results',showWarnings = F)
dir.create('misc',showWarnings = F)

# Smart method "For" loop -----------------------------------------------

#Create a vector of dirs we wish to make
dirs <- c('bin','data','plots','results','misc')

#We can select an element in a vector by it's index, which is it's position starting from 1 


#Paste concatenates multiple arguments into a single string
#By default puts a space between elements
paste('A','B')
#paste0 is a shortcut for no space. 
paste0('A','B')
paste('A','B', sep='')


#long form of a for loop, where we use an index to select each element 
for(i in 1:length(dirs)) {
  print(paste0('Creating ', dirs[i], ' Directory'))
  dir.create(dirs[i], showWarnings = F)
}

## Shorter version where each element is stored in the dir. I tend to name to vector with a plural noun and the element as a singular noun. 
for(dir in dirs) {
  print(paste0('Creating ', dir, ' Directory'))
  dir.create(dir, showWarnings = F)
}

# Smarter method tidyverse way using a custom function and Map ---------------------

#Create a function that makes a directory
mkdir <- function(dir){
  print(paste0('Creating ',dir,' Directory'))
  dir.create(dir,showWarnings = F) 
}
mkdir('bin')

map(dirs,mkdir)


# Smartest save code in separate file and call when we need --------

source('script2.R')




