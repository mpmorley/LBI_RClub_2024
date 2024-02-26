library(tidyverse)
library(patchwork)
#These are various color palette packages
install.packages("MetBrewer")
install.packages("ggsci")
install.packages("wesanderson")
library(wesanderson)
library(ggsci)
library(MetBrewer)




#Lists. Like vectors but more flexible. An element can be anything. 
# Text, numeric, data frame or another list. 

df <- data.frame(colA=1:10,colB=rep('A',10))

#Create a list 

a.list <- list("Blue", 10, c(21,32,11), TRUE,df)
a.list

# Using a single [] we create a list slice
a.list[1]
#or
a.list[c(1,4)]

# Using double [[]] is the reference, this is the same as single [] in a vector
a.list[[1]]



#Lists can have each element with a name. 

b.list <- list(color="Blue", id=10, summary=c(21,32,11), afib=TRUE,data=df)
b.list

b.list['color']
b.list[['color']]


#Retrieve the names of the list
names(b.list)

## Add to out list

b.list$age <- 40
b.list[['sex']] <- 'F'

names(b.list)


#Convert list to a vector
unlist(b.list)


#Loop over lists
d.list <- list('purple','green','red')

#Using a for loop, this works only for list of 3 elements 
for(i in 1:3){
  print(d.list[[i]])
}

#We can modify it so it takes any size list
for(i in 1:length(d.list)){
  print(d.list[[i]])
}

#But what if we wanted to save the results to a new list, How ugly! 
newlist <- list()
for(i in 1:length(d.list)){
  newlist[[i]] <- print(d.list[[i]])
}

newlist

# Loop over list the R way 

newlist <- lapply(d.list,print)
newlist

## Or the tidyverse way, the map family of functions are in the purrr library

newlist <- map(d.list,print)
newlist



## Now with numbers 
c.list <- list(c(runif(10,min=1,max=25)),runif(10,min=1,max=25),runif(10,min=1,max=25))
c.list
#using base R
lapply(c.list,mean)
#Using purrr library from tidyverse
map(c.list,mean)


# ggplot! -----------------------------------------------------------------



#read in data, make sure you alter this to your path
data <- read_csv('scMetaTest.csv') 

# We start with the ggplot function to define plot data and aesthetics. The we add layers, points. lines etc/. 

ggplot(data,aes(x=celltype,y=percent.mito)) + geom_point()

#Let's jitter the points, so they do not lie on top of each other
ggplot(data,aes(x=celltype,y=percent.mito))  + geom_jitter(position=position_jitter(0.2))



# Boxplots ----------------------------------------------------------------
#Boxplot
ggplot(data,aes(x=celltype,y=percent.mito))  + geom_boxplot()
#Points and Boxplots, much better
ggplot(data,aes(x=celltype,y=percent.mito)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2))

#using fill 
ggplot(data,aes(x=celltype,y=percent.mito,fill=Gender))  + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2))
#using color
ggplot(data,aes(x=celltype,y=percent.mito,color=Gender))  + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2))
#using both, but we have an issue with jitter, the points are mixed
ggplot(data,aes(x=celltype,y=percent.mito,color=Gender,fill=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2))
#We can dodge the points so the appear with he proper box, but not jittered
ggplot(data,aes(x=celltype,y=percent.mito,color=Gender,fill=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_dodge(0.2))

#Putting it all together with jitterdodge
ggplot(data,aes(x=celltype,y=percent.mito,color=Gender,fill=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2))


#Let's make our own theme. 

myTheme <- theme_bw() + theme(
  legend.position='bottom',
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black")
)



ggplot(data,aes(x=celltype,y=percent.mito,color=Gender,fill=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2)) +
  myTheme
  


#Use other color pals, this is ggsci 
ggplot(data,aes(x=celltype,y=percent.mito,color=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2)) +
  myTheme +
  scale_color_aaas()

#Wes anderson
ggplot(data,aes(x=celltype,y=percent.mito,color=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2)) +
  myTheme + 
  scale_color_manual(values = wes_palette("Royal1"))

#Metcolorbrewer 
ggplot(data,aes(x=celltype,y=percent.mito,color=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2)) +
  myTheme + 
  scale_color_manual(values=met.brewer("Juarez", 2))



# Let's put multiple plots together with patchwork -------------------------

#set plots to a variable, these will be ggplot objects

p1 <- ggplot(data,aes(x=celltype,y=percent.mito,color=Gender)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2))å

p1

p2 <- ggplot(data,aes(x=nCount_RNA,y=percent.mito)) + geom_point()

#We can add layers to the ggplot2 object as well. Here we formating the scale to display % using the scales package
p2 <- p2 + scale_y_continuous(labels = scales::percent_format(scale = 1))


#Use patchwork operators to put plots togather. 

p1+p2
p1/p2 

#We can add layers to patchwork object too using '*'
p1/p2 * myTheme



#### let's make a figure, back to lists!

plots.list <- list()
plots.list$A <- ggplot(data,aes(x=celltype,y=percent.mito,color=ChemistryVer)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2))
plots.list$B <- ggplot(data,aes(x=celltype,y=hybrid_score,color=ChemistryVer)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitterdodge(0.2))

plots.list$C <- ggplot(data,aes(x=nCount_RNA,y=percent.mito)) + geom_point()



wrap_plots(plots.list) & myTheme & scale_color_aaas() 


wrap_plots(plots.list,guides='collect') & myTheme & scale_color_aaas() 


#We can know apply theme each plot. 

#Add some annotation 
wrap_plots(plots.list,guides='collect') + plot_annotation(tag_levels = 'A')& myTheme & scale_color_aaas() 



#With patchwork we layouts using the "+","/". 

plots.list$A + (plots.list$B/plots.list$C) + 
  plot_annotation(tag_levels = 'A') & myTheme & scale_color_aaas() 

plots.list$A + (plots.list$B/plots.list$C) + 
  plot_annotation(tag_levels = 'A') + plot_layout(guides='collect') & myTheme & scale_color_aaas() 







