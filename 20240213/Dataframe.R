#############################################################
#   20240213 MPM R club Data Frame Basics ANd reading Data
#
#
#
############################################################

#Very basic matrix intro
m1 <- matrix(1:10,nrow=5)
m1
#Matrix is indexed like matrix[row,column]
#First row
m1[1,]
#first col
m1[,1]

#a single element
m1[5,2]

#We can add row or col names. 

colnames(m1) <- c('GeneA','GeneB')
rownames(m1) <- paste('sample',1:5,sep='')
m1

m1['sample1',]
m1[,'GeneA']

#Can perform arithmetic on all elements of matrix
m1+10
m1*2


m2 <- m1*2
m1*m2

#get the size of the matrix
dim(m1)


# R has built in function to sum by row or Cols
rowSums(m1)
colSums(m1)

#Can perform comparisons 
m1 > 5
sum(m1>5)

#mix comparisons and rowSum/colSum to calculate % of samples greater than a value

colSums(m1>4)/dim(m1)[1]*100

summary(m1)

##################################################################
#
#  Data Frames
#
#################################################################


# Data Frames. Similar to matrix but each column can be a different data type

id <- paste('PENNID',1:50,sep='')
sex <- sample(c('F','M'),50,replace = TRUE)
smoker <- sample(c(T,F),50,replace = TRUE)
etiology <- c(rep('BOS',10),rep('COPD',15),rep('NF',25))
SFTPC <-runif(50,min=2,max=10)
PDGFRA<-runif(50,min=5,max=20)


df <- data.frame(
                 'sex'=sex,
                 'smoker'=smoker,
                 'etiology'=etiology,
                 'SFTPC'=SFTPC,
                 "PDGFRA"=PDGFRA,
                 row.names = id
                 )

df
#Use summary to get basic stats on the data.frame
summary(df)

df$sex <- factor(df$sex)
df$etiology <- factor(df$etiology)

summary(df)


#We can use table function to summarize a column
table(df$etiology)

#Two columns
table(df$etiology,df$smoker)

#Three columns!!
table(df$etiology,df$smoker,df$sex)

#What about numeric columns? 
table(df$etiology,df$SFTPC)

#Multiple ways to access cols, rows and "cells". 

#Columns just like matrix

df[,1]
df[,2:3]

#or by column name 

df$SFTPC
df[,'SFTPC']

#More that one column 
df[,c('SFTPC','PDGFRA')]


#We can select rows by names. 

df[c('PENNID39','PENNID39'),]

df[c('PENNID39','PENNID39'),'SFTPC']

#What if forget the ','? Well it selects the columns 

df[c('SFTPC','PDGFRA')]

#But if we give it columns names but mistakenly put the , to reference rows we get something strange. 
df[c('SFTPC','PDGFRA'),]

#How to subset a data frame 

#let's build first with a comparison, these returns a logical vector
#for each row. 

df$smoker == F

#We can then use that to select the rows we want. 

df[df$smoker==F,]

#What about multiple conditions, say all Females with COPD

df[df$etiology=='COPD' & df$sex=='F',]

## Don't forget if you want to match multiple values use %in%

df[df$etiology %in% c('COPD','BOS') & df$sex=='F',]


#Sometimes it's helpful to know what data structure we have. 
is.data.frame(df)
is.data.frame(df$SFTPC)
is.vector(df$SFTPC)


####################
#   read in data using base R
#######################
d.base <- read.csv('testData.csv')
d.base
summary(d.base)

#Let's start using the Tidyverse
library(tidyverse)
library(lubridate) 

d.tidy <- read_csv('testData.csv')

d.tidy

#Let's fix the date, we use lubridate function parse_date_time

d.tidy$DOB <-  parse_date_time(d.tidy$DOB, orders = c('mdy', 'dmy'))



####################################################
#   Write out the results
#
####################################################

#Using Base R 
write.csv(x=df,file='TestDf_base.csv')
#Using Tidy R 
write_csv(x=df,path='TestDf_tidy.csv')



#####################################################
# Tidyverse Preview
#
######################################################

#Column Selction in base R
df[,c('SFTPC','PDGFRA')]

#Tidyverse (dplyr)
select(df,c(SFTPC,PDGFRA))

#filter/Subset in Base R

df[df$etiology=='COPD' & df$sex=='F',]

#Tidyverse
filter(df,etiology=='COPD' & sex=='F')

