# Assignment operators

my_variable_name_1 <- 10
10 -> my_variable_name_2 
my_variable_name_3 = 20


# Data Types -------------------------------------------------------------

#numeric - (10.5, 55, 787)
#integer - (1L, 55L, 100L, where the letter "L" declares this as an integer)
#complex - (9 + 3i, where "i" is the imaginary part)
#character (a.k.a. string) - ("k", "R is exciting", "FALSE", "11.5")
#logical (a.k.a. boolean) - (TRUE or FALSE)

#numeric
var <- 10.5 
class(var)

var <- 10
class(var)

var <- 10L
class(var)

is.numeric(var)
is.integer(var)

var <- 'AT2'
is.numeric(var)
is.character(var)

#Boolen
var <- T
var <- TRUE

var <- "TRUE"

isTRUE(var)
isFALSE(var)


#This will give an error.. sometimes
true <- FALSE
var <- true
class <- 10


# Vectors!!! ---------------------------------------

c(1,2,3,4,5)
1:5
a<-10:15


# how to select or remove elements.  
a[1]
a[6]
a[-2]

i=2
a[i]

#Seq function to create sequences

seq(from=0,to=10,by=2)
seq(1,10,by=2)


rep(c(1:5))
rep(c(1:5),3)
?rep


#Creating vector with random numbers

runif(10,min=1,max=25)

set.seed(100)
runif(10,min=1,max=25)
runif(10,min=1,max=25)
set.seed(100)
runif(10,min=1,max=25)

r50 <-runif(50,min=1,max=25)



# Functions for vectors ---------------------------------------------------

length(r50)
sort(r50)


#Use the summary function to get some details
summary(r50)
mean(r50)
#Other numeric functions for a vector. 
min(r50)
max(r50)
sd(r50)



#Vectors can store numerical, character, logical,dates

c1 <- c('DCM','DCM')
c2 <- c('HCM','HCM','HCM')
c3 <- c(c1,c2)
c3


groups <- c(rep('DCM',10),rep('HCM',15),rep('NF',20),rep('nonfail',5))
groups

#Create a random vector of chars, I'll explain the factor function later. 

groups.a <- sample(c('DCM','HCM','NF'),50,replace=T)
groups.a <- sample(c('DCM','HCM','NF'),4,replace=T)
#Table function works great for a summary
table(groups)
groups.b<-sample(c('DCM','HCM','NF'),50,replace=T, prob = c(.2,.3,.5))
summary(factor(groups.b))



#Create a logical vector         

logical <- c(T,F,F,T,T)
logical

#What if we tried to sum a char vector
sum(groups)

#What happens when we do a logical?
sum(logical)

#We can make logical vectors with a comparison
r50 > 5
r50 >= 5
r50 == 5
#Then 
sum(r50>10)
table(r50>5)

r50[r50>15]

#we can use the table function to get a summary
table(groups)

sum(r50 > mean(r50))

#Can we do the same for chars? yes
groups=="HCM"

sum(groups=='HCM')

#Selecting elements from a char vector. 

groups[2:5]
groups[groups=='HCM' | groups=='NF']

groups=='HCM'


#Using the 'or' and 'and' operator. 
groups[groups=='HCM' | groups=='DCM']
groups[groups=='HCM' & groups=='DCM']

#The 'in' operator 

groups[!groups %in% c('HCM','DCM')]
#### DANGER DO NOT EVER I MNEAN EVER DO THIS!!!!!!!!!!!!!!!!!!
groups[groups == c('HCM','DCM')]

#The not/negation operator 
groups[groups != 'NF']

#Let's use a comparisons to fix the nonfail label

groups[groups=='nonfail']

groups[groups=='nonfail'] <- 'NF'

sum(groups=='nonfail')



#What can't we do,Cannot mix data types, everyting will be cast to char.  
c(1,2,'DCM',F)





# Factors ----------------------------------------------------------------


groups.factor <- factor(groups)
groups.factor

as.numeric(groups)
as.numeric(groups.factor)

plot(groups.factor,r50)

#reorder the factor to have the order we want on the x sxis
groups.factor <- factor(groups,levels = c('NF','HCM','DCM'))
groups.factor
plot(groups.factor,r50)



# NAs vs Nulls ---------------------------------------------------------------------

# NaN  : means 0/0 -- Stands for Not a Number
# NA   : is generally interpreted as a missing, does not exist or undefined
# NULL : is for empty object.it's a reserved word. NULL is perhaps returned by expressions and functions, so that values are undefined.


r50.na <- r50
r50.na[c(3,8,10,12)] <- NA
r50.na

summary(r50.na)
mean(r50.na)

r50.nona <- na.omit(r50.na)
mean(r50.nona)

mean(r50.na,na.rm = T)
