library(tidyverse)
library(lubridate)


df <- read_csv('testData.csv')


# Select: Subset variable (columns) ------------------------------------------------------------------
#single column
select(df,smoker)
#Range of cols
select(df,id:etiology)
#remove a col
select(df,-DOB)

select(df,starts_with('s'))



# filter: Subset observations (rows) ------------------------------------------------------------------

filter(df,sex=='F')

filter(df,sex=='F' & smoker==TRUE)



# Pipes -------------------------------------------------------------------
# %>% is called a pipe 
#Pipes send the output of one function as the first argument to the next function. 

df %>% select(sex,smoker) %>%
  filter(sex=='F')


# Mutate ------------------------------------------------------------------
df <- df %>% 
  mutate(DOB=parse_date_time(DOB, orders = c('mdy', 'dmy'))
  )

df %>%
  mutate(age = today()-as.Date(DOB))

#Let's calculate BMI, first we need to convert inches to meters. 

df <- df %>% 
  mutate(ht_m = ht_in*0.0254
  )

df <- df %>% 
  mutate(BMI = wt_kg/(ht_m^2)
  )


df %>% 
  filter(DOB>=as.Date("1980-01-01"))

# ArrangeS Sort column by observations ------------------------------------------------------------------

df %>% arrange(-wt_kg)
df %>% arrange(wt_kg)
df %>% arrange(etiology,-wt_kg) %>% head(n=5)


# summarize --------------------------------------------------------------

df %>% summarise(mean(ht_in),sd(ht_in))

df %>% summarise(mean_ht=mean(ht_in),sd_ht=sd(ht_in))


# group_by ----------------------------------------------------------------

df %>% group_by(etiology) %>%
  summarise(mean(ht_in),sd(ht_in))

df %>% group_by(etiology,smoker) %>%
  summarise(mean(ht_in),sd(ht_in))

df %>% group_by(sex,etiology,smoker) %>%
  summarise(mean(ht_in),sd(ht_in))





# join --------------------------------------------------------------------

#Open a 2nd file and "merge' results. 

genes <- read_csv('geneData.csv')

dim(genes)
dim(df)



df %>% filter(id %in% genes$sampid)
#We can use the "not" in front and retirn the not matching
df %>% filter(!id %in% genes$sampid)

#using join, there's multiple flavors of join.. We'll look at inner_join and left/right join 


merge.inner <- inner_join(df,genes,by=c('id'='sampid'))
dim(merge.inner)
merge.left <- left_join(df,genes,by=c('id'='sampid'))
dim(merge.left)


# long data ---------------------------------------------------------------
# Part of tidyr package, use to be called gather and spread

#Go from wide to long, the cols defines which varaibles to pivot,in our case the 
#genes, we know have 2 variables, one for genes and one for expression. 
#We know have a dataset that is no longer one observation per sample but rather
#an observation for each gene per sample. 
tmp <- pivot_longer(merge.inner,
                    cols=Gene1:Gene10
)


# Same as before but we give names to the 2 new labels to the variables we created.
merge.long <- pivot_longer(merge.inner,
                           cols=Gene1:Gene10,
                           names_to='gene',
                           values_to='UMI'
)

#now we can summarize the all the genes quite easily. 

summ <- merge.long %>% 
  group_by(etiology,gene) %>%
  summarise(mean=mean(UMI),
            sd=sd(UMI),
            min=min(UMI),
            max=max(UMI)
  )



## Slice family  function to subset rows by "position"

summ %>% slice(1)

summ %>% slice_max(mean,n=1)
summ %>% slice_max(mean,n=2)

#Notice the data is still grouped.. We can ungroup it as well.  

summ %>% ungroup %>% slice_max(mean,n=1)


#We can pivot the summary table to a more human readable version. 

pivot_wider(summ, names_from = c(etiology), values_from = c(mean,sd,max,min))


#We can plot long data nicely as well. 

ggplot(merge.long,aes(x=etiology,y=UMI,fill=etiology)) + geom_boxplot() + facet_grid(~gene) + theme_bw()

#We can even pipe to ggplot as well. 

pivot_longer(merge.inner,
             Gene1:Gene10,
             names_to='gene',
             values_to='UMI') %>%
  ggplot(aes(x=etiology,y=UMI,fill=etiology)) + geom_boxplot() + facet_grid(~gene) + theme_bw()

# Tidyverse fct_reorder function. 


pivot_longer(merge.inner,
             Gene1:Gene10,
             names_to='gene',
             values_to='UMI') %>%
  mutate(gene=fct_reorder(gene,UMI),gene) %>%
  ggplot(aes(x=etiology,y=UMI,fill=etiology)) + geom_boxplot() + facet_grid(~gene) + theme_bw()






