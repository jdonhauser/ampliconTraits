# make a trait specific database ****************************************************************************** ####
# for each trait including only the sequences where the value of the trait is not NA
# create table of sequences to include for each trait
# for continuous traits: bin into intervals and define trait value as interval

library(ggplot2)
## setup for graphics ***************************************************************************************************************************** ####
plot.theme1 <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "white",
                                                     colour = "black",
                                                     size = 0.5, linetype = "solid"),
                     panel.border= element_rect(fill=NA,size = 0.5, linetype = 'solid',colour = "black"),
                     axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),legend.text = element_text(size=13),
                     axis.title = element_text(size=14),
                     legend.title = element_text(color = "black", size = 14),
                     strip.text.x = element_text(size=14),
                     strip.background = element_rect(colour="black", fill="white")
)

## Import file *************************************************************************************************************************************************** ####
# Import table matching silva accessions

dbAll <- read.csv('../2_matchTaxids/dbMadin_merged_Silva.csv', header = T) 

# for categorical traits: replace spaces with '_' to retrieve complete value when selecting fields 
for(i in c("gram_stain","metabolism","pathways","carbon_substrates", "sporulation", "motility", "range_tmp", "range_salinity", "cell_shape", "isolation_source")){
  dbAll[,i] <- gsub(' ','_',dbAll[,i])
}

## Categorical Traits ******************************************************************************************************************************************** ####
### Metabolism *************************************************************************************************************************************************** ####
# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$metabolism), c("Feature.ID", "metabolism")]

# export
dir.create('./metabolism')
write.table(db_single,'./metabolism/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

### gram_stain *************************************************************************************************************************************************** ####
# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$gram_stain), c("Feature.ID", "gram_stain" )]

# export
dir.create('./gram_stain')
write.table(db_single,'./gram_stain/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

### sporulation *************************************************************************************************************************************************** ####
# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$sporulation), c("Feature.ID", "sporulation")]

# export
dir.create('./sporulation')
write.table(db_single,'./sporulation/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

### motility *************************************************************************************************************************************************** ####
# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$motility), c("Feature.ID", "motility")]

# export
dir.create('./motility')
write.table(db_single,'./motility/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

### range_tmp *************************************************************************************************************************************************** ####
# values: "thermophilic" "mesophilic" "psychrophilic" "extreme thermophilic" "facultative psychrophilic" "psychrotolerant" "thermotolerant"           
# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$range_tmp), c("Feature.ID", "range_tmp")]

# export
dir.create('./range_tmp')
write.table(db_single,'./range_tmp/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

### range_salinity *************************************************************************************************************************************************** ####
# values_"moderate-halophilic" "extreme-halophilic"  "non-halophilic"      "halophilic"          "stenohaline"         "halotolerant"        "euryhaline" 
# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$range_salinity), c("Feature.ID", "range_salinity")]

# export
dir.create('./range_salinity')
write.table(db_single,'./range_salinity/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

### cell_shape *************************************************************************************************************************************************** ####
# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$cell_shape), c("Feature.ID", "cell_shape")]

# export
dir.create('./cell_shape')
write.table(db_single,'./cell_shape/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

## Continuous traits ********************************************************************************************************************************************* ####
# bin in intervals of different size with cut
# cut decreases/increases the lower/upper limits of the lowest/highest interval by 0.1% to ensure that the borders 
# are within the interval. 

# make a table with original range (what interval should be) and what the interval becomes to convert back later
Lim <- data.frame()

### Topt ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$optimum_tmp), c("Feature.ID", "optimum_tmp")]

# plot histogram
ggplot(db_single, aes(x=optimum_tmp)) + 
  geom_histogram(binwidth=1)
range(db_single$optimum_tmp)

# bin in 10 intervals (add column with range)
db_single$optimum_tmp10 <- cut(db_single$optimum_tmp, breaks = 10) # if break is one number it defines the numbers of bins (which are of equal size)
db_single$optimum_tmp10 <- gsub('\\(','',db_single$optimum_tmp10)
db_single$optimum_tmp10 <- gsub('\\]','',db_single$optimum_tmp10)
db_single$optimum_tmp10 <- gsub(',','-',db_single$optimum_tmp10)

# add original and enlarged interval to conversion table
unique(db_single$optimum_tmp10)
Lim["optimum_tmp10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("3-13.2","2.9-13.2","94.8-105","94.8-105")

a <- db_single[, c("Feature.ID", "optimum_tmp10")]

# export
dir.create('./optimum_tmp10')
write.table(a,'./optimum_tmp10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals
db_single$optimum_tmp10 <- factor(db_single$optimum_tmp10, 
                                  levels = sort(unique(db_single$optimum_tmp10))[c(2,1,3:length(unique(db_single$optimum_tmp10)))])

pdf('./optimum_tmp10/numberOfSequencesPerInterval.pdf')
ggplot(db_single, aes(x=optimum_tmp10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# bin in 20 intervals
db_single$optimum_tmp20 <- cut(db_single$optimum_tmp, breaks = 20) 
db_single$optimum_tmp20 <- gsub('\\(','',db_single$optimum_tmp20)
db_single$optimum_tmp20 <- gsub('\\]','',db_single$optimum_tmp20)
db_single$optimum_tmp20 <- gsub(',','-',db_single$optimum_tmp20)

# add original and enlarged interval to conversion table
unique(db_single$optimum_tmp20)
Lim["optimum_tmp20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("3-8.1","2.9-8.1","99.9-105","99.9-105")

a <- db_single[, c("Feature.ID", "optimum_tmp20")]

# export
dir.create('./optimum_tmp20')
write.table(a,'./optimum_tmp20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals
db_single$optimum_tmp20 <- factor(db_single$optimum_tmp20, 
                                  levels = sort(unique(db_single$optimum_tmp20))[c(3,16,1,2,4:15,17:20)])

pdf('./optimum_tmp20/numberOfSequencesPerInterval.pdf', width = 10)
ggplot(db_single, aes(x=optimum_tmp20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


### d1_lo ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$d1_lo), c("Feature.ID", "d1_lo")]

# plot histogram
ggplot(db_single, aes(x=d1_lo)) + 
  geom_histogram(binwidth=0.2)
range(db_single$d1_lo)
# histogram with values < 750
range(db_single[db_single$d1_lo<750,'d1_lo'])
ggplot(db_single[db_single$d1_lo<750,], aes(x=d1_lo)) + 
  geom_histogram(binwidth=0.2)

# number of sizes >= 80
length(which(db_single$d1_lo>79)) #6

# histogram with values < 80
range(db_single[db_single$d1_lo<80,'d1_lo'])
ggplot(db_single[db_single$d1_lo<80,], aes(x=d1_lo)) + 
  geom_histogram(binwidth=0.2)

ggplot(db_single[db_single$d1_lo<5,], aes(x=d1_lo)) + 
  geom_histogram(binwidth=0.1)

# there are empty bins because there are gaps in range when split in equal intervals
# split such that there is the desired number of bins that are not empty, use smallest amount of bins possible

# bin in 5 non-empty intervals (add column with range)
db_single$d1_lo5 <- cut(db_single$d1_lo, breaks =51)
db_single$d1_lo5 <- gsub('\\(','',db_single$d1_lo5)
db_single$d1_lo5 <- gsub('\\]','',db_single$d1_lo5)
db_single$d1_lo5 <- gsub(',','-',db_single$d1_lo5)

# add original and enlarged interval to conversion table
unique(db_single$d1_lo5)
Lim["d1_lo5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.1-14.8","-0.65-14.8","735-750","735-751")

a <- db_single[, c("Feature.ID", "d1_lo5")]

# export
dir.create('./d1_lo5')
write.table(a,'./d1_lo5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals
db_single$d1_lo5 <- factor(db_single$d1_lo5, 
                                  levels = sort(unique(db_single$d1_lo5)))

pdf('./d1_lo5/numberOfSequencesPerInterval.pdf')
ggplot(db_single, aes(x=d1_lo5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# bin in 10 non-empty intervals (add column with range)
db_single$d1_lo10 <- cut(db_single$d1_lo, breaks =382) 
db_single$d1_lo10 <- gsub('\\(','',db_single$d1_lo10)
db_single$d1_lo10 <- gsub('\\]','',db_single$d1_lo10)
db_single$d1_lo10 <- gsub(',','-',db_single$d1_lo10)

# add original and enlarged interval to conversion table
unique(db_single$d1_lo10)
Lim["d1_lo10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.1-2.06","-0.65-2.06","748-750","748-751")

a <- db_single[, c("Feature.ID", "d1_lo10")]

# export
dir.create('./d1_lo10')
write.table(a,'./d1_lo10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals
db_single$d1_lo10 <- factor(db_single$d1_lo10, 
                           levels = sort(unique(db_single$d1_lo10))[c(1,3,5:7,10,2,4,9,8)])

pdf('./d1_lo10/numberOfSequencesPerInterval.pdf')
ggplot(db_single, aes(x=d1_lo10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# bin in 20 non-empty intervals (add column with range)
db_single$d1_lo20 <- cut(db_single$d1_lo, breaks =1923)
db_single$d1_lo20 <- gsub('\\(','',db_single$d1_lo20)
db_single$d1_lo20 <- gsub('\\]','',db_single$d1_lo20)
db_single$d1_lo20 <- gsub(',','-',db_single$d1_lo20)

# add original and enlarged interval to conversion table
unique(db_single$d1_lo20)
Lim["d1_lo20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.1-0.49","-0.6499-0.49","749.6-750","749.6-750.7")

a <- db_single[, c("Feature.ID", "d1_lo20")]

# export
dir.create('./d1_lo20')
write.table(a,'./d1_lo20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals
# upper border of interval to sort 
int <- unique(db_single$d1_lo20)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d1_lo20)

db_single$d1_lo20 <- factor(db_single$d1_lo20, 
                            levels = names(int)[order(int)])
pdf('./d1_lo20/numberOfSequencesPerInterval.pdf', width = 10)
ggplot(db_single, aes(x=d1_lo20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 30 non-empty intervals (add column with range)
db_single$d1_lo30 <- cut(db_single$d1_lo, breaks =4233)
db_single$d1_lo30 <- gsub('\\(','',db_single$d1_lo30)
db_single$d1_lo30 <- gsub('\\]','',db_single$d1_lo30)
db_single$d1_lo30 <- gsub(',','-',db_single$d1_lo30)

# add original and enlarged interval to conversion table
unique(db_single$d1_lo30)
Lim["d1_lo30",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.1-0.2772","-0.6499-0.2772","749.8-750","749.8-750.7")

a <- db_single[, c("Feature.ID", "d1_lo30")]

# export
dir.create('./d1_lo30')
write.table(a,'./d1_lo30/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals
# upper border of interal to sort 
int <- unique(db_single$d1_lo30)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d1_lo30)

db_single$d1_lo30 <- factor(db_single$d1_lo30, 
                            levels = names(int)[order(int)])
pdf('./d1_lo30/numberOfSequencesPerInterval.pdf', width = 15)
ggplot(db_single, aes(x=d1_lo30)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()



### d1_up ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$d1_up), c("Feature.ID", "d1_up")]

# plot histogram
ggplot(db_single, aes(x=d1_up)) + 
  geom_histogram(binwidth=0.2)
range(db_single$d1_up)

# bin in 5  intervals (add column with range)
db_single$d1_up5 <- cut(db_single$d1_up, breaks =5) 
db_single$d1_up5 <- gsub('\\(','',db_single$d1_up5)
db_single$d1_up5 <- gsub('\\]','',db_single$d1_up5)
db_single$d1_up5 <- gsub(',','-',db_single$d1_up5)

# add original and enlarged interval to conversion table
unique(db_single$d1_up5)
Lim["d1_up5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.01-2.01","1e-05-2.01","8-10","8-10")

a <- db_single[, c("Feature.ID", "d1_up5")]

# export
dir.create('./d1_up5')
write.table(a,'./d1_up5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals
# upper border of interval to sort 
int <- unique(db_single$d1_up5)
int <- as.numeric(ifelse(substr(int,2,2)=='e',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d1_up5)

db_single$d1_up5 <- factor(db_single$d1_up5, 
                            levels = names(int)[order(int)])
pdf('./d1_up5/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=d1_up5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()



# bin in 10 non-empty intervals (add column with range)
db_single$d1_up10 <- cut(db_single$d1_up, breaks =10) 
db_single$d1_up10 <- gsub('\\(','',db_single$d1_up10)
db_single$d1_up10 <- gsub('\\]','',db_single$d1_up10)
db_single$d1_up10 <- gsub(',','-',db_single$d1_up10)

# add original and enlarged interval to conversion table
unique(db_single$d1_up10)
Lim["d1_up10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.01-1.01","1e-05-1.01","9-10","9-10")

a <- db_single[, c("Feature.ID", "d1_up10")]

# export
dir.create('./d1_up10')
write.table(a,'./d1_up10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals (by upper border of interval)
int <- unique(db_single$d1_up10)
int <- as.numeric(ifelse(substr(int,2,2)=='e',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d1_up10)

db_single$d1_up10 <- factor(db_single$d1_up10, 
                           levels = names(int)[order(int)])
pdf('./d1_up10/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=d1_up10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 20 non-empty intervals (add column with range)
db_single$d1_up20 <- cut(db_single$d1_up, breaks =22) 
db_single$d1_up20 <- gsub('\\(','',db_single$d1_up20)
db_single$d1_up20 <- gsub('\\]','',db_single$d1_up20)
db_single$d1_up20 <- gsub(',','-',db_single$d1_up20)

# add original and enlarged interval to conversion table
unique(db_single$d1_up20)
Lim["d1_up20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.01-0.464","1e-05-0.464","9.55-10","9.55-10")

a <- db_single[, c("Feature.ID", "d1_up20")]

# export
dir.create('./d1_up20')
write.table(a,'./d1_up20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interval
int <- unique(db_single$d1_up20)
int <- as.numeric(ifelse(substr(int,2,2)=='e',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d1_up20)

db_single$d1_up20 <- factor(db_single$d1_up20, 
                            levels = names(int)[order(int)])
pdf('./d1_up20/numberOfSequencesPerInterval.pdf', width = 10)
ggplot(db_single, aes(x=d1_up20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# bin in 30 non-empty intervals (add column with range)
db_single$d1_up30 <- cut(db_single$d1_up, breaks =41) 
db_single$d1_up30 <- gsub('\\(','',db_single$d1_up30)
db_single$d1_up30 <- gsub('\\]','',db_single$d1_up30)
db_single$d1_up30 <- gsub(',','-',db_single$d1_up30)

# add original and enlarged interval to conversion table
unique(db_single$d1_up30)
Lim["d1_up30",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.01-0.254","1e-05-0.254","9.76-10","9.76-10")

a <- db_single[, c("Feature.ID", "d1_up30")]

# export
dir.create('./d1_up30')
write.table(a,'./d1_up30/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interval
int <- unique(db_single$d1_up30)
int <- as.numeric(ifelse(substr(int,2,2)=='e',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d1_up30)

db_single$d1_up30 <- factor(db_single$d1_up30, 
                            levels = names(int)[order(int)])
pdf('./d1_up30/numberOfSequencesPerInterval.pdf', width = 15)
ggplot(db_single, aes(x=d1_up30)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

### d2_lo ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$d2_lo), c("Feature.ID", "d2_lo")]

# plot histogram
ggplot(db_single, aes(x=d2_lo)) + 
  geom_histogram(binwidth=0.2)
range(db_single$d2_lo)
# histogram with the values < 750
range(db_single[db_single$d2_lo<750,'d2_lo'])
ggplot(db_single[db_single$d2_lo<750,], aes(x=d2_lo)) + 
  geom_histogram(binwidth=0.2)

# there are empty bins because there are gaps in range when split in equal intervals
# split such that there is the desired number of bins that are not empty

# bin in 5 non-empty intervals (add column with range)
db_single$d2_lo5 <- cut(db_single$d2_lo, breaks =16) 
db_single$d2_lo5 <- gsub('\\(','',db_single$d2_lo5)
db_single$d2_lo5 <- gsub('\\]','',db_single$d2_lo5)
db_single$d2_lo5 <- gsub(',','-',db_single$d2_lo5)

# add original and enlarged interval to conversion table
unique(db_single$d2_lo5)
Lim["d2_lo5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.012-46.9","-0.738-46.9","703-750","703-751")

a <- db_single[, c("Feature.ID", "d2_lo5")]

# export
dir.create('./d2_lo5')
write.table(a,'./d2_lo5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interval
int <- unique(db_single$d2_lo5)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_lo5)

db_single$d2_lo5 <- factor(db_single$d2_lo5, 
                            levels = names(int)[order(int)])
pdf('./d2_lo5/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=d2_lo5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 10 non-empty intervals (add column with range)
db_single$d2_lo10 <- cut(db_single$d2_lo, breaks =63) 
db_single$d2_lo10 <- gsub('\\(','',db_single$d2_lo10)
db_single$d2_lo10 <- gsub('\\]','',db_single$d2_lo10)
db_single$d2_lo10 <- gsub(',','-',db_single$d2_lo10)


# original and enlarged interval
unique(db_single$d2_lo10)
Lim["d2_lo10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.012-11.9","-0.738-11.9","738-750","738-751")

a <- db_single[, c("Feature.ID", "d2_lo10")]

# export
dir.create('./d2_lo10')
write.table(a,'./d2_lo10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervalsby upper border of interval
int <- unique(db_single$d2_lo10)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_lo10)

db_single$d2_lo10 <- factor(db_single$d2_lo10, 
                           levels = names(int)[order(int)])
pdf('./d2_lo10/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=d2_lo10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 20 non-empty intervals (add column with range)
db_single$d2_lo20 <- cut(db_single$d2_lo, breaks =251) 
db_single$d2_lo20 <- gsub('\\(','',db_single$d2_lo20)
db_single$d2_lo20 <- gsub('\\]','',db_single$d2_lo20)
db_single$d2_lo20 <- gsub(',','-',db_single$d2_lo20)

# add original and enlarged interval to conversion table
unique(db_single$d2_lo20)
Lim["d2_lo20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.012-3","-0.738-3","747-750","747-751")

a <- db_single[, c("Feature.ID", "d2_lo20")]

# export
dir.create('./d2_lo20')
write.table(a,'./d2_lo20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interval
int <- unique(db_single$d2_lo20)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_lo20)

db_single$d2_lo20 <- factor(db_single$d2_lo20, 
                            levels = names(int)[order(int)])
pdf('./d2_lo20/numberOfSequencesPerInterval.pdf', width = 10)
ggplot(db_single, aes(x=d2_lo20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 30 non-empty intervals (add column with range)
db_single$d2_lo30 <- cut(db_single$d2_lo, breaks =515) 
db_single$d2_lo30 <- gsub('\\(','',db_single$d2_lo30)
db_single$d2_lo30 <- gsub('\\]','',db_single$d2_lo30)
db_single$d2_lo30 <- gsub(',','-',db_single$d2_lo30)

# add original and enlarged interval to conversion table
unique(db_single$d2_lo30)
Lim["d2_lo30",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.012-1.47","-0.738-1.47","749-750","749-751")

a <- db_single[, c("Feature.ID", "d2_lo30")]

# export
dir.create('./d2_lo30')
write.table(a,'./d2_lo30/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interval
int <- unique(db_single$d2_lo30)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_lo30)

db_single$d2_lo30 <- factor(db_single$d2_lo30, 
                            levels = names(int)[order(int)])
pdf('./d2_lo30/numberOfSequencesPerInterval.pdf', width = 15)
ggplot(db_single, aes(x=d2_lo30)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

### d2_up ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$d2_up), c("Feature.ID", "d2_up")]

# plot histogram
ggplot(db_single, aes(x=d2_up)) + 
  geom_histogram(binwidth=0.2)
range(db_single$d2_up)

# there are empty bins because there are gaps in range when split in equal intervals
# split such that there is the desired number of bins that are not empty

# bin in 5 non-empty intervals (add column with range)
db_single$d2_up5 <- cut(db_single$d2_up, breaks =5) 
db_single$d2_up5 <- gsub('\\(','',db_single$d2_up5)
db_single$d2_up5 <- gsub('\\]','',db_single$d2_up5)
db_single$d2_up5 <- gsub(',','-',db_single$d2_up5)

# add original and enlarged interval to conversion table
unique(db_single$d2_up5)
Lim["d2_up5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.014-19.2","-0.082-19.2","76.8-96","76.8-96.1")

a <- db_single[, c("Feature.ID", "d2_up5")]

# export
dir.create('./d2_up5')
write.table(a,'./d2_up5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interval 
int <- unique(db_single$d2_up5)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_up5)

db_single$d2_up5 <- factor(db_single$d2_up5, 
                            levels = names(int)[order(int)])
pdf('./d2_up5/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=d2_up5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 10 non-empty intervals (add column with range)
db_single$d2_up10 <- cut(db_single$d2_up, breaks =11) 
db_single$d2_up10 <- gsub('\\(','',db_single$d2_up10)
db_single$d2_up10 <- gsub('\\]','',db_single$d2_up10)
db_single$d2_up10 <- gsub(',','-',db_single$d2_up10)

# add original and enlarged interval to conversion table
unique(db_single$d2_up10)
Lim["d2_up10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.014-8.74","-0.082-8.74","87.3-96","87.3-96.1")

a <- db_single[, c("Feature.ID", "d2_up10")]

# export
dir.create('./d2_up10')
write.table(a,'./d2_up10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interval
int <- unique(db_single$d2_up10)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_up10)

db_single$d2_up10 <- factor(db_single$d2_up10, 
                           levels = names(int)[order(int)])
pdf('./d2_up10/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=d2_up10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 20 non-empty intervals (add column with range)
db_single$d2_up20 <- cut(db_single$d2_up, breaks =24) 
db_single$d2_up20 <- gsub('\\(','',db_single$d2_up20)
db_single$d2_up20 <- gsub('\\]','',db_single$d2_up20)
db_single$d2_up20 <- gsub(',','-',db_single$d2_up20)

# add original and enlarged interval to conversion table
unique(db_single$d2_up20)
Lim["d2_up20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.014-4.01","-0.082-4.01","92-96","92-96.1")

a <- db_single[, c("Feature.ID", "d2_up20")]

# export
dir.create('./d2_up20')
write.table(a,'./d2_up20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervalsbyupper border of interval
int <- unique(db_single$d2_up20)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_up20)

db_single$d2_up20 <- factor(db_single$d2_up20, 
                            levels = names(int)[order(int)])
pdf('./d2_up20/numberOfSequencesPerInterval.pdf', width = 10)
ggplot(db_single, aes(x=d2_up20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 30 non-empty intervals (add column with range)
db_single$d2_up30 <- cut(db_single$d2_up, breaks =59) 
db_single$d2_up30 <- gsub('\\(','',db_single$d2_up30)
db_single$d2_up30 <- gsub('\\]','',db_single$d2_up30)
db_single$d2_up30 <- gsub(',','-',db_single$d2_up30)

# add original and enlarged interval to conversion table
unique(db_single$d2_up30)
Lim["d2_up30",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.014-1.64","-0.082-1.64","94.4-96","94.4-96.1")

a <- db_single[, c("Feature.ID", "d2_up30")]

# export
dir.create('./d2_up30')
write.table(a,'./d2_up30/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border 
int <- unique(db_single$d2_up30)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$d2_up30)

db_single$d2_up30 <- factor(db_single$d2_up30, 
                            levels = names(int)[order(int)])
pdf('./d2_up30/numberOfSequencesPerInterval.pdf', width = 15)
ggplot(db_single, aes(x=d2_up30)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()



### doubling_h ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$doubling_h), c("Feature.ID", "doubling_h")]

# plot histogram
ggplot(db_single, aes(x=doubling_h)) + 
  geom_histogram(binwidth=1)
range(db_single$doubling_h)

# there are empty bins because there are gaps in range when split in equal intervals
# split such that there is the desired number of bins that are not empty

# bin in 5 non-empty intervals (add column with range)
db_single$doubling_h5 <- cut(db_single$doubling_h, breaks =5) 
db_single$doubling_h5 <- gsub('\\(','',db_single$doubling_h5)
db_single$doubling_h5 <- gsub('\\]','',db_single$doubling_h5)
db_single$doubling_h5 <- gsub(',','-',db_single$doubling_h5)

# add original and enlarged interval to conversion table
unique(db_single$doubling_h5)
Lim["doubling_h5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.16-86.7","-0.273-86.7","346-433.01","346-433")

a <- db_single[, c("Feature.ID", "doubling_h5")]

# export
dir.create('./doubling_h5')
write.table(a,'./doubling_h5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$doubling_h5)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$doubling_h5)

db_single$doubling_h5 <- factor(db_single$doubling_h5, 
                            levels = names(int)[order(int)])
pdf('./doubling_h5/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=doubling_h5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 10 non-empty intervals (add column with range)
db_single$doubling_h10 <- cut(db_single$doubling_h, breaks =12)
db_single$doubling_h10 <- gsub('\\(','',db_single$doubling_h10)
db_single$doubling_h10 <- gsub('\\]','',db_single$doubling_h10)
db_single$doubling_h10 <- gsub(',','-',db_single$doubling_h10)
# add original and enlarged interval to conversion table
unique(db_single$doubling_h10)
Lim["doubling_h10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.16-36.2","-0.273-36.2","397-433.01","397-433")

a <- db_single[, c("Feature.ID", "doubling_h10")]

# export
dir.create('./doubling_h10')
write.table(a,'./doubling_h10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$doubling_h10)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$doubling_h10)

db_single$doubling_h10 <- factor(db_single$doubling_h10, 
                                levels = names(int)[order(int)])
pdf('./doubling_h10/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=doubling_h10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# bin in 20 non-empty intervals (add column with range)
db_single$doubling_h20 <- cut(db_single$doubling_h, breaks =38) 
db_single$doubling_h20 <- gsub('\\(','',db_single$doubling_h20)
db_single$doubling_h20 <- gsub('\\]','',db_single$doubling_h20)
db_single$doubling_h20 <- gsub(',','-',db_single$doubling_h20)

# add original and enlarged interval to conversion table
unique(db_single$doubling_h20)
Lim["doubling_h20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.16-11.6","-0.273-11.6","422-433.01","422-433")

a <- db_single[, c("Feature.ID", "doubling_h20")]

# export
dir.create('./doubling_h20')
write.table(a,'./doubling_h20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$doubling_h20)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$doubling_h20)

db_single$doubling_h20 <- factor(db_single$doubling_h20, 
                                 levels = names(int)[order(int)])
pdf('./doubling_h20/numberOfSequencesPerInterval.pdf', width = 10)
ggplot(db_single, aes(x=doubling_h20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 30 non-empty intervals (add column with range)
db_single$doubling_h30 <- cut(db_single$doubling_h, breaks =77) 
db_single$doubling_h30 <- gsub('\\(','',db_single$doubling_h30)
db_single$doubling_h30 <- gsub('\\]','',db_single$doubling_h30)
db_single$doubling_h30 <- gsub(',','-',db_single$doubling_h30)

# add original and enlarged interval to conversion table
unique(db_single$doubling_h30)
Lim["doubling_h30",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.16-5.78","-0.273-5.78","427-433.01","427-433")

a <- db_single[, c("Feature.ID", "doubling_h30")]

# export
dir.create('./doubling_h30')
write.table(a,'./doubling_h30/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$doubling_h30)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$doubling_h30)

db_single$doubling_h30 <- factor(db_single$doubling_h30, 
                                 levels = names(int)[order(int)])
pdf('./doubling_h30/numberOfSequencesPerInterval.pdf', width = 15)
ggplot(db_single, aes(x=doubling_h30)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 40 non-empty intervals (add column with range)
db_single$doubling_h40 <- cut(db_single$doubling_h, breaks =124)
# gives 41 non-empty intervals (123 breaks gives only 39)
db_single$doubling_h40 <- gsub('\\(','',db_single$doubling_h40)
db_single$doubling_h40 <- gsub('\\]','',db_single$doubling_h40)
db_single$doubling_h40 <- gsub(',','-',db_single$doubling_h40)

# add original and enlarged interval to conversion table
unique(db_single$doubling_h40)
Lim["doubling_h40",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.16-3.65","-0.273,3.65","430-433.01","430-433")

a <- db_single[, c("Feature.ID", "doubling_h40")]

# export
dir.create('./doubling_h40')
write.table(a,'./doubling_h40/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border 
int <- unique(db_single$doubling_h40)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$doubling_h40)

db_single$doubling_h40 <- factor(db_single$doubling_h40, 
                                 levels = names(int)[order(int)])
pdf('./doubling_h40/numberOfSequencesPerInterval.pdf', width = 15)
ggplot(db_single, aes(x=doubling_h40)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# bin in 50 non-empty intervals (add column with range)
db_single$doubling_h50 <- cut(db_single$doubling_h, breaks =187) 
# gives 51 non-empty intervals (186 breaks gives only 49)
db_single$doubling_h50 <- gsub('\\(','',db_single$doubling_h50)
db_single$doubling_h50 <- gsub('\\]','',db_single$doubling_h50)
db_single$doubling_h50 <- gsub(',','-',db_single$doubling_h50)

# add original and enlarged interval to conversion table
unique(db_single$doubling_h50)
Lim["doubling_h50",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.16-2.47","-0.273,2.47","431-433.01","431-433")

a <- db_single[, c("Feature.ID", "doubling_h50")]

# export
dir.create('./doubling_h50')
write.table(a,'./doubling_h50/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$doubling_h50)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$doubling_h50)

db_single$doubling_h50 <- factor(db_single$doubling_h50, 
                                 levels = names(int)[order(int)])
pdf('./doubling_h50/numberOfSequencesPerInterval.pdf', width = 15)
ggplot(db_single, aes(x=doubling_h50)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

### genome_size ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$genome_size), c("Feature.ID", "genome_size")]

# plot histogram
ggplot(db_single, aes(x=genome_size)) + 
  geom_histogram(binwidth=10000)
range(db_single$genome_size)

# plot histogram for workflow summary figure
pdf('./genomeSize_histogram.pdf', width = 5, height = 3.5)
ggplot(db_single, aes(x=genome_size)) + 
  geom_histogram(binwidth=500000) + plot.theme1
dev.off()
range(db_single$genome_size)

# there are empty bins because there are gaps in range when split in equal intervals
# split such that there is the desired number of bins that are not empty

# bin in 5 non-empty intervals (add column with range)
db_single$genome_size5 <- cut(db_single$genome_size, breaks =5)
db_single$genome_size5 <- gsub('\\(','',db_single$genome_size5)
db_single$genome_size5 <- gsub('\\]','',db_single$genome_size5)
db_single$genome_size5 <- gsub(',','-',db_single$genome_size5)

# add original and enlarged interval to conversion table
unique(db_single$genome_size5)
Lim["genome_size5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("134927-3.32e+06","1.19e+05-3.32e+06","1.29e+07-16040677","1.29e+07-1.61e+07")

a <- db_single[, c("Feature.ID", "genome_size5")]

# export
dir.create('./genome_size5')
write.table(a,'./genome_size5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interal to sort 
int <- unique(db_single$genome_size5)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$genome_size5)

db_single$genome_size5 <- factor(db_single$genome_size5, 
                                 levels = names(int)[order(int)])
pdf('./genome_size5/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=genome_size5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# bin in 10 non-empty intervals (add column with range)
db_single$genome_size10 <- cut(db_single$genome_size, breaks =10) 
db_single$genome_size10 <- gsub('\\(','',db_single$genome_size10)
db_single$genome_size10 <- gsub('\\]','',db_single$genome_size10)
db_single$genome_size10 <- gsub(',','-',db_single$genome_size10)

# add original and enlarged interval to conversion table
unique(db_single$genome_size10)
Lim["genome_size10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("134927-1.73e+06","1.19e+05-1.73e+06","1.45e+07-16040677","1.45e+07-1.61e+07")

a <- db_single[, c("Feature.ID", "genome_size10")]

# export
dir.create('./genome_size10')
write.table(a,'./genome_size10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$genome_size10)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$genome_size10)

db_single$genome_size10 <- factor(db_single$genome_size10, 
                                 levels = names(int)[order(int)])
pdf('./genome_size10/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=genome_size10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 20 non-empty intervals (add column with range)
db_single$genome_size20 <- cut(db_single$genome_size, breaks =20) 
db_single$genome_size20 <- gsub('\\(','',db_single$genome_size20)
db_single$genome_size20 <- gsub('\\]','',db_single$genome_size20)

# add original and enlarged interval to conversion table
unique(db_single$genome_size20)
Lim["genome_size20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("134927-9.3e+05","1.19e+05-9.3e+05","1.52e+07-16040677","1.52e+07-1.61e+07")

a <- db_single[, c("Feature.ID", "genome_size20")]

# export
dir.create('./genome_size20')
write.table(a,'./genome_size20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$genome_size20)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$genome_size20)

db_single$genome_size20 <- factor(db_single$genome_size20, 
                                  levels = names(int)[order(int)])
pdf('./genome_size20/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=genome_size20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

### rRNA16S_genes ********************************************************************************************************************************************************** ####

# make with 5 bins, 10 bins and with exact numbers 

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$rRNA16S_genes), c("Feature.ID", "rRNA16S_genes")]

# plot histogram
ggplot(db_single, aes(x=rRNA16S_genes)) + 
  geom_histogram(binwidth=1)
range(db_single$rRNA16S_genes)

# bin in 5 non-empty intervals (add column with range)
db_single$rRNA16S_genes5 <- cut(db_single$rRNA16S_genes, breaks =5)
db_single$rRNA16S_genes5 <- gsub('\\(','',db_single$rRNA16S_genes5)
db_single$rRNA16S_genes5 <- gsub('\\]','',db_single$rRNA16S_genes5)
db_single$rRNA16S_genes5 <- gsub(',','-',db_single$rRNA16S_genes5)

# export
unique(db_single$rRNA16S_genes5)
Lim["rRNA16S_genes5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("1-4.2","0.984-4.2","13.8-17","13.8-17")

a <- db_single[, c("Feature.ID", "rRNA16S_genes5")]

# export
dir.create('./rRNA16S_genes5')
write.table(a,'./rRNA16S_genes5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$rRNA16S_genes5)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$rRNA16S_genes5)

db_single$rRNA16S_genes5 <- factor(db_single$rRNA16S_genes5, 
                                  levels = names(int)[order(int)])
pdf('./rRNA16S_genes5/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=rRNA16S_genes5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 10 non-empty intervals (add column with range)
db_single$rRNA16S_genes10 <- cut(db_single$rRNA16S_genes, breaks =10)
db_single$rRNA16S_genes10 <- gsub('\\(','',db_single$rRNA16S_genes10)
db_single$rRNA16S_genes10 <- gsub('\\]','',db_single$rRNA16S_genes10)
db_single$rRNA16S_genes10 <- gsub(',','-',db_single$rRNA16S_genes10)

# add original and enlarged interval to conversion table
unique(db_single$rRNA16S_genes10)
Lim["rRNA16S_genes10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("1-2.6","0.984-2.6","13.8-17","13.8-17")

a <- db_single[, c("Feature.ID", "rRNA16S_genes10")]

# export
dir.create('./rRNA16S_genes10')
write.table(a,'./rRNA16S_genes10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border
int <- unique(db_single$rRNA16S_genes10)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$rRNA16S_genes10)

db_single$rRNA16S_genes10 <- factor(db_single$rRNA16S_genes10, 
                                   levels = names(int)[order(int)])
pdf('./rRNA16S_genes10/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=rRNA16S_genes10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# exact values (rounded to integers)
# trait database is aggregated for species resulting in non-integer values for 16S gene copies; round to integers
a <- db_single[, c("Feature.ID", "rRNA16S_genes")]
a$rRNA16S_genes <- round(a$rRNA16S_genes)

# export
dir.create('./rRNA16S_genes_exact_int')
write.table(a,'./rRNA16S_genes_exact_int/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)


### optimum_ph ********************************************************************************************************************************************************** ####

# subset with only accession number and value of trait
db_single <- dbAll[!is.na(dbAll$optimum_ph), c("Feature.ID", "optimum_ph")]

# plot histogram
ggplot(db_single, aes(x=optimum_ph)) + 
  geom_histogram(binwidth=0.5)
range(db_single$optimum_ph)

# there are empty bins because there are gaps in range when split in equal intervals
# split such that there is the desired number of bins that are not empty

# bin in 5 non-empty intervals (add column with range)
db_single$optimum_ph5 <- cut(db_single$optimum_ph, breaks =5) 
db_single$optimum_ph5 <- gsub('\\(','',db_single$optimum_ph5)
db_single$optimum_ph5 <- gsub('\\]','',db_single$optimum_ph5)
db_single$optimum_ph5 <- gsub(',','-',db_single$optimum_ph5)
# add original and enlarged interval to conversion table
unique(db_single$optimum_ph5)
Lim["optimum_ph5",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.75-2.85","0.74-2.85","9.15-11.25","9.15-11.3")
a <- db_single[, c("Feature.ID", "optimum_ph5")]

# export
dir.create('./optimum_ph5')
write.table(a,'./optimum_ph5/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border of interal to sort 
int <- unique(db_single$optimum_ph5)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$optimum_ph5)

db_single$optimum_ph5 <- factor(db_single$optimum_ph5, 
                                    levels = names(int)[order(int)])
pdf('./optimum_ph5/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=optimum_ph5)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 10 non-empty intervals (add column with range)
db_single$optimum_ph10 <- cut(db_single$optimum_ph, breaks =10) 
db_single$optimum_ph10 <- gsub('\\(','',db_single$optimum_ph10)
db_single$optimum_ph10 <- gsub('\\]','',db_single$optimum_ph10)
db_single$optimum_ph10 <- gsub(',','-',db_single$optimum_ph10)

# add original and enlarged interval to conversion table
unique(db_single$optimum_ph10)
Lim["optimum_ph10",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.75-1.8","0.74-1.8","10.2-11.25","10.2-11.3")

a <- db_single[, c("Feature.ID", "optimum_ph10")]

# export
dir.create('./optimum_ph10')
write.table(a,'./optimum_ph10/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border  
int <- unique(db_single$optimum_ph10)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$optimum_ph10)

db_single$optimum_ph10 <- factor(db_single$optimum_ph10, 
                                levels = names(int)[order(int)])
pdf('./optimum_ph10/numberOfSequencesPerInterval.pdf', width = 7)
ggplot(db_single, aes(x=optimum_ph10)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# bin in 20 non-empty intervals (add column with range)
db_single$optimum_ph20 <- cut(db_single$optimum_ph, breaks =21) 
db_single$optimum_ph20 <- gsub('\\(','',db_single$optimum_ph20)
db_single$optimum_ph20 <- gsub('\\]','',db_single$optimum_ph20)
db_single$optimum_ph20 <- gsub(',','-',db_single$optimum_ph20)

# add original and enlarged interval to conversion table
unique(db_single$optimum_ph20)
Lim["optimum_ph20",c("LI_or","LI_cut","HI_or","HI_cut")] <-
  c("0.75-1.25","0.74-1.25","10.8-11.25","10.8-11.3")

a <- db_single[, c("Feature.ID", "optimum_ph20")]

# export
dir.create('./optimum_ph20')
write.table(a,'./optimum_ph20/traits.tsv', row.names = F, quote = F,sep = '\t',col.names = F)

# plot number of sequences in each interval
# sort Intervals by upper border 
int <- unique(db_single$optimum_ph20)
int <- as.numeric(ifelse(substr(int,1,1)=='-',gsub('(-.*?)(-.*)','\\1', int),gsub('-.*','', int)))
names(int) <- unique(db_single$optimum_ph20)

db_single$optimum_ph20 <- factor(db_single$optimum_ph20, 
                                 levels = names(int)[order(int)])
pdf('./optimum_ph20/numberOfSequencesPerInterval.pdf', width = 10)
ggplot(db_single, aes(x=optimum_ph20)) +
  geom_histogram(stat = "count") +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  ggtitle("Number of sequences per interval") +
  plot.theme1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# export table with highest and lowest intervals to reconvert limits
write.csv(Lim,'./IntervalLimits_Reconversion.csv')




