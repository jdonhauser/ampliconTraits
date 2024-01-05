# merge Silva taxonomy with Madin's database and subset sequences to IDs contained ************************************************** ####
# crossmapping based on ncbi taxids using taxmap_embl-ebi_ena_ssu_ref_138.txt (contains ncbi taxids for silva sequences);
# Downloaded from: https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/taxonomy/taxmap_embl-ebi_ena_ssu_ref_138.txt.gz
# use species name to map additional database entries not matched by taxid


## Import files ******************************************************************************************************************** ####
# working directory should be set to 2_matchTaxids

# read mapping files
taxmap <- read.csv('./taxmap_embl-ebi_ena_ssu_ref_138.txt', sep = '\t', header = T) 
# read silva taxonomy (containing accession.start.stop and taxonmic path separated by; corresponds to fasta headers)
silvatax <- read.csv('../1_prepareSilva/silva-138-ssu-tax-341f-806r-derep-uniq/taxonomy.tsv',sep = '\t')

# split taxonomy into several columns to match by species name
library(tidyr)
silvatax <- separate(data = silvatax, col = Taxon, into = c("kingdom", "phylum","class","order","family","genus","species"),sep = ';')
# prepare column species in silva to match formatting of Madin's database
silvatax$species <- gsub('s__','', silvatax$species)
silvatax$species <- gsub('_',' ', silvatax$species)
silvatax$species <- gsub('^ ','', silvatax$species)

# extract accession number from silva taxonomy (remove part with start position and stop position)
silvatax$primaryAccession <- gsub('\\..*','',silvatax$Feature.ID)

### combine taxonomy with mapping table
# primaryAccession is not unique due to different start stop position for same accession. Therefore use 
# feature.ID (accession.start.stop) to merge tables

# create feature.ID column in taxmap by combining accession.start.stop
taxmap$Feature.ID <- paste(taxmap$primaryAccession,taxmap$start,taxmap$stop, sep = '.')

# merge by Feature.ID
silva2ncbi <- merge(silvatax,taxmap, by = "Feature.ID", all.x = T) 
any(is.na(silva2ncbi$ncbi_taxonid)) # FALSE

# import trait databse by Madin
# downloaded from https://figshare.com/ndownloader/files/26019008
db <- read.csv('./condensed_species_NCBI.csv')

## Combine based on taxid ********************************************************************************************************** ####

# make subset with silva database containing only taxids in Madin db
silvaSub <- silva2ncbi[silva2ncbi$ncbi_taxonid%in%db$species_tax_id,]
# count number of taxids retrieved
length(unique(silvaSub$ncbi_taxonid)) # 13426

# make a subset with the species from Madin db that are not found in Silva and see how many traits this subset has annotated
db_notFound <- db[!(db$species_tax_id%in%silvaSub$ncbi_taxonid),]
# 1467 taxids

# count the number of non NA entries for each column (# traits annotated among the taxa not found)
apply(db_notFound,2, function(x)length(which(!is.na(x))))


## Combine based on species name *************************************************************************************************** ####
# check if some of the not found species can be matched with silva using the species name
db_matchSpName <- db_notFound[db_notFound$species%in%silvatax$species,]


# 497 taxa found 
# as all silva entries have a taxid there must be a disagreement between the taxid in the Madin database and in silva
# combine the matched by species with silva2ncbi to check taxid in each ressource
db_matchSpName <- merge(db_matchSpName,silva2ncbi,by = 'species', suffixes = c('.Madin','.silva')) # if there are multiple matches for same taxid rows become multiplicated
db_matchSpName <- db_matchSpName[,-c(9:79)]
length(unique(db_matchSpName$species))
length(unique(db_matchSpName$species_tax_id))

# export csv to inspect and compare taxids from the two resources manually
write.csv(db_matchSpName, './taxaMatchedSpeciesNameButNotTaxid.csv')

# in these cases the silva taxid belongs to a strain or subspecies within the species of which the taxid is used in Madins database
# exception: 
# unidentified bacterium: taxonomy from genus to class missing in Madin's DB => remove
#  Bacillales bacterium (family and genus are missing in Madin's DB but different genera in silva: => remove)

### check properties for taxids in Madins database not found in Silva with either method (number of traits annotated, species names) ****** ####
db_notFound2 <- db_notFound[!(db_notFound$species_tax_id%in%db_matchSpName$species_tax_id),]
# count the number of non NA entries for each column (traits annotated among the non found taxa)
apply(db_notFound2,2, function(x)length(which(!is.na(x))))
unique(db_notFound2$species)

## merge tables mapped by taxid and by species name with silva and then the two parts together ******************************************** ####

# combine each of the parts of Madin's db (mapped by taxid and species name respectively) with sivla to obtain trait annotations for each sequences
# need to be combined separately because they were merged based on different columns
# if for matching columns there are multiple entries on one side these become multiplicated

# part matched by accession number
db_inSilva <- db[(db$species_tax_id%in%silvaSub$ncbi_taxonid),]
# merge with silva
db_inSilva <- merge(db_inSilva, silvaSub, by.x = 'species_tax_id',by.y = 'ncbi_taxonid', all.x = T,sort = F,suffixes = c('.Madin','.silva'))

# part matched by species name
db_matchSpName <- db_notFound[db_notFound$species%in%silva2ncbi$species,]
db_matchSpName <- db_matchSpName[db_matchSpName$species!='unidentified bacterium',]
db_matchSpName <- db_matchSpName[db_matchSpName$species!='Bacillales bacterium',]
# merge with silva
db_matchSpName <- merge(db_matchSpName,silva2ncbi[!(silva2ncbi$ncbi_taxonid%in%db_inSilva$species_tax_id),],by = 'species',sort = F,suffixes = c('.Madin','.silva'))

# duplicate columns that were united with merging with their original name and sort columns in same way in both data frames
# db_insilva: re-add column ncbi_taxonid
db_inSilva$ncbi_taxonid <- db_inSilva$species_tax_id

# db_matchSpName: rename column species to species.Madin and re-add species silva
colnames(db_matchSpName)[1] <- 'species.Madin'
db_matchSpName$species.silva <- db_matchSpName$species.Madin

all(colnames(db_inSilva)%in%colnames(db_matchSpName)) # all column names present in both dataframes

# order db_matchSpName in the same way as db_inSilva
db_matchSpName <- db_matchSpName[,colnames(db_inSilva)]
all(colnames(db_inSilva)==colnames(db_matchSpName))

# combine parts of Madins database matched by taxid and by species names
db_inSilva <- rbind(db_inSilva,db_matchSpName)
                    
# total number of species from trait database that could be matched with SILVA
length(unique(db_inSilva$species_tax_id)) # 13921

## Export output ********************************************************************************************************************* ####
# export accession numbers (with start and stop to match fasta headers) that are found in Madin's DB
# to subset fasta file of Silva database based on list of accession number
write.table(db_inSilva$Feature.ID, './id_list.txt',row.names = F,col.names = F,quote = F)

# export whole database merged with Silva information to make a database for each trait
write.csv(db_inSilva, './dbMadin_merged_Silva.csv')

## number of species and sequences for each trait ************************************************************************************ ####
specSeq_perTrait <- data.frame(numberSequences=1:94)
rownames(specSeq_perTrait) <- colnames(db_inSilva)

# number of sequences per trait
specSeq_perTrait$numberSequences <- apply(db_inSilva,2,function(x){length(which(!(is.na(x))))})

#number of species per traits
a <- db_inSilva[!duplicated(db_inSilva$species_tax_id),]
specSeq_perTrait$numberSpecies <- apply(a,2,function(x){length(which(!(is.na(x))))})
  
# export
write.csv(specSeq_perTrait, './numberSpeciesAndSequences_perTrait.csv')
