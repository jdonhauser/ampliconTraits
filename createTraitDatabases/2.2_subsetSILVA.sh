
# subset silva database sequences to the accession numbers from Madin's database ******************************************************************************************************************
# list with accession numbers idsList.txt created in silva2ncbi.R
# working directory should be set to 2_matchTaxids

# path to sequences
path="../1_prepareSilva/silva-138-ssu-seqs-341f-806r-uniq"

# install and activate seqkit to subset sequences based on list
conda create --name seqkit -c bioconda
conda activate seqkit

# rename fasta header to contain only the accession.start.stop 
# remove taxonomy from fasta header
awk '/^>/{print $1; next}{print}' $path/dna-sequences.fasta > ./dna-sequences-renamed.fasta

# use seqkit grep to find sequences corresponding to accession in id_list.txt
seqkit grep -n -f id_list.txt dna-sequences-renamed.fasta > silva138_SSURef_inMadinDB.fasta

