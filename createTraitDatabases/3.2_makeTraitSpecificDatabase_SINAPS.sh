# subset sequences and format for SINAPS *****************************************************************************************
# SINAPS format: fasta file	with trait and its value in the header so the header takes the form ID;trait=value
# create a fasta file for each trait, containing only the sequences for which the trait is annotated
# working directory is 3_makeTraitSpecificDB/trait

# activate seqkit to subset sequences
conda activate seqkit

# create variable for path to fasta file with all sequences
path="../../2_matchTaxids"

# metabolism ***************************************************************************************************************
cd ./metabolism

# create list of ids with the sequence IDs annotated for trait to filter fasta file
awk '{ print $1 }' traits.tsv > id_list.txt
# subset fasta
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_metabolism.fasta

# format for sinaps
# create fasta headers from traits file to replace the headers in fasta file afterwards
# in output file: 1.column: current fasta header (to be used as index); 2. column: replacement (sequenceID;trait=value)
awk '{print ">" $1 "\t>" $1 ";metabolism=" $2}' traits.tsv > seq_headers
# replace fasta headers
# sequences in fasta and seq_headers are not in the same order; use indexed array to add correct header to each sequence
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_metabolism.fasta > silva138_SSURef_inMadinDB_metabolism_renamed.fasta 

# cell shape *****************************************************************************************************
cd ../cell_shape

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_cellShape.fasta

# format for sinaps
# file with two columns: 1. current fasta header (to be used as index); 2. replacement
awk '{print ">" $1 "\t>" $1 ";cellshape=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_cellShape.fasta > silva138_SSURef_inMadinDB_cellShape_renamed.fasta 

# gram stain *******************************************************************************************************************
cd ../gram_stain

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_gramstain.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";gramstain=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_gramstain.fasta > silva138_SSURef_inMadinDB_gramstain_renamed.fasta 

# motility **************************************************************************************************************
cd ../motility

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_motility.fasta

# format for sinaps
# file with two columns: 1. current fasta header (to be used as index); 2. replacement
awk '{print ">" $1 "\t>" $1 ";motility=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_motility.fasta > silva138_SSURef_inMadinDB_motility_renamed.fasta 

# range_salinity ***************************************************************************************************
cd ../range_salinity

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_rangesalinity.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";rangesalinity=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_rangesalinity.fasta > silva138_SSURef_inMadinDB_rangesalinity_renamed.fasta 

# range_tmp ************************************************************************************************************
cd ../range_tmp

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_rangetmp.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";rangetmp=" $2}' traits.tsv > seq_headers
# remove taxonomy from fasta header
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_rangetmp.fasta > silva138_SSURef_inMadinDB_rangetmp_renamed.fasta 

# sporulation ***************************************************************************************************************
cd ../sporulation

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_sporulation.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";sporulation=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_sporulation.fasta > silva138_SSURef_inMadinDB_sporulation_renamed.fasta 


# d1_lo5 *************************************************************************************************************
cd ../d1_lo5

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_lo5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_lo5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_lo5.fasta > silva138_SSURef_inMadinDB_d1_lo5_renamed.fasta 

# d1_lo10 ***************************************************************************************************************
cd ../d1_lo10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_lo10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_lo10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_lo10.fasta > silva138_SSURef_inMadinDB_d1_lo10_renamed.fasta 

# d1_lo20 ***********************************************************************************************************
cd ../d1_lo20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_lo20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_lo20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_lo20.fasta > silva138_SSURef_inMadinDB_d1_lo20_renamed.fasta 

# d1_lo30 *******************************************************************************************************
cd ../d1_lo30

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_lo30.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_lo30=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_lo30.fasta > silva138_SSURef_inMadinDB_d1_lo30_renamed.fasta 

#  d1_up5 ****************************************************************************************************************
cd ../d1_up5

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_up5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_up5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_up5.fasta > silva138_SSURef_inMadinDB_d1_up5_renamed.fasta 

# d1_up10 ******************************************************************************************************************
cd ../d1_up10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_up10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_up10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_up10.fasta > silva138_SSURef_inMadinDB_d1_up10_renamed.fasta 

# d1_up20 **********************************************************************************************************
cd ../d1_up20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_up20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_up20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_up20.fasta > silva138_SSURef_inMadinDB_d1_up20_renamed.fasta 

# d1_up30 ************************************************************************************************************
cd ../d1_up30

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d1_up30.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d1_up30=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d1_up30.fasta > silva138_SSURef_inMadinDB_d1_up30_renamed.fasta 

#  d2_lo5 ******************************************************************************************************************
cd ../d2_lo5

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_lo5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_lo5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_lo5.fasta > silva138_SSURef_inMadinDB_d2_lo5_renamed.fasta 

# d2_lo10 ******************************************************************************************************************
cd ../d2_lo10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_lo10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_lo10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_lo10.fasta > silva138_SSURef_inMadinDB_d2_lo10_renamed.fasta 

# d2_lo20 *******************************************************************************************************
cd ../d2_lo20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_lo20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_lo20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_lo20.fasta > silva138_SSURef_inMadinDB_d2_lo20_renamed.fasta 

# d2_lo30 *********************************************************************************************************
cd ../d2_lo30

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_lo30.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_lo30=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_lo30.fasta > silva138_SSURef_inMadinDB_d2_lo30_renamed.fasta 


#  d2_up5 ***************************************************************************************************************
cd ../d2_up5

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_up5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_up5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_up5.fasta > silva138_SSURef_inMadinDB_d2_up5_renamed.fasta 

# d2_up10 *********************************************************************************************************
cd ../d2_up10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_up10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_up10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_up10.fasta > silva138_SSURef_inMadinDB_d2_up10_renamed.fasta 

# d2_up20 *************************************************************************************************************
cd ../d2_up20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_up20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_up20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_up20.fasta > silva138_SSURef_inMadinDB_d2_up20_renamed.fasta 

# d2_up30 **********************************************************************************************************
cd ../d2_up30

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_d2_up30.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";d2_up30=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_d2_up30.fasta > silva138_SSURef_inMadinDB_d2_up30_renamed.fasta 

# doubling_h5 ****************************************************************************************************************
cd ../doubling_h5

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_doubling_h5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";doubling_h5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_doubling_h5.fasta > silva138_SSURef_inMadinDB_doubling_h5_renamed.fasta 

# doubling_h10 ***********************************************************************************************************
cd ../doubling_h10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_doubling_h10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";doubling_h10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_doubling_h10.fasta > silva138_SSURef_inMadinDB_doubling_h10_renamed.fasta 

# doubling_h20 ****************************************************************************************************
cd ../doubling_h20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_doubling_h20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";doubling_h20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_doubling_h20.fasta >silva138_SSURef_inMadinDB_doubling_h20_renamed.fasta 

# doubling_h30 **************************************************************************************************************
cd ../doubling_h30

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_doubling_h30.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";doubling_h30=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_doubling_h30.fasta > silva138_SSURef_inMadinDB_doubling_h30_renamed.fasta 

#  doubling_h40 *********************************************************************************************************
cd ../doubling_h40

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_doubling_h40.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";doubling_h40=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_doubling_h40.fasta > silva138_SSURef_inMadinDB_doubling_h40_renamed.fasta 

# doubling_h50 **********************************************************************************************************
cd ../doubling_h50

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_doubling_h50.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";doubling_h50=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_doubling_h50.fasta > silva138_SSURef_inMadinDB_doubling_h50_renamed.fasta 

#  genome_size5 ***************************************************************************************************************
cd ../genome_size5

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_genome_size5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";genome_size5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_genome_size5.fasta > silva138_SSURef_inMadinDB_genome_size5_renamed.fasta 


# genome_size10 *********************************************************************************************************
cd ../genome_size10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_genome_size10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";genome_size10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_genome_size10.fasta > silva138_SSURef_inMadinDB_genome_size10_renamed.fasta 

# genome_size20 ************************************************************************************************************
cd ../genome_size20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_genome_size20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";genome_size20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_genome_size20.fasta > silva138_SSURef_inMadinDB_genome_size20_renamed.fasta 

# optimum_ph5 ***************************************************************************************************************
cd ../optimum_ph5

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_optimum_ph5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";optimum_ph5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_optimum_ph5.fasta > silva138_SSURef_inMadinDB_optimum_ph5_renamed.fasta 

# optimum_ph10 ****************************************************************************************************************
cd ../optimum_ph10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_optimum_ph10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";optimum_ph10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_optimum_ph10.fasta > silva138_SSURef_inMadinDB_optimum_ph10_renamed.fasta 

# optimum_ph20 **************************************************************************************************************
cd ../optimum_ph20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_optimum_ph20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";optimum_ph20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_optimum_ph20.fasta > silva138_SSURef_inMadinDB_optimum_ph20_renamed.fasta 

# optimum_tmp10 ************************************************************************************************************
cd ../optimum_tmp10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_optimum_tmp10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";optimum_tmp10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_optimum_tmp10.fasta > silva138_SSURef_inMadinDB_optimum_tmp10_renamed.fasta 


#  optimum_tmp20 ***********************************************************************************************************
cd ../optimum_tmp20

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_optimum_tmp20.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";optimum_tmp20=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_optimum_tmp20.fasta > silva138_SSURef_inMadinDB_optimum_tmp20_renamed.fasta 


# rRNA16S_genes5 ***********************************************************************************************************
cd ../rRNA16S_genes5
# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_rRNA16S_genes5.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";rRNA16S_genes5=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_rRNA16S_genes5.fasta > silva138_SSURef_inMadinDB_rRNA16S_genes5_renamed.fasta 


# rRNA16S_genes10 ***********************************************************************************************************
cd ../rRNA16S_genes10

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_rRNA16S_genes10.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";rRNA16S_genes10=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_rRNA16S_genes10.fasta > silva138_SSURef_inMadinDB_rRNA16S_genes10_renamed.fasta 

# rRNA16S_genes_exact_int *********************************************************************************************
cd ../rRNA16S_genes_exact_int

# subset fasta
awk '{ print $1 }' traits.tsv > id_list.txt
seqkit grep -n -f id_list.txt $path/silva138_SSURef_inMadinDB.fasta > silva138_SSURef_inMadinDB_rRNA16S_genes_exact_int.fasta

# format for sinaps
awk '{print ">" $1 "\t>" $1 ";rRNA16S_genes_exact_int=" $2}' traits.tsv > seq_headers
awk 'FNR==NR {f2[$1]=$2;next}/^>/{$1=f2[$1]}1' seq_headers  silva138_SSURef_inMadinDB_rRNA16S_genes_exact_int.fasta > silva138_SSURef_inMadinDB_rRNA16S_genes_exact_int_renamed.fasta 