######### Make SILVA SSURef 138 (NOT clustered at 99%) database trimmed to 341F 806R regions  ##########################

# installation of rescript in an existing qiime environment as described here https://github.com/bokulich-lab/RESCRIPt/ under option 2

# DAtabase preparation ******************************************************************************************************
# Download taxonomy files and import into qiime
# silva 138.1 does not have the corresponding files (taxmap_embl-ebi_ena_ssu_ref_138.txt.gz) to map the accesions to the taxids
# therefore v138 used
qiime rescript get-silva-data \
    --p-version '138' \
    --p-target 'SSURef' \
    --p-include-species-labels \
    --o-silva-sequences silva-138-ssu-rna-seqs.qza \
    --o-silva-taxonomy silva-138-ssu-tax.qza \
	--p-no-rank-propagation # this is not default; prevent usage of higher level ranks if lower levels are unknown
	
# silva sequences are rna => reverse transcribe to DNA
qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138-ssu-rna-seqs.qza \
    --o-dna-sequences silva-138-ssu-seqs.qza
	
# remove sequences with 5 or more ambiguous bases and any homopolymers > 8 bases
qiime rescript cull-seqs \
    --i-sequences silva-138-ssu-seqs.qza \
    --o-clean-sequences silva-138-ssu-seqs-cleaned.qza

#remove sequences that are to short; minimum length depending on taxonomy (Archaea <900; Bacteria < 1200, Eukaryota < 1400)
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138-ssu-seqs-cleaned.qza \
    --i-taxonomy silva-138-ssu-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138-ssu-seqs-filt.qza \
    --o-discarded-seqs silva-138-ssu-seqs-discard.qza

# dereplicate using --p-mode 'uniq' (default), keeps identical sequencess with different taxonomy
qiime rescript dereplicate \
    --i-sequences silva-138-ssu-seqs-filt.qza  \
    --i-taxa silva-138-ssu-tax.qza \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138-ssu-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138-ssu-tax-derep-uniq.qza
	
# make amplicon specific database (for 341F 806 R, primers should be 5' to 3')
# Fw CCTAYGGGRBGCASCAG
# Rev: GGACTACHVGGGTWTCTAAT
qiime feature-classifier extract-reads \
    --i-sequences silva-138-ssu-seqs-derep-uniq.qza \
    --p-f-primer CCTAYGGGRBGCASCAG \
    --p-r-primer GGACTACHVGGGTWTCTAAT \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads silva-138-ssu-seqs-341f-806r.qza

#dereplicate extracted seqs (shorter sequences may be non-unique again)
qiime rescript dereplicate \
    --i-sequences silva-138-ssu-seqs-341f-806r.qza \
    --i-taxa silva-138-ssu-tax-derep-uniq.qza \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138-ssu-seqs-341f-806r-uniq.qza \
    --o-dereplicated-taxa  silva-138-ssu-tax-341f-806r-derep-uniq.qza
	
# export sequences and taxonomy from qiime
# taxonomy
qiime tools export  --input-path silva-138-ssu-tax-341f-806r-derep-uniq.qza  --output-path silva-138-ssu-tax-341f-806r-derep-uniq
# sequences
qiime tools export  --input-path silva-138-ssu-seqs-341f-806r-uniq.qza --output-path silva-138-ssu-seqs-341f-806r-uniq