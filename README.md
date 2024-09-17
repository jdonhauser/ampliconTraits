# ampliconTraits

ampliconTraits is a trait sequence database and a workflow to classify environmental ASVs. It uses a phenotypic trait database (Madin et al. 2020) with sequences from the SILVA SSU Ref v138 database. Trait databases are available for the V3-V4 region (341F - 806R) from the prokaryotic 16S rRNA gene. ampliconTraits uses SINAPS (Edgar 2017) to classify environmental ASVs. See [Make your own databases](#make-your-own-databases)
for creation of customized trait sequence databases, e.g databases for a different region of the 16S rRNA gene. Currently available traits are: d1_lo (cell diameter, lower), d1_up (cell diameter, upper) , d2_lo (cell length, lower), d2_up (cell length, upper), doubling_h (doubling time), genome_size, optimum_pH, optimum_tmp, rRNA16S_genes (16S rRNA gene copy number), metabolism (oxygen preference), motility, salinity_range (salinity preference) and sporulation. 

Check out the paper: https://doi.org/10.1016/j.ecoinf.2024.102817

For questions please contact: jdonhauser93@gmail.com

## Get the databases
Download the databases

```
wget https://erda.ku.dk/archives/f5d4b1d41f74ba3d6f73b212dbb11591/traitDatabases.tar
```
The trait databases are called `trait.fasta`. For continuous traits, annotations have been binned into discrete intervals to enable classification. The number after the trait name in the file name indicates the number of equally sized intervals, the trait was binned into. E.g. `genome_size10.fasta` contains genome size annotations, obtained from binning the whole range of genome sizes into 10 intervals. For 16S rRNA gene copy number, moreover exact copy numbers are available. We tested interval numbers between 5 and 30 for most traits, resulting in similar accuracy. The choice of interval number is a trade-off between the precision and the confidence of a classification. If the number of intervals is too low, differences in the trait composition of environmental samples may not be resolved. The best choice may depend on the dataset. From our experience we recommend:
* d1_lo: 20 - 30 intervals
* d1_up: 20 - 30 intervals
* d2_lo: 30 intervals
* d2_lo: 30 intervals
* doubling_h: 20 - 50 intervals
* genome_size: 10 - 20 intervals 
* optimum_pH: 10 - 20 intervals 
* optimum_tmp: 10 - 20 intervals 
* rRNA16S_genes: exact copy numbers

## Trait classifications of ASVs
Classifying environmental ASVs with SINAPS requires the installation of usearch as described here: https://www.drive5.com/usearch/manual/install.html. 
The SINAPS command and algorithm are described in more detail here: https://drive5.com/usearch/manual/cmd_sinaps.html and here: https://drive5.com/usearch/manual/sinaps_algo.html 

When using SINAPS please cite:
Edgar, R. C. (2017). SINAPS: Prediction of microbial traits from marker gene sequences (p. 124156). bioRxiv. https://doi.org/10.1101/124156

To classify ASVs from your dataset for e.g. genome size you need the sequences for the ASVs in a fasta file, e.g. `ASVs.fasta`
```
>000274f8057c3721b573e5175f08ad47
GCGCGAAACCTTTACACTGCACGACAGTGCGATAAGGGGACTCCGAGTGCGAGGACATACTAGTCCTCGCTTTTACCGACCGTAAGGTGG
>00085ae998183f5005393eed40a1adb6
TCGAGGATCTTCGGCAATGGACGCAAGTCTGACCGAGCGACGCCGCGTGCGGGATGAAGGCCTTCGGGTTGTAAACCGCTGTCAGTGGGG
>000842d488485e1241ee02c02e29b6cf
TAGGGAATTTTCCACAATGGGCGAAAGCCTGATGGAGCAACGCCGCGTGCAGGATGAAGGCCTTCGGGTTGTAAACTGCTTTTATGTATG
```
and the database for genome size `genome_size10.fasta`. The headers of the sequences indicate sequenceID;trait=value. For continuous traits like genome size, value is a range.
```
>AB001777.1.1508;genome_size10=1.19e+05-1.73e+06
TCGAGAATCTTTCGCAATGGACGAAAGTCTGACGAAGCGACGCCGCGTGTGTGATGAAGG
>L33689.1.1458;genome_size10=4.91e+06-6.5e+06
TGGGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATGCCGCGTGAGTGAAGAAGG
>AB001815.1.1507;genome_size10=1.19e+05-1.73e+06
TCGAGAATCTTTCGCAATGGACGAAAGTCTGACGAAGCGACGCCGCGTGTGTGATGAAGG
```
Then use the following command to classify your sequences

```
usearch -sinaps ASVs.fasta -db genome_size10.fasta -attr genome_size10 -tabbedout genomeclassification.txt -strand plus
```
| Option | Description |
| --- | --- |
| -db | sequence trait database |
| -attr | name of the trait |
| -tabbedout | output file |

The columns of the output `genomeclassification.txt` are query sequence id; annotation; bootstrap value; strand
```
000274f8057c3721b573e5175f08ad47	3.32e+06-4.91e+06	77	+
00085ae998183f5005393eed40a1adb6	6.5e+06-8.09e+06	49	+
000842d488485e1241ee02c02e29b6cf	1.73e+06-3.32e+06	49	+
```
Note that the left border of the lowest interval were slightly decreased and the right border of the highest interval was slightly increased to include the minimum and maximum value of the trait, respectively, when creating the databases (this explains why some values become negative). Use the file `IntervalLimits_Reconversion.csv` to convert annotations corresponding to the lowest or highest interval back to their original value. Columns are "lowest interval original value" ,"lowest interval new value", "highest interval original value", "highest interval new value"
```
		LI_or		LI_cut		HI_or		HI_cut
optimum_tmp10	3-13.2		2.9-13.2	94.8-105	94.8-105
optimum_tmp20	3-8.1		2.9-8.1		99.9-105	99.9-105
d1_lo5		0.1-14.8	-0.65-14.8	735-750		735-751
```
For instance in R this could be done with something like

```
genomeclassification[genomeclassification$Value==Lim["genome_size10",'LI_cut'],'Value'] <- Lim["genome_size10",'LI_or'] # reconvert lowest interval
genomeclassification[genomeclassification$Value==Lim["genome_size10",'HI_cut'],'Value'] <- Lim["genome_size10",'HI_or'] # reconvert highest interval
```
where `genomeclassification` is a data frame with the SINAPS output where the column with the trait annotation is called `Value` and
`Lim` is a data frame with the conversion table


## Make your own databases
The scripts to make the databases (a combination of R and command line) are in `createDatabases` and were run interactively. The following sections give some explanation. 
The scripts expect the following direcory structure:
```
base_directory/
├── 1_prepareSilva
├── 2_matchTaxids
└── 3_makeTraitSpecificDB
```
Here is an overview of the files used and produced. The files `condensed_species_NCBI.csv` and `taxmap_embl-ebi_ena_ssu_ref_138.txt` need to be downloaded, all other files are created by the scripts.
```
  base_directory/
├── 1_prepareSilva
│   ├── silva-138-ssu-rna-seqs.qza
│   ├── silva-138-ssu-seqs-341f-806r-uniq
│   │   └── dna-sequences.fasta
│   ├── silva-138-ssu-seqs-341f-806r-uniq.qza
│   ├── silva-138-ssu-seqs-341f-806r.qza
│   ├── silva-138-ssu-seqs-cleaned.qza
│   ├── silva-138-ssu-seqs-derep-uniq.qza
│   ├── silva-138-ssu-seqs-discard.qza
│   ├── silva-138-ssu-seqs-filt.qza
│   ├── silva-138-ssu-seqs.qza
│   ├── silva-138-ssu-tax-341f-806r-derep-uniq
│   │   └── taxonomy.tsv
│   ├── silva-138-ssu-tax-341f-806r-derep-uniq.qza
│   ├── silva-138-ssu-tax-derep-uniq.qza
│   └── silva-138-ssu-tax.qza
├── 2_matchTaxids
│   ├── condensed_species_NCBI.csv
│   ├── dbMadin_merged_Silva.csv
│   ├── dna-sequences-renamed.fasta
│   ├── id_list.txt
│   ├── numberSpeciesAndSequences_perTrait.csv
│   ├── silva138_SSURef_inMadinDB.fasta
│   ├── taxaMatchedSpeciesNameButNotTaxid.csv
│   └── taxmap_embl-ebi_ena_ssu_ref_138.txt
└── 3_makeTraitSpecificDB
    ├── trait1 
    │   ├── id_list.txt
    │   ├── seq_headers
    │   ├── silva138_SSURef_inMadinDB_trait1.fasta
    │   ├── silva138_SSURef_inMadinDB_trait1_renamed.fasta
    │   └── traits.tsv
    └── trait2 
        ├── id_list.txt
        ├── seq_headers
        ├── silva138_SSURef_inMadinDB_trait2.fasta
        ├── silva138_SSURef_inMadinDB_trait2_renamed.fasta
        └── traits.tsv
```


### Download SILVA SSURef_138 and make amplicon specific version of the database
This part downloads sequences from SILVA to be used for trait classifications and creates an amplicon specific version of the database. We created a database for the 341F-806R of the V3-V4 region of the 16S rRNA gene. If you use a different amplicon, change accordingly. Opposed to most applications for taxonomy classification we did NOT use the clustered NR99 version of the database. Using all sequences allows to cross-map more sequences to the trait database. The
directory for these steps is `1_prepareSilva`. Execute all steps in the script 1_silva138ref_341F806R_build.sh. The scripts uses RESCRIPt to download SILVA and extract the region of our amplicon. Instructions for the installation of Rescript can be found here: https://github.com/bokulich-lab/RESCRIPt/. The script follows the workflow described here: https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494

When using SILVA please cite

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J,
Glockner FO (2013) The SILVA ribosomal RNA gene database project: improved
data processing and web-based tools. Nucleic Acids Research 41:D590-D596

and when using RESCRIPt cite

Ii, M. S. R., O’Rourke, D. R., Kaehler, B. D., Ziemski, M., Dillon, M. R., Foster, J. T., & Bokulich, N. A. (2021). RESCRIPt: Reproducible sequence taxonomy reference database management. PLOS Computational Biology, 17(11), e1009581. https://doi.org/10.1371/journal.pcbi.1009581

The end product is an amplicon specific fasta file (in our case for 341F-806R of the V3-V4 region): 
```
>AB000106.1.1343 Bacteria;Proteobacteria;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae;Sphingobium;Sphingomonas sp.
TAGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCGATTTAAGTCAGAGGTGAAAGCCCGGGGCTCAACCCCGGAATAGCCTTTGAGACTGGATTGCTTGAATCCGGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGATCACTGGACCGGCATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG
>HL182401.2.1459 Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;unidentified
TAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTTTTCGGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTGCCGTTCGAATAGGGCGGCACCTTGACGGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGGTTTCTTTAAGTCTGATGTGAAAGCCCCCGGGCTCAACCGGGGAGGGGTCATTGGAAACTGGGGAAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTTGTAGCCGGGTGAAATGCGTAGAGGATGTGGAGGAACACCAGTGGGGCGAAGGCGACTCTCTGGGTCTGTAACTGACGCTGAGGGAGCGAAAGCGTGGGGGAGCGAACAGG
```
and a corresponding taxonomy file: 
```
AB000106.1.1343	d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Sphingobium; s__Sphingomonas_sp.
HL182401.2.1459	d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus; s__unidentified
```


### Match taxids from the trait database with SILVA sequences
The second part matches species from the trait database with species from SILVA using the ncbi taxid. The working directory for this part is `2_matchTaxids`. You need the trait database file `condensed_species_NCBI.csv`. This file contains the ncbi taxid, the taxonomy and trait annotations:
```
  species_tax_id                         species                  genus               family               order                 class         phylum superkingdom gram_stain      metabolism pathways carbon_substrates sporulation motility
1        1243001 Acidipropionibacterium damnosum Acidipropionibacterium Propionibacteriaceae Propionibacteriales        Actinobacteria Actinobacteria     Bacteria       <NA> microaerophilic     <NA>              <NA>        <NA>     <NA>
2        1679466            Apibacter adventoris              Apibacter    Flavobacteriaceae    Flavobacteriales        Flavobacteriia  Bacteroidetes     Bacteria       <NA> microaerophilic     <NA>              <NA>        <NA>     <NA>
3        1591092              Aquaspirillum soli          Aquaspirillum   Chromobacteriaceae        Neisseriales    Betaproteobacteria Proteobacteria     Bacteria       <NA> microaerophilic     <NA>              <NA>        <NA>     <NA>
```
when using the trait database please cite:

Madin, J. S., Nielsen, D. A., Brbic, M., Corkrey, R., Danko, D., Edwards, K., Engqvist, M. K. M., Fierer, N., Geoghegan, J. L., Gillings, M., Kyrpides, N. C., Litchman, E., Mason, C. E., Moore, L., Nielsen, S. L., Paulsen, I. T., Price, N. D., Reddy, T. B. K., Richards, M. A., … Westoby, M. (2020). A synthesis of bacterial and archaeal phenotypic trait data. Scientific Data, 7(1), Article 1. https://doi.org/10.1038/s41597-020-0497-4

This file was published under the Creative Commons Attribution 4.0 (CC-BY 4.0) license (https://creativecommons.org/licenses/by/4.0/) and can be downloaded here under the same license:

```
wget https://erda.ku.dk/archives/f5d4b1d41f74ba3d6f73b212dbb11591/ampliconTraits/condensed_species_NCBI.csv.gz
wget https://erda.ku.dk/archives/f5d4b1d41f74ba3d6f73b212dbb11591/ampliconTraits/condensed_species_NCBI.csv.gz.md5
```

Moreover, you need the file `taxmap_embl-ebi_ena_ssu_ref_138.txt` allowing to link the taxid with the accession number in the headers of the fasta file from SILVA. 
```
primaryAccession	start	stop	submitted_path	submitted_name	ncbi_taxonid
AB000106	1	1343	Bacteria;Proteobacteria;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae;Sphingomonas;	Sphingomonas sp.	28214
AB000278	1	1410	Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Photobacterium;	Photobacterium iliopiscarium	56192
AB000389	1	1508	Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Pseudoalteromonadaceae;Pseudoalteromonas;	Pseudoalteromonas elyakovii	81037
```
This file was downloaded from the SILVA web site where it was published under the Creative Commons Attribution 4.0 (CC-BY 4.0) license (https://creativecommons.org/licenses/by/4.0/). It can be downloaded here under the same license:

```
wget https://erda.ku.dk/archives/f5d4b1d41f74ba3d6f73b212dbb11591/ampliconTraits/taxmap_embl-ebi_ena_ssu_ref_138.txt.gz
wget https://erda.ku.dk/archives/f5d4b1d41f74ba3d6f73b212dbb11591/ampliconTraits/taxmap_embl-ebi_ena_ssu_ref_138.txt.gz.md5
```

Then, in R follow the steps in `2.1_silva2ncbi.R`. This script creates a table `dbMadin_merged_Silva.csv` with sequence IDs corresponding to the species in the trait database. The column Feature.ID corresponds to the sequence identifier and is exported as a separate file `id_list.txt` to select corresponding sequences in the fasta files in the next script. 
```
	species_tax_id	species.Madin	genus.Madin	family.Madin	order.Madin	class.Madin	phylum.Madin	superkingdom	gram_stain	metabolism	pathways	carbon_substrates	etc.	data_source	ref_id	Feature.ID	kingdom	phylum.silva	class.silva	order.silva	family.silva	genus.silva	species.silva	primaryAccession.x	primaryAccession.y	start	stop	submitted_path	submitted_name	ncbi_taxonid
1	1243001	Acidipropionibacterium damnosum	Acidipropionibacterium	Propionibacteriaceae	Propionibacteriales	Actinobacteria	Actinobacteria	Bacteria	NA	microaerophilic	NA	NA	NA	bacdive-microa	151	JQ283461.1.1380	d__Bacteria	 p__Actinobacteriota	 c__Actinobacteria	 o__Propionibacteriales	 f__Propionibacteriaceae	 g__Acidipropionibacterium	Acidipropionibacterium damnosum	JQ283461	JQ283461	1	1380	Bacteria;Actinobacteria;Propionibacteriales;Propionibacteriaceae;Acidipropionibacterium;	Acidipropionibacterium damnosum	1243001
2	1679466	Apibacter adventoris	Apibacter	Flavobacteriaceae	Flavobacteriales	Flavobacteriia	Bacteroidetes	Bacteria	NA	microaerophilic	NA	NA	NA	bacdive-microa	151	KT149222.1.1529	d__Bacteria	 p__Bacteroidota	 c__Bacteroidia	 o__Flavobacteriales	 f__Weeksellaceae	 g__Apibacter	Apibacter adventoris	KT149222	KT149222	1	1529	Bacteria;Bacteroidetes;Flavobacteriia;Flavobacteriales;Flavobacteriaceae;Apibacter;	Apibacter adventoris	1679466
3	1679466	Apibacter adventoris	Apibacter	Flavobacteriaceae	Flavobacteriales	Flavobacteriia	Bacteroidetes	Bacteria	NA	microaerophilic	NA	NA	NA	bacdive-microa	151	KT149221.1.1528	d__Bacteria	 p__Bacteroidota	 c__Bacteroidia	 o__Flavobacteriales	 f__Weeksellaceae	 g__Apibacter	Apibacter adventoris	KT149221	KT149221	1	1528	Bacteria;Bacteroidetes;Flavobacteriia;Flavobacteriales;Flavobacteriaceae;Apibacter;	Apibacter adventoris	1679466
```
Essentially the steps are:
* combine taxonomy from previous part with `taxmap_embl-ebi_ena_ssu_ref_138.txt` to obtain the ncbi taxid for each sequence
* merge with trait database based on taxid (majority of species)
* cross-map additional species by species name (in some cases the taxid corresponds to a strain in one of the resources and therefore does not match the species taxid)

Then with the script `subsetSILVA.sh` subset the fasta file with the SILVA sequences to those the species in the trait database. This script requires installation of SeqKit (https://bioinf.shenwei.me/seqkit/). When using seqkit please cite:

Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), e0163962. https://doi.org/10.1371/journal.pone.0163962

The output file is `silva138_SSURef_inMadinDB.fasta`

### Make trait specific databases and format for SINAPS
Not all traits are available for all species. This part subsets the database for each trait and adds the trait and its value to the header of the fasta. Continuous traits moreover have to be binned in discrete intervals, so the trait value becomes a range. The working directory for this part is `3_makeTraitSpecificDB`

Use the script `makeTraitSpecificDatabase.R` to create a list of sequences available for each trait and create trait values. For categorical values the trait values correspond to the entry in the database, for continuous traits the value is a range. The R function `cut` divides values in x bins of equal size (which may result in bins without values if there are gaps in the data). We selected x to obtain the desired number of non-empty bins. Note that the function slightly decreases the left border of the lowest interval and slightly increases the right border of the highest interval to include the lowest and highest value in the interval. The file `IntervalLimits_Reconversion.csv` can be used to reconvert trait annotations to the original value. For each trait a directory is created. Number suffixes for continuous traits indicate the number of intervals, e.g. genome_size10 and genome_size20 indicate genome size split in 10 and 20 intervals, respectively. For each trait-interval combination the script creates the file `traits.tsv` containing the sequence identifier and the corresponding trait value, e.g.
```
CP029206.393659.395211	1.73e+06-3.32e+06
EF173408.1.1387	1.73e+06-3.32e+06
JF311443.1.1345	3.32e+06-4.91e+06
```
Use the script `makeTraitSpecificDatabase_SINAPS.sh` to subset the SILVA sequences to the sequences available for each trait and write the trait annotation to the fasta header as required for classification with SINAPS.

The final output is a fasta file for each trait with trait annotation in the values `silva138_SSURef_inMadinDB_trait_renamed.fasta`, ready to be used with SINAPS, e.g. for genome size
```
>AB001777.1.1508;genome_size10=1.19e+05-1.73e+06
TCGAGAATCTTTCGCAATGGACGAAAGTCTGACGAAGCGACGCCGCGTGTGTGATGAAGG
```

## References
Donhauser, J., Doménech-Pascual, A., Han, X., Jordaan, K., Ramond, J.-B., Frossard, A., Romaní, A.M., Priemé, A., 2024. Modelling soil prokaryotic traits across environments with the trait sequence database *ampliconTraits* and the R package *MicEnvMod*. Ecological Informatics 83, 102817. doi:10.1016/j.ecoinf.2024.102817

Edgar, R. C. (2017). SINAPS: Prediction of microbial traits from marker gene sequences (p. 124156). bioRxiv. https://doi.org/10.1101/124156

Ii, M. S. R., O’Rourke, D. R., Kaehler, B. D., Ziemski, M., Dillon, M. R., Foster, J. T., & Bokulich, N. A. (2021). RESCRIPt: Reproducible sequence taxonomy reference database management. PLOS Computational Biology, 17(11), e1009581. https://doi.org/10.1371/journal.pcbi.1009581

Madin, J. S., Nielsen, D. A., Brbic, M., Corkrey, R., Danko, D., Edwards, K., Engqvist, M. K. M., Fierer, N., Geoghegan, J. L., Gillings, M., Kyrpides, N. C., Litchman, E., Mason, C. E., Moore, L., Nielsen, S. L., Paulsen, I. T., Price, N. D., Reddy, T. B. K., Richards, M. A., … Westoby, M. (2020). A synthesis of bacterial and archaeal phenotypic trait data. Scientific Data, 7(1), Article 1. https://doi.org/10.1038/s41597-020-0497-4

Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research, 41(Database issue), D590–D596. https://doi.org/10.1093/nar/gks1219


