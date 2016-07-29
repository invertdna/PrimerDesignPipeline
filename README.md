This R script functions as a wrapper to automate some of the functions of Coissac et al.'s ecoPCR and ecoPrimers software (see http://metabarcoding.org/obitools/doc/scripts/ecoPrimers.html). The idea is that you should be able to iteratively design and test (in-silico) your primers in a rapid way using the R interface.

The script works by using R to call ecoPrimers and related software, asking the user to define key variables and paths along the way. It also edits a perl script template to get reference sequences from NCBI via entrez.

The output is a csv file of candidate primers meeting user design criteria. The script calls ecoPCR at the end to see if the candidate primers (by default, the first set, but this can be changed) are likely to amplify the target taxa out of the downloaded genbank dataset. 



#dependencies: 
#ecoPCR https://git.metabarcoding.org/obitools/ecopcr/wikis/home   
#ecoPrimers https://git.metabarcoding.org/obitools/ecoprimers/wikis/home
#NCBI taxdump.tar.gz ftp://ftp.cbi.edu.cn/pub/biomirror/taxonomy/ncbi/
#entrez query template document 
#taxize R package   https://cran.r-project.org/web/packages/taxize/index.html
