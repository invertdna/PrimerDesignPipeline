#RPK August 2014; edited July/Aug 2016

#dependencies: 
#ecoPCR https://git.metabarcoding.org/obitools/ecopcr/wikis/home   
#ecoPrimers https://git.metabarcoding.org/obitools/ecoprimers/wikis/home
#NCBI taxdump.tar.gz ftp://ftp.cbi.edu.cn/pub/biomirror/taxonomy/ncbi/
#entrez query template document
#taxize R package   https://cran.r-project.org/web/packages/taxize/index.html

library(taxize)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol   #function from R help file ("integer"); useful for reporting progress in loops


WorkingDir="/Users/rpk/GoogleDrive/Kelly_Lab/Bioinformatics/PrimerDesign/MBON"
setwd(WorkingDir)

projectTitle="CommonMurre"

ecoPrimerspath="/Users/rpk/ecoPrimers/src"
ecoPCRpath="/Users/rpk/ecoPCR/src"
taxdumpPath="/Users/rpk/taxdump"
seqRequest="Alcidae" #what taxonomic group do you want to download sequences for, that will include both ingroup and outgroup?
target_taxon= as.numeric(get_ids("Uria aalge", db="ncbi")$ncbi)  #what taxon are you trying to amplify?
exclude_taxon= as.numeric(get_ids("Alle alle", db="ncbi")$ncbi)  #what taxonomic group are you trying NOT to amplify?
gene="COI" #COI, 16S, etc
nonTargetGene="cytochrome b"

#download relevant dataset from nucleotide db, containing example taxa and counter-example taxa, from genbank
#create and execute perl script for entrez query
##NOTE: because of the way the ecoPrimers quorum and false positive parameters work, it makes sense to specify a particular gene region of interest here rather than to search for all records for a particular taxon
query=paste0("$query = '\"", seqRequest,"\"[organism] NOT Bacteria[organism] AND mitochondrion[filter] AND (\"",gene,"\") NOT (\"", nonTargetGene,"\")';")
ncbidatabase="nucleotide"
format="gb"  #"gb" , "fasta" , etc
outfilename=paste0(projectTitle,"_",gsub(" ","_", gene),".gb")
qfile=readLines("/Users/rpk/GoogleDrive/Kelly_Lab/Bioinformatics/PrimerDesign/EntrezQuery_template.pl")
qfile[4]=query
qfile[9]=paste("$url = $base . \"esearch.fcgi?db=", ncbidatabase,"&term=$query&usehistory=y\";", sep="")
qfile[22]=paste0("open(OUT, \">", WorkingDir,"/",outfilename,"\") || die \"Can't open file!\\n\";", sep="")
qfile[30]=paste("        $efetch_url .= \"&retmax=$retmax&rettype=",format,"&retmode=text\";", sep="")
writeLines(qfile, "EntrezQuery.pl") #create perl script for query
system ("perl EntrezQuery.pl", intern=T)  #run query script and download from ncbi


#housekeeping
infilename=paste(WorkingDir, "/", outfilename, sep="")
foldername=strsplit(basename(infilename), "\\.")[[1]][1]
dir.create(foldername)  #create folder for db files (if necessary) and move gb file into that folder
system(paste("mv", infilename, foldername, sep=" "))
#create logfile
sink(file=paste0(foldername, "/",format(Sys.time(), "%Y_%m_%d_%X"),"_logfile.txt"), type="output", split=T)
#copy this script into new dir, so there's a record
thisScript="/Users/rpk/GoogleDrive/Kelly_Lab/Bioinformatics/PrimerDesign/PrimerDesignPipeline.R" #NOTE: this is hard-coded, because I couldn't figure out a better way to do it.  So you have to change it manually here.
system(paste0("cp ", thisScript, " ", foldername))
#rename script
system(paste0("mv ", foldername,"/", basename(thisScript)," ", foldername, "/",strsplit(basename(thisScript),"\\.")[[1]][1],"_", format(Sys.Date(), "%Y_%m_%d"),"_", projectTitle,".R"))

#format database using ecoPCRFormat
system(paste("cd ", strsplit(ecoPCRpath, "/src")[[1]][1],"/tools;./ecoPCRFormat.py  -g -n ", WorkingDir,"/",foldername, "/", outfilename," -t ",taxdumpPath," ", WorkingDir,"/",foldername, "/", outfilename,"*", sep=""))


#use ecoPrimers to design primers of relevant characteristics; Note odd notation of file location -- ecoPrimers wants the path and then the prefix of its files (e.g., those with extensions .sdx, .ndx, etc)
specificity=.90 #the proportion of the target sequence records that must be good primer matches
quorum=0.95 #the proportion of the sequence records in which a strict match between the primers and their targets occurs (default: 0.7) [not obvious to me what this does if errors_allowed==0...]
falsepositive=0.1 #the maximum proportion of the counterexample sequence records that fulfill the specified parameters for designing the barcodes and the primers
primer_length=22
errors_allowed=2 #note this cuts both ways: for target-group taxa and for exclusion-group taxa, such that allowing 0 errors isn't actually necessarily the way to create the most stringent primers... an exclusion-group taxon could be 1bp away from your primers and you wouldn't know it.
#note option -c considers the circularity of the genome (i.e., for mtDNA primers); I've put this in the call by default.
##note option -3 asks for the min number of perfect nucleotide matches on the 3prime end.  I've set this at 6 by defaut.

ecoPrimersheader=c("serial_number","primer1","primer2","Tm_primer_1","Tm_primer_1_w_mismatches","Tm_primer_2","Tm_primer_2_w_mismatches","C+G_primer_1","C+G_primer_2","good/bad","in_sequence_count","out_sequence_count","yule","in_taxa_count","out_taxa_count","coverage","Number_of_well_identified_taxa","specificity","min_amplified_length","max_amplified_length","avg_amplified_length")
database=paste(WorkingDir,"/",foldername,"/", outfilename, sep="")
EcoPrimersoutfile=paste(WorkingDir,"/",foldername, "/ecoPrimer_results_",format(Sys.time(), "%b_%d_%H:%M:%S"),".txt", sep="")
system(paste("cd ", ecoPrimerspath,";./ecoPrimers -d ",database," -l 100 -L 500 -e ",errors_allowed," -r ", target_taxon," -E ", exclude_taxon," -t species -s ",specificity," -q ",quorum," -x ", falsepositive," -O ", primer_length," -c -3 6 > ", EcoPrimersoutfile, sep=""))   

#read in results, with header, for easy visual inspection  [the default output from ecoPrimers isn't easily human-readable, because it doesn't have headers, but does save the params with which the program was run, which is helpful]
primerResults=read.table(EcoPrimersoutfile)
names(primerResults)=ecoPrimersheader

#FILTER primers for large differences in Tm between primers and for homopolymers
primerResults<-primerResults[abs(primerResults$Tm_primer_1-primerResults$Tm_primer_2)<5,]  #filter out if difference in Tm is greater than 5degrees
primerResults<-primerResults[primerResults$Tm_primer_1>50&primerResults$Tm_primer_1<60,]  #filter by Tm 
primerResults<-primerResults[primerResults$Tm_primer_2>50&primerResults$Tm_primer_2<60,]  #filter by Tm
if(isTRUE(grep("A{4,}|C{4,}|T{4,}|G{4,}", primerResults$primer1))){primerResults<-primerResults[-grep("A{4,}|C{4,}|T{4,}|G{4,}", primerResults$primer1),]} #filter out homopolymers of length greater than 4
if(isTRUE(grep("A{4,}|C{4,}|T{4,}|G{4,}", primerResults$primer2))){primerResults<-primerResults[-grep("A{4,}|C{4,}|T{4,}|G{4,}", primerResults$primer2),]} #filter out homopolymers of length greater than 4
primerResults<-primerResults[order(primerResults[,"good/bad"], decreasing=T),] #sort so Good/Good primers will be first
primerResults<-primerResults[primerResults$avg_amplified_length<300,] #filter for shorter amplicon lengths

#identify unique primer sequences for forward and reverse candidate primers
primerResults<-data.frame(match(primerResults$primer1, unique(primerResults$primer1)), primerResults); names(primerResults)[1:2]<-c("Primer1_ID","Primer2_ID")
primerResults[,2]<-match(primerResults$primer2, unique(primerResults$primer2))
head(primerResults); dim(primerResults)


#use ecoPCR to test for amplification of non-target taxa in database
ecoPCRheader=c("accession_number","sequence_length","taxonomic_id","rank","species_taxonomic_id","scientific_name","genus_taxonomic_id","genus_name","family_taxonomic_id","family_name","super_kingdom_taxonomic_id","super_kingdom_name","strand_(direct_or_reverse)","first_oligonucleotide","number_of_errors_for_the_first_strand","second_oligonucleotide","number_of_errors_for_the_second_strand","amplification_length","sequence_description")
database= database  #using database from above
ecoPCRoutfile=paste0(WorkingDir,"/",foldername,"/ecopcr.out")
#-e Maximum number of errors (mismatches) allowed per primer (default: 0)
#-c Considers that the sequences of the database are circular (e.g. mitochondrial or chloroplast DNA)

#NOTE: customize this script to reflect your design needs; as written, it tries the first bunch of primer sets against the database and stores the number of genera amplified 
primerResults$N_Genera_amplified<-NA 
primerResults$Genera_amplified<-NA
if(nrow(primerResults)>100) rowmax<-100 else rowmax=nrow(primerResults)
for (i in 1: rowmax){
primer1= as.character(primerResults$primer1[i])
primer2= as.character(primerResults$primer2[i])
try(system(paste("cd ",ecoPCRpath,";./ecoPCR -d ",database," -l100 -L500 -e2 -k -c ",primer1," ",primer2," > ", ecoPCRoutfile, sep=""), ignore.stderr = T))
#this is awkward, because ecoPCR writes a file I can't seem to store as an object in R... so I've got to write it out and then read it back in...
	temp=readLines(ecoPCRoutfile); temp=gsub("###","<none>",temp);writeLines(temp, ecoPCRoutfile)
	results=read.table(ecoPCRoutfile, sep="|"); names(results)=ecoPCRheader
primerResults$N_Genera_amplified[i]<-length(unique(results$genus_name))
if(length(unique(results$genus_name))==1) primerResults$Genera_amplified[i]<-as.character(unique(results$genus_name)[1]) else  primerResults$Genera_amplified[i]<-gsub(" +","",paste(unique(results$genus_name), collapse=","))
if(is.wholenumber(i/10)){print(paste0("Trying Primer Set ",i," of ",rowmax))} #report progress every 10th try
}
primerResults$Genera_amplified<-gsub(" +","", primerResults$Genera_amplified)
primerResults_singleGenus<-primerResults[!is.na(primerResults$N_Genera_amplified)&primerResults$N_Genera_amplified==1,] #filter for genus-specific primers, if desired

write.csv(primerResults_singleGenus, paste0(foldername,"/primerResults_filtered_", primer_length,"bp_",gsub(" ","_", gene),"_",format(Sys.time(), "%b_%d_%H:%M:%S"),".csv"), row.names=F)  #write out results

###optional: create fasta file of primer results
# writeLines(paste(">",foldername, "fwd", seq(1:dim(primerResults)[1]), "\r", as.character(primerResults[seq(1:dim(primerResults)[1]),2]),"\r", ">",foldername, "rev", seq(1:dim(primerResults)[1]), "\r", as.character(primerResults[seq(1:dim(primerResults)[1]),3]),"\r",sep=""), paste(WorkingDir,"/",foldername, "/ecoPrimer_results.fasta", sep=""))

sink()


#Double-check for reasonable GC content, 3' end binding affinity, etc.

#THEN, do alignment of a handful of target taxa from ecosystem of interest, and make sure primers that work in theory look like they'll work in practice.  Add ambiguities if absolutely necessary.  Then re-use ecoPCR w degenerate primers

