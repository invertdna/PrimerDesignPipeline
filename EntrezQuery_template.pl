#!/usr/bin/perl

use LWP::Simple;
$query = '16S[All Fields] AND "Metazoa"[Organism] NOT ("Bacteria"[Organism] OR "Bacteria"[Organism]) NOT "Panulirus argus"[Organism] NOT nuclear[All Fields] NOT "similar to"[All Fields] NOT predicted[All Fields] AND (mitochondrion[filter] AND ("100"[SLEN] : "100000"[SLEN]))';


#assemble the esearch URL
$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y";


#post the esearch URL
$output = get($url);

#parse WebEnv, QueryKey and Count (# records retrieved)
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);


#open output file for writing
open(OUT, ">/Volumes/Nautilus/BlastDB/16Smetazoa/Metazoa16S_20150731.fasta") || die "Can't open file!\n";


#retrieve data in batches of 500
$retmax = 500;
for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
        $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        $efetch_out = get($efetch_url);
        print OUT "$efetch_out";
}
close OUT;

