#modified from https://gist.github.com/al2na/8540391

DesignProbe<-function(TargetSeq,name,primer1, primer2, startPos, fragLength, optTM,
        primer3="/Applications/primer3-2.3.7/src/primer3_core",
        thermo.param="/Applications/primer3-2.3.7/src/primer3_config/",
        settings="/Users/rpk/GoogleDrive/Kelly_Lab/Bioinformatics/PrimerDesign/primer3_settings.txt"){
#size_range='151-300',Tm=c(57,60,63),
 #/Applications/primer3-2.3.7/primer3web_hybridOligo_settings_rpk.txt
  p3.input="Primer3inputFile.txt"
  p3.output="Primer3outputFile.txt"
  write(
    paste0( sprintf("SEQUENCE_ID=%s\n",name  ),
            sprintf("SEQUENCE_TEMPLATE=%s\n",TargetSeq),
#            sprintf("SEQUENCE_TARGET=%s\n",sequence_target),
#            sequence_exclude_region,
            "PRIMER_TASK=pick_detection_primers\n",
            "PRIMER_PICK_LEFT_PRIMER=0\n" ,
            "PRIMER_PICK_INTERNAL_OLIGO=1\n",
            "PRIMER_PICK_RIGHT_PRIMER=0\n"  ,
            "PRIMER_EXPLAIN_FLAG=1\n"  ,
#            "PRIMER_PAIR_MAX_DIFF_TM=3\n",
#            sprintf("PRIMER_MIN_TM=%s\n" ,Tm[1]),
#            sprintf("PRIMER_OPT_TM=%s\n" ,Tm[2]),
#            sprintf("PRIMER_MAX_TM=%s\n" ,Tm[3]),
#            sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n" ,size_range),
            sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n" ,thermo.param),
			sprintf("SEQUENCE_PRIMER=%s\n" ,primer1),
			sprintf("SEQUENCE_PRIMER_REVCOMP=%s\n" ,primer2),
            sprintf(paste0("SEQUENCE_INCLUDED_REGION=",startPos, ",", fragLength, "\n")),
            sprintf("PRIMER_INTERNAL_OPT_TM =%s\n" ,optTM),
            "=")
    ,  p3.input
  )

try(system2(primer3 ,args=c(sprintf("-p3_settings_file=%s -output=%s",settings, p3.output), p3.input)))

  #import and parse the output into a dataframe named designed.primers
  out=read.delim(p3.output, sep='=', header=FALSE)

  
  returned.primers=as.numeric(as.vector(out[out[,1]=='PRIMER_INTERNAL_NUM_RETURNED',][,2]))
  if (length(returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if ((returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if (returned.primers>0){
    designed.primers=data.frame(matrix(NA, nrow=returned.primers, ncol=9))
    for (i in seq(0,returned.primers-1,1)){
      #IMPORT SEQUENCES
      id=sprintf(  'PRIMER_INTERNAL_%i_SEQUENCE',i)
      PROBE_SEQUENCE=as.character(out[out[,1]==id,][,2])

      #IMPORT PRIMING POSITIONS
      id=sprintf(  'PRIMER_INTERNAL_%i',i)
      POSITION=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #IMPORT Tm
      id=sprintf(  'PRIMER_INTERNAL_%i_TM',i)
      PROBE_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      #IMPORT penalty
      id=sprintf(  'PRIMER_INTERNAL_%i_PENALTY',i)
      PENALTY=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      #IMPORT GC CONTENT
      id=sprintf(  'PRIMER_INTERNAL_%i_GC_PERCENT',i)
      GC=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      #IMPORT SELF ANY
      id=sprintf(  'PRIMER_INTERNAL_%i_SELF_ANY_TH',i)
      SELF=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      #IMPORT HAIRPIN
      id=sprintf(  'PRIMER_INTERNAL_%i_HAIRPIN_TH',i)
      HAIRPIN=as.numeric(as.vector((out[out[,1]==id,][,2])),',')

      #Aggegate in a dataframe
      designed.primers[i+1,] =c(i,
                             PROBE_SEQUENCE,
                             POSITION[1], POSITION[2], PROBE_TM, PENALTY, GC, SELF, HAIRPIN
      					)
    }
	colnames(designed.primers)<-c("Index", "ProbeSequence", "StartPosition", "Length", "Tm", "PositionPenalty", "GC", "SelfPenalty", "Hairpin")
}

return(designed.primers)
}








