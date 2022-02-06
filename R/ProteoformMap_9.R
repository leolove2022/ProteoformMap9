
#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
ProteoformMap_9 = function(Peptide_FC_P){
  Peptide_FC_P[,"Fold"]=0 # label all = 0 and add the FC>1.5 as 1 in next
  
  for(n in 1:nrow(Peptide_FC_P)){
    Peptide_FC_P[n,"Protein.Id"]# get protein.Id
    test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]>1.5)&&(Peptide_FC_P[n,"Plasma_p"]<0.05)){# now small means <1
      ###print("OK")
      Peptide_FC_P[n,"selected"]=1
      Peptide_FC_P[n,"Fold"]=1 ## new add to label the fold change for further N,C,M analysis
      ###print(Peptide_FC_P[n,1])
      
    }
    # use to label the FC<1 as -1
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]<1)){# now small means <1
      
      
      Peptide_FC_P[n,"Fold"]=-2 ## new add to label the fold change for further N,C,M analysis
      
      
    }
    # use to label the FC>1 but P>0.05 as -1
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]>1)&&(Peptide_FC_P[n,"Plasma_p"]>0.05)){# now small means <1
      
      
      Peptide_FC_P[n,"Fold"]=-1 ## new add to label the fold change for further N,C,M analysis
      
      
    }
    ###print(n)
  }
  # add all contain 1 row
  # the Sequence may change, 
  for(n in 1:nrow(Peptide_FC_P)){
    { 
      test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
      if((max(test$selected)==1)){# not know why have a ".1"
        
        Peptide_FC_P[n,"selected"]=1
        ##print(Peptide_FC_P[n,"Protein.Id"])
      }
      
    }
  }
  
  for(n in 1:nrow(Peptide_FC_P)){
    Peptide_FC_P[n,"Protein.Id"]# get protein.Id
    test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]>1.5)&&(Peptide_FC_P[n,"Plasma_p"]<0.05)){# now small means <1
      ###print("OK")
      Peptide_FC_P[n,"selected"]=1
      #Peptide_FC_P[n,"Fold"]=0 ## new add to label the fold change for further N,C,M analysis
      ###print(Peptide_FC_P[n,1])
      
    }
    ###print(n)
  }
  # add all contain 1 row
  # the Sequence may change, 
  for(n in 1:nrow(Peptide_FC_P)){
    { 
      test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
      if((max(test$selected)==1)){# not know why have a ".1"
        
        Peptide_FC_P[n,"selected"]=1
        ###print(Peptide_FC_P[n,"Protein.Id"])
      }
      
    }
  }
  
  # run above to label all contain select =1 proteins, the other peptide is alo =1
  #address_Threshold_peptide=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_complete_1_all_thresthod0.csv")
  #address_Threshold_peptide_2=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_completed_Selected_thresthod0.csv")
  #write.csv(Peptide_FC_P,file=address_Threshold_peptide,row.names=FALSE)# transform ok
  Peptide_FC_P_Selected=subset(Peptide_FC_P,selected==1) ############### This will use in deduplicate
  #write.csv(Peptide_FC_P_Selected,file=address_Threshold_peptide_2,row.names=FALSE)# transform ok
  # get the selected peptide result
  # get Peptide_FC_P_Selected
  
  #### deduplicate to get the list of query protein sequence
  Peptide_duplicate_2<-Peptide_FC_P_Selected[!duplicated(Peptide_FC_P_Selected$Protein.Id),] #### this is the deduplicate list use to get protein seq
  #address2="~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_complete_1_selected_threstho0_manual_deduplicate.csv"
  #write.csv(Peptide_duplicate_2,file=address2,row.names=FALSE)# transform ok
  # produce a no redundent Protein.Id list to get the protein sequence
  
  # set the query list
  #Peptide_query <- data.frame(Protein.Id= character(), age= numeric(), stringsAsFactors=FALSE)
  Peptide_query <- data.frame(Protein.Id= Peptide_duplicate_2$Protein.Id)
  
  ############get the Peptide Query to get the qequence
  #up <- UniProt.ws(taxId=9606) # human
  #species(up)
  mouseUp <- UniProt.ws(10090) # mouse
  #mouseUp
  
  
  #address_ProteinUniprotID=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/second analysis/Unipro_test.csv")
  #address_ProteinSeq=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/second analysis/Uniprot_ProteinSeq.csv")
  
  #data_UniprotAccess<-read.table(address_ProteinUniprotID,header=TRUE,sep=',')
  keys <- Peptide_query_test_p$Protein.Id
  columns <- c("PROTEIN-NAMES","SEQUENCE")
  kt <- "UNIPROTKB"
  res <- AnnotationDbi::select(mouseUp, keys, columns, kt)
  Peptide_ProteinSeq=res
  
  Peptide_ProteinSeq_2<-na.omit(Peptide_ProteinSeq) # delete the NA value ##########This is the Peptide_ProteinSeq_2 is the protein seq
  #write.csv(res,file=address_ProteinSeq,row.names=FALSE)# transform ok
  ##
  ##### begin to labe the position
  
  indexOf = function(str,str2){
    cd=nchar(str);
    cd2=nchar(str2);
    if(cd==0||cd2==0){
      return(0);
    }
    for(i in 1:cd){
      t=substr(str,i,i);
      for(j in 1:cd2){
        if(t==substr(str2,j,j)&&j==1){
          if(cd2==1){
            return(i);
          }else{
            c=TRUE;
            for(k in 1:(cd2-1)){
              if(substr(str,i+k,i+k)!=substr(str2,j+k,j+k)){
                c=FALSE;
                break;
              }
            }
            if(c==TRUE){
              return(i);
            }
          }
        }else{
          break;
        }
      }
    }
    return(0);
  }
  
  ## wait to align
  # data_peptide= Peptide_FC_P_Selected, data_protein_seq=Peptide_ProteinSeq_2
  # example
  #address_1=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_complete_1_selected_threstho0_manual_align.xlsx")
  #address_2=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_complete_1_selected_threstho0_manual_deduplicate_seq.xlsx")
  
  #address_4=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_complete_1_selected_threstho0_manual_deduplicate_align_info_p.csv")
  # example
  #data_peptide<-read.xlsx(address_1,rowNames = F,colNames = T)
  data_peptide<-Peptide_FC_P_Selected
  #data_protein_seq<-read.xlsx(address_2,rowNames = F,colNames = T)
  data_protein_seq<-Peptide_ProteinSeq_2
  for(n in 1: nrow(data_peptide)){# begin 0
    
    #n=1
    data_peptide[n,"Protein.Id"]
    data_protein_seq_test<-subset(data_protein_seq,UNIPROTKB==data_peptide[n,"Protein.Id"])
    if(!is.na(data_protein_seq_test[1,"SEQUENCE"])&&nchar(data_protein_seq_test[1,"SEQUENCE"])>10){# begin1
      a=data_protein_seq_test[1,"SEQUENCE"]
      data_peptide[n,"PeptideLocation"]=indexOf(a,data_peptide[n,"PeptideSeq"])
      data_peptide[n,"PeptideLen"]=nchar(data_peptide[n,"PeptideSeq"])
      data_peptide[n,"ProteinLen"]=nchar(data_protein_seq_test[1,"SEQUENCE"])
      data_peptide[n,"RTPosition"]=indexOf(a,data_peptide[n,"PeptideSeq"])/nchar(data_protein_seq_test[1,"SEQUENCE"])
      ##print(data_peptide[n,1])
      ##print(n)
      
    }# final1
  }# final 0
  #write.csv(data_peptide,file=address_4,row.names=FALSE)
  
  # ########analysis the RTPosition and the Fold use the data_peptide from above section
  for(n in 1: nrow(data_peptide)){# begin 0
    
    #n=1
    data_peptide[n,"Protein.Id"]
    data_protein_seq_test<-subset(data_protein_seq,UNIPROTKB==data_peptide[n,"Protein.Id"])
    if(!is.na(data_protein_seq_test[1,"SEQUENCE"])&&nchar(data_protein_seq_test[1,"SEQUENCE"])>10){# begin1
      a=data_protein_seq_test[1,"SEQUENCE"]
      data_peptide[n,"PeptideLocation"]=indexOf(a,data_peptide[n,"PeptideSeq"])
      data_peptide[n,"PeptideLen"]=nchar(data_peptide[n,"PeptideSeq"])
      data_peptide[n,"ProteinLen"]=nchar(data_protein_seq_test[1,"SEQUENCE"])
      data_peptide[n,"RTPosition"]=indexOf(a,data_peptide[n,"PeptideSeq"])/nchar(data_protein_seq_test[1,"SEQUENCE"])
      #print(data_peptide[n,1])
     # print(n)
      
    }# final1
  }# final 0
  #write.csv(data_peptide,file=address_4,row.names=FALSE)
  
  # ########analysis the RTPosition and the Fold use the data_peptide from above section
  
  # 20220206 version 3 
  for(n in 1:nrow(data_protein_seq)){
    #print(n)
    
    #n=2 # test # This method may meet mistake wen RT position more than 1
    test2<-subset(data_peptide,Protein.Id==data_protein_seq[n,"UNIPROTKB"])
    #print(test2[1,1])
    firstpeptideposition=min(test2$RTPosition)
    lastpeptideposition=max(test2$RTPosition)
    test2_first<-subset(test2,RTPosition==firstpeptideposition)
    test2_last<-subset(test2,RTPosition==lastpeptideposition)
    test2_N=test2_first[1,"Fold"]
    test2_C=test2_last[nrow(test2_last),"Fold"]
    #print(test2[1,1])
    # This method may meet mistake wen RT position more than 1
    # 1 C case
    
    if((test2_N<0)&&(test2_C>0)){
      print(n)####this print not close
      test2_small<-subset(test2,Fold==-2)
      test2_large<-subset(test2,Fold==1)
      if(min(test2_large$RTPosition)>max(test2_small$RTPosition)){
        for(n in 1: nrow(data_peptide)){
          if(data_peptide[n,"Protein.Id"]==test2[1,1]){
            data_peptide[n,"Proteoform"]="C"
          }# all change as C
          
        }# This is for to change the all peptide select as C
        
        
        
      }# this is the test small and test large compare
      
    }# this is the for the list in proteinseq list all name
    #2 N case
    if((test2_N>0)&&(test2_C<0)){
      test2_small<-subset(test2,Fold==-2)
      test2_large<-subset(test2,Fold==1)
      if(min(test2_small$RTPosition)>max(test2_large$RTPosition)){
        for(n in 1: nrow(data_peptide)){# all change as C
          if(data_peptide[n,"Protein.Id"]==test2[1,1]){
            data_peptide[n,"Proteoform"]="N"
          }# all change as C
          
        }# This is for to change the all peptide select as C
      }# this is the test small and test large compare
    }# this is the for the list in proteinseq list all name
    
    ##############
    # test use talbe for M case
    
    # 3 M case   method 3
    # test for fga wrong
    #test2<-subset(Peptide_proteoform,Protein.Id=="E9PV24")
    firstpeptideposition=min(test2$RTPosition)
    lastpeptideposition=max(test2$RTPosition)
    test2_first<-subset(test2,RTPosition==firstpeptideposition)
    test2_last<-subset(test2,RTPosition==lastpeptideposition)
    
    if((test2_N<0)&&(test2_C<0)){
      proteoM="M"
      # first order the test2
      #TEST4<-test2[order(test2[,"RTPosition"]),test2[,"F"]),]
      test2_M_large=subset(test2,Fold>0)
      test2_M_small=subset(test2,Fold<0)
      test2_M_smallRTP<-min(test2_M_large$RTPosition)
      test2_M_largeRTP<-max(test2_M_large$RTPosition)
      for(n in 1: nrow(test2_M_small)){
        test2_M_samll_value=test2_M_small[n,"RTPosition"]
        if((test2_M_samll_value>test2_M_smallRTP)&&(test2_M_samll_value<test2_M_largeRTP)){
          proteoM=NA
          break
        }
        #print(proteoM)
        #print(n)
      }
      for(n in 1: nrow(data_peptide)){# all change as C
        if(data_peptide[n,"Protein.Id"]==test2[1,1]){
          data_peptide[n,"Proteoform"]=proteoM
          
        }# all change as C
      }
      #test2$Proteoform=proteoM # test
      ############# 3 M case   method 3
    }# this is for for proteoM if end } not use for test
    
    
    ##    test2_small<-subset(test2,Fold==-2)
    ##    test2_large<-subset(test2,Fold==1)
    ##    if(min(test2_small$RTPosition)>max(test2_large$RTPosition)){
    # #     for(n in 1: nrow(data_peptide)){
    ##       if(data_peptide[n,"Protein.Id"]==test2[1,1]){
    ##          data_peptide[n,"Proteoform"]="N"
    ##       }# all change as C
    ##      }# This is for to change the all peptide select as C
    ##      }# this is the test small and test large compare
    
  }# this is the for the list in proteinseq list all name This is the for cycle for the C,N,M adugement
  #above is judge N and C 
  #####################  20220206 version 3  change test2 adjuget for M, set the firt and last peptide and change the TEST4, and delete the wrong test2 value
  
  # output the preotorm not NA from data_peptide above
  
  Peptide_proteoform=data_peptide[complete.cases(data_peptide[,"Proteoform"]),]
  Peptide_proteoform_list=Peptide_proteoform[!duplicated(Peptide_proteoform$Protein.Id),]
  
  
  # final1
  # final 0
  return(Peptide_proteoform);
}


### test proteoformmap8
