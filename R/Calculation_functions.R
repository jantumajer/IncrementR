###########################################
### Function for calculation of eccentricity index according to Schweingruber (1996), Braam et al. (1987) and Alestalo et al. (1971)
### For details about eqations used for calculation of indexes, see notes below or summary in Tumajer, Treml (2013), Geochronometria 40(1): 59-76.
### Procedure: [1] creates new output tabs (private functions) [2] splits table with RWL series into 4 tables according aspect (data taken from .IDformat$Aspect - private function) [3] calculates indexes and stores them in three separate tables -> list
### Stores a results (i.e., eccentricity indexes) in RWL (in R) format - other functions of eg., dplR can be used in further calculations
### Title of the eccentricity series name: IDPlot_IDTree_IDLevelOrientation, where Orientation means Direction from which eccentricity index is calculated
###########################################

Eccentricity <- function (trw.series, complete=FALSE)
{
  ID<-.IDdistinct(trw.series, complete)

  tab.ecc.schw <-.CreateTableTRW(trw.series,ID$IDAspect)
  colnames.store.schw <- .CreateTableMeta(ID$IDAspect)

  tab.ecc.braam <-.CreateTableTRW(trw.series,ID$IDAspect)
  colnames.store.braam <- .CreateTableMeta(ID$IDAspect)

  tab.ecc.alestalo <-.CreateTableTRW(trw.series,ID$IDAspect)
  colnames.store.alestalo <- .CreateTableMeta(ID$IDAspect)
  ###################################################################################################

  for (i in 1:nrow(ID$IDAspect))
  {

    ser.n <- as.character(ID$IDAspect[i,"N.OC"]) # Extracts names of series from the same height level from IDformat$ID Aspect
    ser.s <- as.character(ID$IDAspect[i,"S.OC"])
    ser.w <- as.character(ID$IDAspect[i,"W.OC"])
    ser.e <- as.character(ID$IDAspect[i,"E.OC"])

    sub <- trw.series[,c(ser.n, ser.s, ser.w, ser.e)] # Selects series of those names

    tab.ecc.schw[,c(4*i-3, 4*i-2, 4*i-1, 4*i)]<-.Schwein(ser.n,ser.s,ser.w,ser.e, sub)
    tab.ecc.braam[,c(4*i-3, 4*i-2, 4*i-1, 4*i)]<-.Braam(ser.n,ser.s,ser.w,ser.e, sub)
    tab.ecc.alestalo[,c(4*i-3, 4*i-2, 4*i-1, 4*i)]<-.Alestalo(ser.n,ser.s,ser.w,ser.e, sub)

    colnames.store.schw[,c(4*i-3, 4*i-2, 4*i-1, 4*i)] <- c(ser.n,ser.s,ser.w,ser.e, sub)
    colnames.store.braam[,c(4*i-3, 4*i-2, 4*i-1, 4*i)] <- c(ser.n,ser.s,ser.w,ser.e, sub)
    colnames.store.alestalo[,c(4*i-3, 4*i-2, 4*i-1, 4*i)] <- c(ser.n,ser.s,ser.w,ser.e, sub)

    rm(ser.e, ser.s, ser.w, ser.n) # Clearing memory
  }
  colnames(tab.ecc.schw) <- colnames.store.schw[1,]
  colnames(tab.ecc.braam) <- colnames.store.braam[1,]
  colnames(tab.ecc.alestalo) <- colnames.store.alestalo[1,]

  return(list(Schweingruber=tab.ecc.schw, Braam=tab.ecc.braam, Alestalo=tab.ecc.alestalo))
}

###########################################
### Function for calculation of BAI index,
### A new approach for BAI calculation baseo on approximation of tree ring to ellipse
### Used arguments of this function are stored TRW (data.frame name), plot and tree ID (number)
###########################################

BAIcalculation<-function(trw.series){
  IDs<-.IDdistinct(trw.series)
  result<-data.frame(cambAge=c(1:length(rownames(trw.series))))
  rownames(result)<-rownames(trw.series)
  coln<-"cambAge"
  for(i in 1:length(IDs$IDAspect$IDPlot)){
    plot<-IDs$IDAspect$IDPlot[i]
    tree<-IDs$IDAspect$IDTree[i]
    level<-IDs$IDAspect$IDLevel[i]
    w <- .seriesTRW_one(trw.series,plot,tree,level)
    if(length(w[,1])!=length(w[,2]) | length(w[,1])!=length(w[,3]) | length(w[,1])!=length(w[,4])){w <- .replace(w)}
    widths<-.widthsCalculation(w)
    baiCalculation<-data.frame(date=widths$date,cambAge=widths$cambAge,bai=rep(0,length(widths$date)))
    area<-((pi*abs(widths$W)*abs(widths$N))/4)+((pi*abs(widths$N)*abs(widths$E))/4)+((pi*abs(widths$E)*abs(widths$S))/4)+((pi*abs(widths$S)*abs(widths$W))/4)
    baiCalculation$bai<-area[1]
    for(i in 2:length(baiCalculation$date)){
      j<-i-1
      baiCalculation$bai[i]<-area[i]-area[j]
    }
    coln<-c(coln,paste(plot, tree, level, sep="_"))
    b<-c(rep("<NA>",length(result[,1])-length(baiCalculation[,3])),baiCalculation[,3])
    result<-cbind(result, b)
  }
  colnames(result)<-coln
  result$cambAge <- NULL
  return(result)
}

###########################################
### Function for calculation taper
### Calculates taper for each level
### Used arguments of this function are stored TRW (data.frame name), plot and tree ID (number) and direction of tree core, default is North-South ("N-S"), second Argument is East-West ("E-W")
###########################################

taperCalcul<-function(trw.series,meta){
  
  #trw.series<-trw
  
  IDa<-.IDdistinct(trw.series)
  ID<-IDa$IDAspect
  res<-data.frame(plot=numeric(0),tree=numeric(0),level=numeric(0),level.height=numeric(0),taper=numeric(0),taper.angle=numeric(0))
  
  trees<-paste(ID$IDPlot, ID$IDTree, sep="_")
  not<-unique(trees)
  
  for(i in 1:length(not)){
    split <- strsplit(not[i], "_")
    unl <- unlist(split)
    plot<-as.numeric(unl[1])
    tree<-as.numeric(unl[2])
    sbTree<-subset(ID,IDPlot==plot & IDTree==tree)
    sbTree<-sbTree[order(sbTree$IDLevel),] 
    print(sbTree)
    for(j in 1:(length(sbTree[,1]))){
      print(j)
      ID_1<-subset(sbTree,IDLevel==sbTree$IDLevel[j])
      ID_2<-subset(sbTree,IDLevel==sbTree$IDLevel[(j+1)])
      #print(ID_1)
      #print(ID_2)
      subMeta1<-subset(meta,Plot_ID==ID_1$IDPlot & Tree_ID==ID_1$IDTree & Level_ID==ID_1$IDLevel)
      subMeta2<-subset(meta,Plot_ID==ID_1$IDPlot & Tree_ID==ID_1$IDTree & Level_ID==ID_2$IDLevel)
      if(nrow(subMeta2) == 0){
        subMeta2<-subset(meta,Plot_ID==ID_1$IDPlot & Tree_ID==ID_1$IDTree & Level_ID==999)
      }
      treeHeight1<-subMeta1$Level_cm
      treeHeight2<-subMeta2$Level_cm
      treeHeight<-as.numeric(treeHeight2)-as.numeric(treeHeight1)
      ser.n <- toString(ID_1$N.OC)
      ser.s <- toString(ID_1$S.OC)
      ser.e <- toString(ID_1$E.OC)
      ser.w <- toString(ID_1$W.OC)
      w <- trw.series[,c(ser.n, ser.s, ser.e, ser.w)]
      if(length(w[,1])!=length(w[,2]) | length(w[,1])!=length(w[,3]) | length(w[,1])!=length(w[,4])){w <- .replace(w)}
      widths<-.widthsCalculation(w)
      last<-length(widths$cambAge)
      oNW<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$W[last]*widths$W[last])))/4
      oSW<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$W[last]*widths$W[last])))/4
      oNE<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$E[last]*widths$E[last])))/4
      oSE<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$E[last]*widths$E[last])))/4
      oElipse1<-oNW+oSW+oNE+oSE
      #print(subMeta1)
      #print(subMeta2)
      if(subMeta2$Level_ID!=999){
        ser.n <- toString(ID_2$N.OC)
        ser.s <- toString(ID_2$S.OC)
        ser.e <- toString(ID_2$E.OC)
        ser.w <- toString(ID_2$W.OC)
        w <- trw.series[,c(ser.n, ser.s, ser.e, ser.w)]
        if(length(w[,1])!=length(w[,2]) | length(w[,1])!=length(w[,3]) | length(w[,1])!=length(w[,4])){w <- .replace(w)}
        widths<-.widthsCalculation(w)
        last<-length(widths$cambAge)
        oNW<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$W[last]*widths$W[last])))/4
        oSW<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$W[last]*widths$W[last])))/4
        oNE<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$E[last]*widths$E[last])))/4
        oSE<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$E[last]*widths$E[last])))/4
        oElipse2<-oNW+oSW+oNE+oSE
      }else{
        oElipse2<-0
      }
      rElipse1<-oElipse1/pi
      rElipse2<-oElipse2/pi
      
      r<-data.frame(plot=plot,
                    tree=tree,
                    level=j,
                    level.height=treeHeight1,
                    taper=((rElipse1-rElipse2)/(treeHeight*10)),
                    taper.angle=atan(0.5*((rElipse1-rElipse2)/(treeHeight*10)))*(180/pi))
      
      
      res<-rbind(res,r)
    }
  }
  return(res)
} #Druhá verze téhle funkce

###########################################
### Function for estimating the number of missing rings from pith offset (as metadata tabel)
### There are three options of algorithm for estimating number of missing rings (according to Altman et al. (2016): Forest Ecology and Management, 380:82-89)
### 1] pith offset is divided by the average TRW of last nyrs tree-rings
### 2] pith offset is converted to area (pi*PO^2) and then divided by average BAI of last nyrs tree-rings
### 3] mean of 1] and 2] (adviced by Altman et al. (2016): Forest Ecology and Management, 380:82-89 as the aproach with the lowest bias)
###########################################

EMR <- function(trw.series, p.off, nyrs, method="TRW"){
  j<-1
  tab <- data.frame(MissingRings=NA, Series=NA)

  ### TRW based estimate of number of missing rings
  #########
  if (method=="TRW") {
    for (i in (1:ncol(trw.series))){
      df <- data.frame(na.omit(trw.series[,i])) # trw series with removed NAs
      mean <- mean(df[(1:nyrs),1]) # mean of last nyrs tree-rings

      p.off.2 <- subset(p.off, subset=ID==as.character(data.frame(colnames(trw.series))[i,])) # Subset of pit.offset table to contain only row with tree, which is currently analysed

      if (as.numeric(p.off.2["P.OFFSET"])==0 && !is.na(p.off.2["P.OFFSET"])) {tab[j,1] <- 0; tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,])} # If core goes through pith -> 0 missing rings
      else
      {if (!is.na(p.off.2["P.OFFSET"])) {value <- round(as.numeric(p.off.2["P.OFFSET"])/mean, digits=0); # If core misses the pith -> rounded(offset/meanTRW) missing rings
      tab[j,1] <- value;
      }}

      tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
      j <- j+1
    }
  }

  ### BAI based estimate of number of missing rings
  #########
  if (method=="BAI") {
    bai <- bai.in(trw.series, subset(p.off, select=c(ID,P.OFFSET))) # Calculation of BAI

    # Almost the same script as in case of TRW
    for (i in (1:ncol(bai))){
      df <- data.frame(na.omit(bai[,i])) # bai series with removed NAs
      mean <- mean(df[(1:nyrs),1]) # mean of last nyrs tree-rings

      p.off.2 <- subset(p.off, subset=ID==as.character(data.frame(colnames(bai))[i,])) # Subset of pit.offset table to contain only row with tree, which is currently analysed

      if (as.numeric(p.off.2["P.OFFSET"])==0 && !is.na(p.off.2["P.OFFSET"])) {tab[j,1] <- 0; tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,])}
      else
      {if (!is.na(p.off.2["P.OFFSET"])) {value <- round(pi*(as.numeric(p.off.2["P.OFFSET"])^2)/mean, digits=0); # The area delimited by the pith offset is divided by mean BAI of last nyrs tree-rings
      tab[j,1] <- value;
      }}

      tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
      j <- j+1
    }
  }

  ### Altman et al. (2016) approach - mean of missing rings estimated based on TRW and BAI
  #########
  if (method=="Both") {
    bai <- bai.in(trw.series, subset(p.off, select=c(ID,P.OFFSET))) # Calculation of BAI

    for (i in (1:ncol(bai))){
      df.trw <- data.frame(na.omit(trw.series[,i])) # trw and BAI series with removed NAs
      df.bai <- data.frame(na.omit(bai[,i]))
      mean.trw <- mean(df.trw[(1:nyrs),1]) # mean of last nyrs tree-rings
      mean.bai <- mean(df.bai[(1:nyrs),1])

      p.off.2 <- subset(p.off, subset=ID==as.character(data.frame(colnames(trw))[i,])) # Subset of pit.offset table to contain only row with tree, which is currently analysed

      if (as.numeric(p.off.2["P.OFFSET"])==0 && !is.na(p.off.2["P.OFFSET"])) {tab[j,1] <- 0; tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,])}
      else
      {if (!is.na(p.off.2["P.OFFSET"])) {value <- round(((pi*(as.numeric(p.off.2["P.OFFSET"])^2)/mean.bai)+(as.numeric(p.off.2["P.OFFSET"])/mean.trw))/2, digits=0); # BAI and TRW approaches are calculated simultaneously and averaged
      tab[j,1] <- value;
      }}

      tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
      j <- j+1
    }
  }

  return(tab)

}

###########################################
### Function for replacing missing tree rings
### Several options is necessary to fill
### trw.series - file that contains TRW series
### nyrs - number of years for interpolation
### p.off - pith offset - an estimated distance to the pith
### method - interpolation method
### no.series.per.height - number of cores from the same tree height for interpolation
###########################################

RMR <- function (trw.series, mr.estimate, nyrs=5, nsph=4) {
  
  ID <- .IDdistinct(trw.series)
  trw.series.descend<-trw.series[ order(-as.numeric(row.names(trw.series))), ]
  mr.estimate[is.na(mr.estimate)]<-0
  
  #One core per sampling height
  if (nsph==1){
    trw.result<-data.frame(matrix(NA, ncol = length(trw.series), nrow = (length(trw.series)+max(mr.estimate$MissingRings))))
    names(trw.result)<-names(trw.series)
    for (i in 1:length(trw.series)){
      ser<-trw.series[,i]
      a<-length(ser);nas<-0
      while (is.na(ser[a])==T){
        nas<-nas+1
        a<-a-1
      }
      ser<-na.omit(trw.series[,i])
      avg<-mean(ser[1:nyrs])
      miss.i<-which(mr.estimate$Series==names(trw.series[i]))
      miss<-mr.estimate$MissingRings[miss.i]
      ser2<-c(rep(NA,nas),na.omit(trw.series.descend[,i]),rep(avg,miss))
      trw.result[1:length(ser2),i]<-ser2
    }
  }
  
  #Two cores per sampling height
  if(nsph==2){
    n <- substring(names(trw.series),1,(nchar(names(trw.series))-1))
    mr.estimate$IDCore<-n
    n <- unique(n)
    len<-NULL
    naFront<-NULL
    for(i in 1:length(trw.series)){
      ser<-trw.series[,i]
      a<-length(ser);nas<-0
      while (is.na(ser[a])==T){
        nas<-nas+1
        a<-a-1
      }
      naFront<-c(naFront,nas)
      len<-c(len,(length(na.omit(trw.series[,i]))+nas))
    }
    mr.estimate$Length<-len
    mr.estimate$nasf<-naFront
    print(mr.estimate)
    
    for(i in 1:length(n)){
      ind<-which(n[i]==mr.estimate$IDCore)
      miss<-mr.estimate$MissingRings[ind]
      len<-mr.estimate$Length[ind]
      if(miss[1]==0 & miss[2]==0) m<-len
      if(miss[1]!=0 & miss[2]==0) m<-max(len)
      if(miss[1]==0 & miss[2]!=0) m<-max(len)
      if(miss[1]!=0 & miss[2]!=0) m<-round(mean(len+miss))
      mr.estimate$LengthTotal[ind]<-m
    }
    print(mr.estimate)
    
    trw.result<-data.frame(matrix(NA, ncol = length(trw.series), nrow = max(mr.estimate$LengthTotal)))
    names(trw.result)<-names(trw.series)
    for(i in 1:length(n)){
      ind<-which(n[i]==mr.estimate$IDCore)
      sA<-mr.estimate$Series[ind[1]]
      nasA<-mr.estimate$nasf[ind[1]]
      sB<-mr.estimate$Series[ind[2]]
      nasB<-mr.estimate$nasf[ind[2]]
      serA<-na.omit(trw.series[,sA])
      serB<-na.omit(trw.series[,sB])
      avgA<-mean(serA[1:nyrs])
      avgB<-mean(serB[1:nyrs])
      mA<-mr.estimate$LengthTotal[ind[1]]-mr.estimate$Length[ind[1]]
      mB<-mr.estimate$LengthTotal[ind[2]]-mr.estimate$Length[ind[2]]
      serA2<-c(rep(NA,nasA),na.omit(trw.series.descend[,sA]),rep(avgA,mA))
      serB2<-c(rep(NA,nasB),na.omit(trw.series.descend[,sB]),rep(avgB,mB))
      trw.result[1:length(serA2),sA]<-serA2
      trw.result[1:length(serB2),sB]<-serB2
    }
  }
  
  #Four cores per sampling height
  if(nsph==4){
    n <- substring(names(trw.series),1,(nchar(names(trw.series))-1))
    mr.estimate$IDCore<-n
    n <- unique(n)
    len<-NULL
    naFront<-NULL
    for(i in 1:length(trw.series)){
      ser<-trw.series[,i]
      a<-length(ser);nas<-0
      while (is.na(ser[a])==T){
        nas<-nas+1
        a<-a-1
      }
      naFront<-c(naFront,nas)
      len<-c(len,(length(na.omit(trw.series[,i]))+nas))
    }
    mr.estimate$Length<-len
    mr.estimate$nasf<-naFront
    mr.estimate$LengthTotal<-rep(0,length(mr.estimate$nasf))
    for(i in 1:length(n)){
      ind<-which(n[i]==mr.estimate$IDCore)
      if(length(ind)==4){
        miss<-mr.estimate$MissingRings[ind]
        len<-mr.estimate$Length[ind]
        m<-c(0,0,0,0)
        if(miss[1]==0 & miss[2]==0 & miss[3]==0 & miss[4]==0) m<-len
        if(miss[1]!=0 & miss[2]!=0 & miss[3]!=0 & miss[4]!=0) m[1:4]<-round(median(len+miss))
        if(miss[1]==0 | miss[2]==0 | miss[3]==0 | miss[4]==0) {
          abc<-subset(mr.estimate,IDCore==n[i] & MissingRings==0)
          m[1:4]<-max(abc$Length)
        }
        if(miss[1]==0) m[1]<-len[1]
        if(miss[2]==0) m[2]<-len[2]
        if(miss[3]==0) m[3]<-len[3]
        if(miss[4]==0) m[4]<-len[4]
        mr.estimate$LengthTotal[ind]<-m
      }else{;}
    }
    
    trw.result<-data.frame(matrix(NA, ncol = length(trw.series), nrow = max(mr.estimate$LengthTotal)))
    names(trw.result)<-names(trw.series)
    for(i in 1:length(n)){
      ind<-which(n[i]==mr.estimate$IDCore)
      if(length(ind)==4){
        sA<-mr.estimate$Series[ind[1]]
        nasA<-mr.estimate$nasf[ind[1]]
        sB<-mr.estimate$Series[ind[2]]
        nasB<-mr.estimate$nasf[ind[2]]
        sC<-mr.estimate$Series[ind[3]]
        nasC<-mr.estimate$nasf[ind[3]]
        sD<-mr.estimate$Series[ind[4]]
        nasD<-mr.estimate$nasf[ind[4]]
        serA<-na.omit(trw.series[,sA])
        serB<-na.omit(trw.series[,sB])
        serC<-na.omit(trw.series[,sC])
        serD<-na.omit(trw.series[,sD])
        avgA<-mean(serA[1:nyrs])
        avgB<-mean(serB[1:nyrs])
        avgC<-mean(serC[1:nyrs])
        avgD<-mean(serD[1:nyrs])
        mA<-mr.estimate$LengthTotal[ind[1]]-mr.estimate$Length[ind[1]]
        mB<-mr.estimate$LengthTotal[ind[2]]-mr.estimate$Length[ind[2]]
        mC<-mr.estimate$LengthTotal[ind[3]]-mr.estimate$Length[ind[3]]
        mD<-mr.estimate$LengthTotal[ind[4]]-mr.estimate$Length[ind[4]]
        if(mA<0) mA<-0;if(mB<0) mB<-0;if(mC<0) mC<-0;if(mD<0) mD<-0;
        serA2<-c(rep(NA,nasA),na.omit(trw.series.descend[,sA]),rep(avgA,mA))
        serB2<-c(rep(NA,nasB),na.omit(trw.series.descend[,sB]),rep(avgB,mB))
        serC2<-c(rep(NA,nasC),na.omit(trw.series.descend[,sC]),rep(avgC,mC))
        serD2<-c(rep(NA,nasD),na.omit(trw.series.descend[,sD]),rep(avgD,mD))
        trw.result[1:length(serA2),sA]<-serA2
        trw.result[1:length(serB2),sB]<-serB2
        trw.result[1:length(serC2),sC]<-serC2
        trw.result[1:length(serD2),sD]<-serD2
      }else{
        ind<-which(n[i]==mr.estimate$IDCore)
        for(j in 1:length(ind)){
          k<-ind[j]
          ser<-na.omit(trw.series.descend[,k])
          trw.result[1:length(ser),k]<-ser
        }
      }
    }
  }
  
  start<-as.numeric(rownames(trw.series.descend)[1])
  end<-as.numeric(rownames(trw.series.descend)[1])-length(trw.result[,1])+1
  years<-c(start:end)
  rownames(trw.result)<-years
  trw.result<-trw.result[ order(as.numeric(row.names(trw.result))), ]
  return(trw.result)
}


############################################
### Apical growth chronologies
############################################

apical <- function (trw.series, meta, mr.estimate) {
  j <- 1
  n.ring <- data.frame(ObservedRing=NA, MissingRing=NA, TotalRing=NA, Series=NA)
  n.ring.elev <- data.frame(TotalRing=NA, IDPlot=NA, IDTree=NA, IDLevel=NA, Pith=NA, Height.cm=NA, Speed.cmyr=NA, MeanSpeedError.cmyr=NA)
  n.ring.elev.apex <- data.frame(TotalRing=NA, IDPlot=NA, IDTree=NA, IDLevel=NA, Pith=NA, Height.cm=NA, Speed.cmyr=NA, MeanSpeedError.cmyr=NA)

  ID <- .IDdistinct(trw.series)

  for (i in (1:ncol(trw.series))){
    df.trw <- data.frame(na.omit(trw.series[,i])) # Getting number of non-NA rings from trw series
    miss.ring <- subset(mr.estimate, subset=Series==as.character(data.frame(colnames(trw))[i,])) # Number of missing rings is taken from dataframe previously created using EMR

    # Appending results to output dataframe
    n.ring[j,1] <- nrow(df.trw) # number of measured rings
    if (nrow(miss.ring)==1) {n.ring[j,2] <- miss.ring[1,1]} # number of missing rings
    if (!is.na(miss.ring[1,1])) {n.ring[j,3] <- n.ring[j,2] + n.ring[j,1]} else {n.ring[j,3] <- n.ring[j,1]}  # total number of rings (measured+missing)
    n.ring[j,4] <- as.character(data.frame(colnames(trw.series))[i,]) # name of series
    j <- j+1
  }

  for (k in (1:nrow(ID$IDAspect))) {
    ser.n <- as.character((ID$IDAspect)[k,"N.OC"])
    ser.s <- as.character((ID$IDAspect)[k,"S.OC"])
    ser.e <- as.character((ID$IDAspect)[k,"E.OC"])
    ser.w <- as.character((ID$IDAspect)[k,"W.OC"])

    subs <- subset(n.ring, subset=(Series==ser.n|Series==ser.s|Series==ser.e|Series==ser.w)) # subset of estimated number of tree-rings for cores from the same level
    pith <- subset(subs, subset=MissingRing==0) # subset of number of tree-rings from cores from the same level, which hit the pit
    ### If at least one core goes through the pit, then it is used. Othervise, we use median of estimated number of rings based on all available cores.
    if (nrow(pith)>0) {ring.estimate <- median(pith[,"TotalRing"]);
    n.ring.elev[k,"Pith"] <- "Yes"}
    if (nrow(pith)==0) {ring.estimate <- median(subs[,"TotalRing"]);
    n.ring.elev[k,"Pith"] <- "No"}

    n.ring.elev[k,"TotalRing"] <- round(ring.estimate, 0)
    n.ring.elev[k,c("IDPlot", "IDTree", "IDLevel")] <- ((ID$IDAspect)[k,c("IDPlot", "IDTree", "IDLevel")])

    meta.subs <- subset(meta, subset=(meta[,1]==n.ring.elev[k,"IDPlot"] & meta[,2]==n.ring.elev[k,"IDTree"] & meta[,3]==n.ring.elev[k,"IDLevel"]))

    if (nrow(meta.subs)==1) {n.ring.elev[k,"Height.cm"] <- meta.subs[1,4]}

  }

  # Subsetting tree heights
  apex <- subset(meta, subset=meta$Level_ID==999)
  for (j in 1:nrow(apex)) {
    n.ring.elev.apex[j,c("IDPlot", "IDTree", "IDLevel", "Height.cm")] <- apex[j,c("Plot_ID", "Tree_ID", "Level_ID", "Level_cm")]
    n.ring.elev.apex[j,"TotalRing"] <- 0 # tree top has no rings
    n.ring.elev.apex[j,"Pith"] <- "Apex"


  }

  n.ring.elev <- rbind(n.ring.elev, n.ring.elev.apex)
  n.ring.elev <- n.ring.elev[order(n.ring.elev$IDPlot, n.ring.elev$IDTree, n.ring.elev$IDLevel),] # Sorting of "Level" table

  for (l in (2:nrow(n.ring.elev))) {
    if ( n.ring.elev[l,"IDPlot"]==n.ring.elev[(l-1),"IDPlot"] & n.ring.elev[l,"IDTree"]==n.ring.elev[(l-1),"IDTree"] ) {

      n.ring.elev[l,"Speed.cmyr"] <- (n.ring.elev[l,"Height.cm"]-n.ring.elev[(l-1),"Height.cm"]) / (n.ring.elev[(l-1),"TotalRing"]-n.ring.elev[l,"TotalRing"])
      n.ring.elev[l,"MeanSpeedError.cmyr"] <- n.ring.elev[l,"Speed.cmyr"] - (n.ring.elev[l,"Height.cm"]-n.ring.elev[(l-1),"Height.cm"]) / (n.ring.elev[(l-1),"TotalRing"]-n.ring.elev[l,"TotalRing"] + 1)

    }


  }

  lst <- list(N.ring_Core=n.ring, N.ring_Level=n.ring.elev)
  return(lst)
}
