###########################################
### Set of function designed for easier data managing
###
### Function .IDdistinct         - separates codes of PlotID, TreeID, Elevation level of coring and aspect (direction of coring) from series name
###          .IDdistinct_simple  - ligth version of function .IDdistinct
###          .IDdistinct_medium  - formating of ID code in meta file
###          .extractTreeHeights - extract tree heights from metadata file
###          .extractTRW         - extract tree ring data for specific tree
###
### Important for future linking of series with metadata + calculation of excentricity indexes for different elevation levels
### The format of series code should be IDPlots_IDTrees_IDElevationAspect
###########################################

.IDdistinct <- function(trw.series, complete=FALSE)
{
  tab <- data.frame(colnames(trw.series), NA, NA, NA, NA) # Creates table with original series names as rows
  colnames(tab) <- c("Original.Code", "IDPlot", "IDTree", "IDLevel", "Aspect") # Appends 4 new columns (IDPlot, ...) to the table
  series.names <- (data.frame(colnames(trw.series)))

  for (i in 1:ncol(trw.series)) # Looping one-by-one through the list of series
  {
    an.series <- as.character(series.names[i,])
    Aspect <- (substr(an.series, nchar(an.series),nchar(an.series))) # Extracts the last character of original series name (i.e., aspect)...
    tab[i,5] <- Aspect # ... and apends it to the output table.
    split <- strsplit(substr(an.series,1,nchar(an.series)-1), "_") # Splits the rest of the string using _ separator and ...
    unlist <- t(data.frame(as.numeric(unlist(split)))) # ... creates a data frame from it.
    tab[i,2] <- unlist[1,1] # Different parts of data frame are appendend to output table.
    tab[i,3] <- unlist[1,2]
    tab[i,4] <- unlist[1,3]
    rm(an.series, Aspect, split, unlist) # Memory clearing.
  }

  N <- subset(tab, subset=Aspect=="N"); colnames(N) <- c("N.OC", "IDPlot", "IDTree", "IDLevel", "N") # Subsetting of table to 4 tables for different aspects
  S <- subset(tab, subset=Aspect=="S"); colnames(S) <- c("S.OC", "S.IDPlot", "S.IDTree", "S.IDLevel", "S")
  E <- subset(tab, subset=Aspect=="E"); colnames(E) <- c("E.OC", "E.IDPlot", "E.IDTree", "E.IDLevel", "E")
  W <- subset(tab, subset=Aspect=="W"); colnames(W) <- c("W.OC", "W.IDPlot", "W.IDTree", "W.IDLevel", "W")

  if (complete==TRUE)  {
    merge <- merge(N[,1:4],S[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("S.IDPlot", "S.IDTree", "S.IDLevel"), all=T) # Merging different aspects into one file
    merge <- merge(merge,W[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("W.IDPlot", "W.IDTree", "W.IDLevel"), all=T)
    merge <- merge(merge,E[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("E.IDPlot", "E.IDTree", "E.IDLevel"), all=T) }
  else {
    merge <- merge(N[,1:4],S[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("S.IDPlot", "S.IDTree", "S.IDLevel"), all=F) # Merging different aspects into one file
    merge <- merge(merge,W[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("W.IDPlot", "W.IDTree", "W.IDLevel"), all=F)
    merge <- merge(merge,E[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("E.IDPlot", "E.IDTree", "E.IDLevel"), all=F) }

  lst <- list(IDCore=tab, IDAspect=merge)
  return(lst)
}

.IDdistinct_medium<-function(file){
  tab <- data.frame(file["ID"], NA, NA, NA)
  colnames(tab) <- c("Code", "Plot_ID", "Tree_ID", "Level_ID") 
  series.names <- tab$Code
  for (i in 1:length(series.names)) # Looping one-by-one through the list of series
  {
    an.series <- as.character(series.names[i])
    split <- strsplit(an.series, "_") # Splits the rest of the string using _ separator and ...
    unlist <- t(data.frame(as.numeric(unlist(split)))) # ... creates a data frame from it.
    tab[i,2] <- unlist[1,1] # Different parts of data frame are append to output table.
    tab[i,3] <- unlist[1,2]
    tab[i,4] <- unlist[1,3]
    rm(an.series, split, unlist) # Memory clearing.
  }
  mer <- merge(file, tab, by.x="ID", by.y="Code")
  mer <- mer[,-c(1)]
  mer <- mer[c("Plot_ID", "Tree_ID", "Level_ID", "Level_cm", "Type.of.the.sample")]
  return(mer)
}

.IDdistinct_simple<-function(file){
  tab <- data.frame(colnames(file), NA, NA) # Creates table with original series names as rows
  tab<-tab[-c(1),]
  colnames(tab) <- c("Code", "IDPlot", "IDTree") # Appends 4 new columns (IDPlot, ...) to the table
  series.names <- tab$Code
  for (i in 1:length(series.names)) # Looping one-by-one through the list of series
  {
    an.series <- as.character(series.names[i])
    split <- strsplit(substr(an.series,1,nchar(an.series)-1), "_") # Splits the rest of the string using _ separator and ...
    unlist <- t(data.frame(as.numeric(unlist(split)))) # ... creates a data frame from it.
    tab[i,2] <- unlist[1,1] # Different parts of data frame are append to output table.
    tab[i,3] <- unlist[1,2]
    #tab[i,4] <- unlist[1,3]
    rm(an.series, split, unlist) # Memory clearing.
  }
  return(tab)
}#ID distinct for plot and tree only (no aspects etc.)

.extractTreeHeights<-function(meta,plot=1,tree=1){
  
  meta <- .IDdistinct_medium(meta)
  
  tree<-subset(subset(meta,IDPlot==plot & IDTree==tree))
  return(tree$Height)
} #extract tree heights from meta data

.extractTRW<-function(trw,plot,tree,direction=""){
  tree<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  ser<-.seriesTRW_one(trw)
  for(i in 2:length(tree$IDLevel)){
    a<-.seriesTRW_one(trw,level=i)
    ser<-cbind(ser,a)
  }
  if(direction==""){
    serie<- ser
  }
  if(direction=="North"){
    serie<- subset(ser, select=tree$N.OC)
  }
  if(direction=="South"){
    serie<- subset(ser, select=tree$S.OC)
  }
  if(direction=="East"){
    serie<- subset(ser, select=tree$E.OC)
  }
  if(direction=="West"){
    serie<- subset(ser, select=tree$W.OC)
  }
  return(serie)
}# extracts TRW for one chosen tree

######################
### Group of not user-available functions used in excentricity calculation
###	CreateTable to save (i) excentricity indexes and (ii) storing names of input series
###	CalCUlation of excentricity indexes according to Schweingruber (1996), Braam et al. (1987) and Alestalo et al. (1971)
######################

.CreateTableTRW<-function(ser,ID){
  table<- data.frame(matrix(nrow=nrow(ser), ncol=(4*nrow(ID)))) # Dataframe to store results ...
  rownames(table) <- rownames(ser) # ... number of rows equal to maximum number of tree-rings in input series
  return(table)
}
.CreateTableMeta<-function(ID){
  table <- data.frame(matrix(nrow=1, ncol=(4*nrow(ID)))) # ... number of columns is 4*number of elevation levels
  return(table)
}

######################
### Group of functions for calculation of excentricity idexes
######################
.Schwein<-function(ser.n,ser.s,ser.w,ser.e, sub){
  Schwein_ind<- c(sub[,1]/sub[,2], sub[,2]/sub[,1], sub[,3]/sub[,4], sub[,4]/sub[,3]) # Calculates excentrecity index and assign it to the table
  return(Schwein_ind)
}
.Braam<-function(ser.n,ser.s,ser.w,ser.e, sub){
  # Because we have got 2 perpendicular measurements, I use average of them !
  Braam_ind <- c((sub[,1]-(sub[,3]+sub[,4])/2)/(sub[,1]+(sub[,3]+sub[,4])/2),
                 (sub[,2]-(sub[,3]+sub[,4])/2)/(sub[,2]+(sub[,3]+sub[,4])/2),
                 (sub[,3]-(sub[,1]+sub[,2])/2)/(sub[,3]+(sub[,1]+sub[,2])/2),
                 (sub[,4]-(sub[,1]+sub[,2])/2)/(sub[,4]+(sub[,1]+sub[,2])/2))
  return(Braam_ind)
}
.Alestalo<-function(ser.n,ser.s,ser.w,ser.e, sub){
  # Because we have got 2 perpendicular measurements, I use average of them !
  Alest_ind<- c(sub[,1]/(sub[,1]+sub[,2]),
                sub[,2]/(sub[,1]+sub[,2]),
                sub[,3]/(sub[,3]+sub[,4]),
                sub[,4]/(sub[,3]+sub[,4]))
  return(Alest_ind)
}

######################
### Functions which select names of series or data (excentricity or TRW) for one particular tree
###	Tree selection is made in arguments
### Other functions in this group extracts data from individual data frames or lists
######################

.treeSelect<-function(plot,tree,sourc){
  subset(.IDdistinct(method)$IDAspect,IDPlot==plot & IDTree==tree)
} # Select names of series for particullar tree
.seriesTRW_one<-function(t,plot=1,tree=1,level=1){
  ID<-.IDdistinct(t)
  ID<-subset(ID$IDAspect,IDPlot==plot & IDTree==tree & IDLevel==level)
  ser.n <- toString(ID$N.OC[1])
  ser.s <- toString(ID$S.OC[1])
  ser.e <- toString(ID$E.OC[1])
  ser.w <- toString(ID$W.OC[1])
  w <- t[,c(ser.n, ser.s, ser.e, ser.w)]
  return(w)
}# returns TRW for selected tree including NA values
.seriesTRW_two<-function(t,tree=1,plot=1,level=1){
  w <- .seriesTRW_one(t,plot,tree,level)
  l1<-na.omit(w[,1]);l2<-na.omit(w[,2]);l3<-na.omit(w[,3]);l4<-na.omit(w[,4])
  cut<-c((length(w[,1])-length(l1)),(length(w[,1])-length(l2)),(length(w[,1])-length(l3)),(length(w[,1])-length(l4)))
  cut<-min(cut)
  if(cut>0){cut<-cut+1}else{cut<-0}
  w<-w[cut:length(w[,1]),]
  return(w)
}# returns TRW for selected tree without NA values

.seriesLength<-function(serie){
  lengths<-c()
  for(i in 1:(length(serie))){
    lengths[i]<-length(na.omit(serie[,i]))
  }
  i<-length(lengths)+1
  lengths[i]<-0
  return(lengths)
}# calculates and return series length
.calcSumTRW<-function(serie,lengths){
  trws<-c()
  apex<-data.frame(apx=rep("NA",length(serie[,1])))
  serie2<-cbind(serie,apex)
  serie3<-cbind(serie,apex)
  for(i in 1:(length(serie2)-1)){
    serie3<-cbind(serie3,serie2)
  }

  for(i in 1:length(lengths)){
    x<-length(serie3[,i])
    y<-x-lengths[i]
    if(x==y){
      e<-0
    }else{
      e<-sum(serie3[y:x,i], na.rm = T)
    }
    trws<-c(trws,e)
  }
  return(trws)
}# calculates series sum for selected series and length

.widthCalculation<-function(s){
  s[is.na(s)] <- 0
  vypocet<-rep(0,length(s))
  vypocet[1]<-s[1]
  for(i in 2:length(s)){
    vypocet[i]<-vypocet[i-1]+s[i]
  }
  return(vypocet)
} #calculates tree-ring distance from the pith for one serie
.widthsCalculation<-function(width){
  widthsN<-.widthCalculation(width[,1])
  widthsS<-.widthCalculation(width[,2])
  widthsE<-.widthCalculation(width[,3])
  widthsW<-.widthCalculation(width[,4])
  w<-data.frame(date=rownames(width), cambAge=c(1:length(widthsN)),N=widthsN,S=widthsS,E=widthsE,W=widthsW)
  return(w)
}# returns data.frame of distance of each tree ring from the pith - serves mainly for BAI calculation and cross section vizualization
.replace<-function(repl){
  NAN<-0;NAS<-0;NAE<-0;NAW<-0
  ISNAN<-F;ISNAS<-F;ISNAE<-F;ISNAW<-F
  lntghs<-c(na.omit(length(repl[,1])),na.omit(length(repl[,2])),na.omit(length(repl[,3])),na.omit(length(repl[,4])))
  for(i in 1:max(lntghs)){
    if(ISNAN==F){if(is.na(repl[i,1])){NAN<-NAN+1}else{ISNAN<-T}}
    if(ISNAS==F){if(is.na(repl[i,2])){NAS<-NAS+1}else{ISNAS<-T}}
    if(ISNAE==F){if(is.na(repl[i,3])){NAE<-NAE+1}else{ISNAE<-T}}
    if(ISNAW==F){if(is.na(repl[i,4])){NAW<-NAW+1}else{ISNAW<-T}}
  }
  r<-repl
  if(NAN>0){a<-NAN+1;b<-NAN+6;r[1:NAN,1]<-mean(repl[a:b,1])}
  if(NAS>0){a<-NAS+1;b<-NAS+6;r[1:NAS,2]<-mean(repl[a:b,2])}
  if(NAE>0){a<-NAE+1;b<-NAE+6;r[1:NAE,3]<-mean(repl[a:b,3])}
  if(NAW>0){a<-NAW+1;b<-NAW+6;r[1:NAW,4]<-mean(repl[a:b,4])}
  return(r)
}# function for filling missing data - complete tree ring series (with an average of five last known tree rings) - serves mainly for BAI calculation and cross section vizualization

######################
### Private functionn for calculation of points for cross-section profile graph
###	Cross section is approximate as an conjunction of four quaters of four diferent elipses
### for each quater of Cartesian coordinate system is used different elipse given calculated from tvo specific cores
######################
.elipseFun <- function(sirky=widths[1,], npoints=100){
  center<-c(0,0)
  t <- seq(0*pi, 2*pi, length.out=npoints)
  q1 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))
  q2 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))
  q3 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))
  q4 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))

  limA<-1;limB<-npoints*0.25
  q1$x <- center[1] + abs(sirky$W)*cos(t[limA:limB])
  q1$y <- center[1] + abs(sirky$N)*sin(t[limA:limB])
  limA<-(npoints*0.25)+1;limB<-npoints*0.5
  q2$x <- center[1] + abs(sirky$E)*cos(t[limA:limB])
  q2$y <- center[1] + abs(sirky$N)*sin(t[limA:limB])
  limA<-(npoints*0.75)+1;limB<-npoints
  q3$x <- center[1] + abs(sirky$W)*cos(t[limA:limB])
  q3$y <- center[1] + abs(sirky$S)*sin(t[limA:limB])
  limA<-(npoints*0.5)+1;limB<-npoints*0.75
  q4$x <- center[1] + abs(sirky$E)*cos(t[limA:limB])
  q4$y <- center[1] + abs(sirky$S)*sin(t[limA:limB])

  vystup<-rbind(q1,q2,q4,q3)
  return(vystup)
}

######################
### Set of private functions for creating graph of stem allometry
###	Twqo main functions are .dataGraphAlometry draw graph of stem alometry
###                         .seriesCalc calculates supporting data
### Result of these functions is one graph of stem alometry
### Funcrions are called twice by function drawEccentricityGraph to create complex graph of stem alometry from all sampled directions
######################
.dataGraphAlometry<-function(trw,vysky,plot=1,tree=1,dir="N-S"){
  if(dir=="N-S"){
    a<-1
    b<-2
  }else{
    a<-3
    b<-4
  }
  tr<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  ser<-.seriesTRW_one(trw,tree=tree,plot=plot,level=1)
  left<-data.frame(ser[a])
  right<-data.frame(ser[b])
  for(i in 2:length(tr$IDLevel)){
    lvl<-i
    ser<-.seriesTRW_one(trw,tree=tree,plot=plot,level=lvl)
    left<-cbind(left,ser[a])
    right<-cbind(right,ser[b])
  }
  l<-.seriesCalc(left,treeHeights);r<-.seriesCalc(right,treeHeights)
  l$trwSum<-l$trwSum*-1;l$length<-l$length*-1
  e<-l$height*-1;e<-c(e,r$heigh)
  LR<-rbind(l,r);LR<-cbind(LR,e)
  return(LR)
}
#set of drawing functions
.linie<-function(graf, vyska, dataLeft, dataRight){
  line1<-.calcSerieLeft(dataLeft,vyska)
  line2<-.calcSerieRight(dataRight,vyska)
  graf<-graf+geom_line(data=line1, aes_string(line1$roky, line1$hodnoty), size=1, colour="#CC0000")
  graf<-graf+geom_line(data=line2, aes_string(line2$roky, line2$hodnoty), size=1, colour="#CC0000")
  return(graf)
}
.plotGraph<-function(plot,tree,heights,trw,ecc,withEccentricity=F,direction="North-South",rangeX){
  ID<-.IDdistinct(trw, F)
  subID<-subset(ID$IDAspect,IDPlot==plot & IDTree==tree)
  subID<-subID[order(subID$IDLevel),]
  subTRW<-trw[,subID$N.OC]
  datetTo<-max(as.numeric(rownames(subTRW)))
  
  TotalRing<-NULL
  for(i in 1:length(subID$IDPlot)){
    tr<-c(as.character(subID$N.OC[i]),
          as.character(subID$S.OC[i]),
          as.character(subID$W.OC[i]),
          as.character(subID$E.OC[i]))
    t<-trw[,tr]
    lngths<-c(length(na.omit(t[,1])),
              length(na.omit(t[,2])),
              length(na.omit(t[,3])),
              length(na.omit(t[,4])))
    l<-round(median(lngths)+0.01)
    TotalRing<-c(TotalRing,l)
  }
  calendarYears<-rep(datetTo,length(TotalRing))
  calendarYears<-calendarYears-TotalRing
  calendarYears<-c(calendarYears,datetTo)
  left<-0
  right<-0
  met<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  left<-met$N.OC;right<-met$S.OC
  if(direction!="North-South"){left<-met$W.OC;right<-met$E.OC}
  serieLeft<-.dataLine(ecc,left)
  serieRight<-.dataLine(ecc,right)
  val<-.widths(trw,left,right,heights)
  graph<-ggplot(val, aes(width, height)) + geom_point(size = 4, color="#336600")+ geom_path(size=1.5,color="#336600")
  graph<-graph+geom_vline(xintercept =0, linetype='dashed', size=1.5, color = "#663300")
  graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  if(withEccentricity==T){
    for(i in 1:ncol(serieLeft)){
      graph<-.linie(graph,heights[i],serieLeft[,i],serieRight[,i])
    }
  }else{
    ;
  }
  graph<-graph+scale_x_continuous(name="Stem width (cm)",limits=rangeX, breaks=seq(rangeX[1],rangeX[2],((abs(rangeX[1])+abs(rangeX[2]))/10)), labels=abs(seq(rangeX[1],rangeX[2],((abs(rangeX[1])+abs(rangeX[2]))/10))))
  graph<-graph+scale_y_continuous(name="Stem height (cm)",limits=c(0,ceiling(max(heights)/100)*100),breaks=heights,labels=heights,
                                  sec.axis=sec_axis(~.,breaks=heights,labels=calendarYears,name="Calendar year"))
  return(graph)
}
.dataLine<-function(ecc,met){
  values<-data.frame(exc[as.character(met[1])])
  for(i in 2:length(met)){
    values<-cbind(values,ecc[as.character(met[i])])
  }
  colnames(values)<-met
  return(values)
} #Gets data for tree excentricity


.seriesCalc<-function(serie,heights){
  lengths<-.seriesLength(serie)
  length2<-c()
  for(i in 1:(length(serie)+1)){
    l<-c(lengths[i:length(lengths)],rep(0,times=i-1))
    length2<-c(length2,l)
  }
  year<-c()
  lab<-as.numeric(rownames(serie))
  yr<-length2[0:(length(serie))+1]
  for(i in length(yr):1){
    a<-length(serie[,1])-yr[i]
    z<-rep(lab[a],times=length(yr))
    year<-c(year,z)
  }

  trws<-.calcSumTRW(serie,length2)
  curve<-data.frame(curve=rep(0:(length(serie)), rep(length(serie)+1,length(serie)+1)),
                    year=year,
                    level=rep(1:(length(serie)+1),length(serie)+1),
                    height=rep(heights, times = length(serie)+1),
                    length=length2,
                    trwSum=trws)
  return(curve)
} # return data for drawing eccentricity index within alometry graph

######################
### Set of other calculation functions
### Supportive functions of various different functions
######################
.widthGraphEccentricity<-function(trw,meta,vysky){
  
  meta <- .IDdistinct_medium(meta)
  
  sirky<-rep(0,length(vysky))
  for(i in 1:length(meta)){
    sirky[i]<-sum(trw[,as.character(meta[i])], na.rm = T)
  }
  sirky[length(vysky)]<-0
  return(data.frame(height=vysky,width=sirky))
}
.dataLine<-function(exct,meta){
  values<-data.frame(exct[,as.character(meta[1])])
  for(i in 2:length(meta)){
    values<-cbind(values,exct[,as.character(meta[i])])
  }
  colnames(values)<-meta
  return(values)
}
.widths<-function(trw,left,right,heights){
  widthsLeft<-.widthGraphEccentricity(trw,left,heights)
  widthsRight<-.widthGraphEccentricity(trw,right,heights)
  widthsLeft$width<-widthsLeft$width*-2
  widthsRight$width<-widthsRight$width*2
  widthsRight<-widthsRight[order(-widthsRight$height),]
  val<-rbind(widthsLeft,widthsRight)
  return(val)
}
.calcSerieRight<-function(serie,vyska){
  serie<-na.omit(serie)
  linie<-data.frame(roky=c(1:length(serie)),hodnoty=serie*10)
  linie$roky<-(linie$roky-1)
  linie$hodnoty<-linie$hodnoty+(vyska-(max(linie$hodnoty)-min(linie$hodnoty))/2)
  return(linie)
}
.calcSerieLeft<-function(serie,vyska){
  serie<-na.omit(serie)
  linie<-data.frame(roky=c(1:length(serie)),hodnoty=serie*10)
  linie$roky<-(linie$roky-1)*-1
  linie$hodnoty<-linie$hodnoty+(vyska-(max(linie$hodnoty)-min(linie$hodnoty))/2)
  return(linie)
}



