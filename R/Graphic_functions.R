######################
### Function drawing graph with stem excentricity in each level
###	To ensure the function, library ggplot2 needs to be installed
###	Graph can be displayed with or without eccentricity curves
### Arguments of this function are: trw (variable containing trii-ring series),
###                                 ecc (variable containing calculated eccentricity - result of function Eccentricity),
###                                 arguments plot and tree (specify which tree data will be vizualized),
###                                  meta (variable containing metadata including height of sampling points),
###                                 withEccentricity (specify if stem eccentricity curves will be vizualized) ,
###                                 and argument method (specify according which method will be vizualized eccentricity curves).
######################

drawEccentricityGraph<-function(trw,ecc,meta,plot,tree,withEccentricity=F,method="Schweingruber"){
  
  meta <- .IDdistinct_medium(meta)

  Schw<-ecc$Schweingruber
  Braa<-ecc$Braam
  Ales<-ecc$Alestalo

  height<-subset(meta,Plot_ID==plot & Tree_ID==tree)
  height<-height$Level_cm

  m<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  north<-m$N.OC;south<-m$S.OC
  val1<-.widths(trw,north,south,height)
  east<-m$E.OC;west<-m$W.OC
  val2<-.widths(trw,east,west,height)
  range<-c(min(val1$width),min(val2$width),max(val1$width),max(val2$width))
  rangeX<-c((ceiling(max(range)/10)*-10),(ceiling(max(range)/10)*10))

  if(method=="Schweingruber"){
    g1<-.plotGraph(plot,tree,height,trw,Schw,withEccentricity,direction="North-South",rangeX)
    g2<-.plotGraph(plot,tree,height,trw,Schw,withEccentricity,direction="East-West",rangeX)
    lab<-c("Schweingruber - N-S","Schweingruber W-E")
  }else{;}
  if(method=="Braam"){
    g1<-.plotGraph(plot,tree,height,trw,Braa,withEccentricity,direction="North-South",rangeX)
    g2<-.plotGraph(plot,tree,height,trw,Braa,withEccentricity,direction="East-West",rangeX)
    lab<-c("Braam - N-S","Braam W-E")
  }else{;}
  if(method=="Alestalo"){
    g1<-.plotGraph(plot,tree,height,trw,Ales,withEccentricity,direction="North-South",rangeX)
    g2<-.plotGraph(plot,tree,height,trw,Ales,withEccentricity,direction="East-West",rangeX)
    lab<-c("Alestalo - N-S","Alestalo W-E")
  }else{;}
  plot_grid(g1, g2, labels=lab, ncol = 2, nrow = 1)
}

######################
### Function drawing graph of cross-section profile
###	Arguments of this function are: trw (variable containing tree-ring series),
###                                 arguments plot and tree (specify which tree data will be vizualized),
###                                 evel (specify which sampling point of given tree will be vizualized),
###                                 show.legent (specify if the legend will be vizualized).
######################

drawCrossSectionProfile<-function(trw,plot=1,tree=1,level=1,show.legend=T){
  w <- .seriesTRW_two(trw,tree,plot,level)
  widths<-.widthsCalculation(w)
  mx<-max(widths[,2:5])
  lim<-c((ceiling(mx/10)*-10),(ceiling(mx/10)*10))
  widths$S<-widths$S*-1;widths$E<-widths$E*-1
  elipse<-.elipseFun(widths[1,])
  for(i in 2:length(widths$cambAge)){
    elipse2<-.elipseFun(widths[i,])
    elipse<-rbind(elipse,elipse2)
  }
  graph<-ggplot(elipse,aes(x,y)) + geom_point(aes(colour = factor(calYear)),size=0.01)+geom_path(aes(colour = factor(calYear)),size=0.5)+scale_fill_brewer(palette="Spectral")
  graph<-graph + scale_x_continuous(limits = lim, breaks=seq(lim[1],lim[2],10), labels=abs(seq(lim[1],lim[2],10))) + scale_y_continuous(limits = lim, breaks=seq(lim[1],lim[2],10), labels=abs(seq(lim[1],lim[2],10)))
  graph<-graph + labs(colour="Calendar year", x="East     <-  ->     West     (mm)", y="South     <-  ->     North     (mm)")
  if(show.legend==T){graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))}
  else{graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")}
  return(graph)
}

######################
### Function drawing graph of BAI in different levels
###	To ensure the function, library ggplot2 needs to be installed
###	Arguments of this function are: baiFile (variable containing calculated data of basal area increment - result of fulction BAIcalculation),
###                                 arguments plot and tree (specify which tree data will be vizualized),
###                                 logscale (specify if logarithmical scalwill be used for vizualization data on y-axis),
###                                 show.legent (specify if the legend will be vizualized).
######################

drawBai<-function(baiFile,plot=1,tree=1,logscale=F,show.legend=T){
  treeForPlot<-subset(.IDdistinct_simple(baiFile),IDPlot==plot & IDTree==tree)
  baiForPlot<-subset(baiFile,select=treeForPlot$Code)

  forPlot <- data.frame(year=integer(),no=double(),x=double())
  for(i in 1:length(baiForPlot[1,])){
    a<-na.omit(baiForPlot[i])
    a[a == 0] <- NA
    a<-na.omit(a)
    df1<-data.frame(year=rownames(a),no=rep(i,length(a[,1])),bai=a[,1])
    forPlot<-rbind(forPlot,df1)
  }

  if(logscale==F){graph<-ggplot(data=forPlot,aes(x=as.numeric(as.character(year)), y=as.numeric(as.character(bai)),group=no,colour=no)) + geom_line(size=1.25)}
  else{graph<-ggplot(data=forPlot,aes(x=as.numeric(as.character(year)), y=as.numeric(as.character(bai)),group=no,colour=no)) + geom_line(size=1.25)+ coord_trans(y = "log10")}
  graph<-graph + labs(colour="Sampling level", x="Calendar year", y="BAI (sq. cm)")
  if(show.legend==T){
    graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }else{
    graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
  }
  y<-as.numeric(levels(forPlot$year))
  lim<-c(min(y),max(y))
  m<-(round((lim[2]-lim[1])/10)+1)*10
  lim[1]<-lim[2]-m
  graph<-graph+scale_x_continuous(limits=lim,breaks=seq(lim[1],lim[2],m/10))
  return(graph)
}

######################
### Function drawing graph of taper for entire tree
###	To ensure the function, library ggplot2 needs to be installed
###	Arguments of this function are: taperFile (variable containing calculated data of basal area increment - result of fulction taperCalcul),
###                                 arguments plot and tree (specify which tree data will be vizualized),
###                                 variant (specify if will be vizualized taper or taper angle - both variant are analogous).
######################

drawTaper<-function(taperFile,plot=1,tree=1,variant="Taper"){
  whatDraw<-taperFile[taperFile$plot==plot,]
  whatDraw<-whatDraw[whatDraw$tree==tree,]
  if(variant=="Taper"){
    p<-ggplot(data=whatDraw, aes(x=level.height, y=taper, group=1)) + geom_line(size=1.3) + geom_point(size=5)
    p<-p+ labs(colour="Stem taper", x="Sampling height", y="Taper [%]")
  }
  if(variant=="Angle"){
    p<-ggplot(data=whatDraw, aes(x=level.height, y=taper.angle, group=1)) + geom_line(size=1.3) + geom_point(size=5)
    p<-p+ labs(colour="Stem taper", x="Sampling height", y="Taper angle [°]")
  }
  range<-c(0,(ceiling(max(whatDraw$level.height)/10)*10))
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p<-p + scale_x_continuous(limits = range, breaks=whatDraw$level.height, labels=whatDraw$level.height)
  return(p)
}

######################
### Function drawing graph of apical growth
###	Arguments of this function are: taperFile (apicalData containing calculated data of apical growth - result of fulction apical),
###                                 arguments plot and tree (specify which tree data will be vizualized).
######################

drawApicalData<-function(trw.series,apicalData,plot=1,tree=1){
  
  ID<-.IDdistinct(trw.series, F)
  subID<-subset(ID$IDAspect,IDPlot==plot & IDTree==tree)
  subTRW<-trw.series[,subID$N.OC]
  datetTo<-max(as.numeric(rownames(subTRW)))-1
  
  dta<-subset(apicalData$N.ring_Level, IDPlot==plot & IDTree==tree)
  dta$Speed.cmyr[dta$Speed.cmyr < 0] <- 0
  dots<-data.frame(where=rep(0,(dta$TotalRing[1]-dta$TotalRing[length(dta$TotalRing)]+1)),
                   howmany=seq(dta$TotalRing[1],dta$TotalRing[length(dta$TotalRing)],-1))
  calendarYears<-rep(datetTo,length(dta$IDLevel))
  calendarYears<-calendarYears-dta$TotalRing+1
  calendarYears<-c("Base",calendarYears)
  
  d<-NULL
  for(i in 2:length(dta$TotalRing)){
    j<-i-1
    
    d1<-seq(from=dta$Height.cm[j],to=(dta$Height.cm[i]-dta$Speed.cmyr[i]),by=dta$Speed.cmyr[i])
    d<-c(d,d1)
  }
  
  d<-c(d,dta$Height.cm[length(dta$TotalRing)])
  d<-round(d,1)
  dots$where<-d
  
  koef<-floor(max(dots$howmany)*0.1+1)*10
  maxError<-floor(((max(dta$Speed.cmyr,na.rm=T)+max(dta$MeanSpeedError.cmyr,na.rm=T))/10)+1)*10
  maxX<-floor(max(dta$Height.cm)*0.01+1)*100
  
  
  g1<-ggplot(data=dots,aes(x=howmany,y=where))+geom_area()
  g1<-g1 + scale_x_continuous(limits = c(0,koef), breaks=seq(0,koef,(koef/5)), labels=seq(0,koef,(koef/5)),"Number of tree rings")
  g1<-g1 + scale_y_continuous(limits = c(0,maxX), breaks=c(0,dta$Height.cm), labels=c(0,dta$Height.cm),"Stem height (cm)",
                              sec.axis=sec_axis(~.,breaks=c(0,dta$Height.cm),labels=calendarYears,name="Calendar year"))
  g1
  
  grp<-dta$Speed.cmyr
  grp<-replace(grp, grp<0, "B")
  grp<-replace(grp, grp>=0, "A")
  dta_ap<-data.frame(height=dta$Height.cm,
                     speed=dta$Speed.cmyr,
                     error=dta$MeanSpeedError.cmyr,
                     group=grp)
  group.colors <- c(A = "#000000", B = "#990000")
  
  g2<-ggplot(data=dta_ap,aes((speed),height,group=group,color=group))+ geom_point(size=3) + geom_errorbarh(aes(xmin=(speed)-error, xmax=(speed)), size=1.25, height=0)+scale_color_manual(values=c(A = "#000000", B = "#990000"))
  g2<-g2 + scale_x_continuous(limits = c(0,maxError), breaks=seq(0,maxError,(maxError/5)), labels=seq(0,maxError,(maxError/5)),"Height growth (cm.year-1)")
  g2<-g2 + scale_y_continuous(limits = c(0,maxX), breaks=c(0,dta$Height.cm), labels=c(0,dta$Height.cm),"Stem height (cm)",position = "right",
                              sec.axis=sec_axis(~.,breaks=c(0,dta$Height.cm),labels=calendarYears))
  g2<-g2 + theme(legend.position = "none")
  g2
  
  
  plot_grid(g1, g2, ncol = 2, nrow = 1)
}
