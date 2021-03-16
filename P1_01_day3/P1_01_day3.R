pdf(file='P1_01_day3.pdf',width=4.5,height=4.5);
gstable=read.table('P1_01_day3.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("TCF4","MYC","KRT20","hSMARCC1","mKate","CCR5","LGR5","KRT20scrm","mKatescrm","TCF20")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1_vs_0 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(550.5560309239661,346.457300592677),c(26.718160324251294,14.381246439696026))
targetgene="TCF4"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1041.1986114238534,862.8747863817616),c(480.1172446145763,393.52319803168217),c(66.39058019965474,56.21759971881174),c(46.14954965097951,65.3693019986183),c(183.78855738197103,151.65678063679445))
targetgene="MYC"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(85.01232830443594,91.51702279806563),c(0.0,0.0),c(72.05806875328379,57.524985758784105),c(24.289236658410267,33.99203703928151),c(57.484526758237635,67.98407407856303),c(1120.5434511746603,1065.5196225774782),c(4197.989735795241,4358.825057267868),c(11875.007802296781,8218.228647266293),c(4935.572888988966,4630.76135358212),c(230.74774825489754,270.62891027427975),c(2755.2090782856712,3193.94409565249))
targetgene="KRT20"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(336.0011071080087,434.0521652708255),c(684.1468325452225,509.88055558922275),c(214.55492381595738,217.02608263541276))
targetgene="hSMARCC1"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(127.92331306762742,103.28349715781691),c(72.05806875328379,58.83237179875647),c(105.25335885311117,105.89826923776165),c(394.29527508819336,322.92435187317443),c(0.0,0.0),c(0.0,0.0),c(29.147083990092323,30.06987891936442),c(42.10134354124447,48.37328347897754),c(98.77622907753509,105.89826923776165),c(43.72062598513848,22.22556267953022),c(108.4919237408992,120.27951567745768),c(75.29663364107184,91.51702279806563))
targetgene="mKate"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(572.4163439165353,618.3935969069291),c(63.1520153118667,53.602827638867005),c(46.95919087292652,27.455106839419685),c(11928.444122945282,11135.006902444642),c(131.16187795541543,173.8823433163247),c(0.0,0.0),c(918.9427869098552,946.547492939993),c(67.20022142160174,77.1357763583696),c(371.6253208736771,377.83456555201377),c(216.9838474817984,266.7067521543627))
targetgene="CCR5"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(4015.820460857164,3973.1461754760203),c(28.337442768145312,31.377264959336785),c(99.5858702994821,105.89826923776165),c(203.21994670869924,171.26757123637995),c(331.14325977632666,315.08003563334023),c(74.48699241912482,122.89428775740241),c(189.4560459356001,129.43121795726424),c(72.05806875328379,49.68066951894991),c(3441.784834496735,3656.758753802708),c(448.54123695864297,507.26578350927804))
targetgene="LGR5"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(27.527801546198305,54.91021367883937),c(140.87757261877957,172.5749572763523),c(85.82196952638294,86.28747863817615),c(89.06053441417099,60.13975783872884),c(32.38564887788036,41.83635327911571),c(1349.671916985664,1202.7951567745768),c(225.0802597012685,307.23571939350603),c(164.35716805524282,224.87039887524696),c(2343.911337536591,2308.8437465911984),c(571.6067026945883,857.6452422218721))
targetgene="KRT20scrm"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(64.77129775576071,86.28747863817615),c(44.530267207085494,60.13975783872884),c(256.65626735720184,171.26757123637995),c(365.14819109810105,413.13398863126764),c(68.81950386549576,69.2914601185354),c(13.763900773099152,27.455106839419685),c(35.624213765668394,31.377264959336785),c(41.29170231929746,44.45112535906044),c(51.00739698266156,77.1357763583696),c(27.527801546198305,33.99203703928151))
targetgene="mKatescrm"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(80.96412219470089,115.04997151756821),c(237.22487803047363,224.87039887524696))
targetgene="TCF20"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("KRT20scrm","hSATB1","hSMARCC2","SOX9scrm","SOX9","POLRD2scrm","hSMARCD1","mKatescrm","TCF20","hSMARCA4")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1_vs_0 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(27.527801546198305,54.91021367883937),c(140.87757261877957,172.5749572763523),c(85.82196952638294,86.28747863817615),c(89.06053441417099,60.13975783872884),c(32.38564887788036,41.83635327911571),c(1349.671916985664,1202.7951567745768),c(225.0802597012685,307.23571939350603),c(164.35716805524282,224.87039887524696),c(2343.911337536591,2308.8437465911984),c(571.6067026945883,857.6452422218721))
targetgene="KRT20scrm"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(28.337442768145312,64.06191595864594),c(71.24842753133679,86.28747863817615),c(29.95672521203933,35.299423079253884))
targetgene="hSATB1"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(123.06546573594535,121.58690171743004),c(1768.2564287322675,2077.4364175160895),c(91.48945808001201,126.81644587731951))
targetgene="hSMARCC2"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(454.208725512272,655.0004060261554),c(115.77869473842227,130.7386039972366),c(175.69214516250094,145.11985043693264),c(262.32375591083087,257.5550498745561),c(193.50425204533514,237.94425927497062),c(431.5387712977558,417.05614675118477),c(0.0,0.0),c(9262.295579073781,10150.54521434545),c(584.5609622457405,597.4754202673713),c(7438.173906027171,8133.248554668089))
targetgene="SOX9scrm"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.0,0.0),c(140.06793139683253,105.89826923776165),c(144.11613750656758,150.3493945968221),c(59.91345042407866,61.447143878701205),c(0.0,0.0),c(101.20515274337612,125.50905983734714),c(4294.337041206935,5106.649872132061),c(1399.8696727463785,1679.9910613644904),c(2186.8409404788713,2446.1192807882967),c(1256.5631764617578,1219.7911752942175),c(369.19639720783607,345.1499145527046))
targetgene="SOX9"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(84.20268708248894,86.28747863817615),c(119.82690084815732,183.03404559613125),c(140.06793139683253,156.88632479668394),c(131.97151917736247,159.50109687662865),c(295.51904601065826,355.60900287248353))
targetgene="POLRD2scrm"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(686.5757562110635,750.4395869441381),c(89.87017563611799,74.52100427842487),c(116.58833596036929,167.34541311646285),c(94.72802296780004,112.43519943762348))
targetgene="hSMARCD1"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(64.77129775576071,86.28747863817615),c(44.530267207085494,60.13975783872884),c(256.65626735720184,171.26757123637995),c(365.14819109810105,413.13398863126764),c(68.81950386549576,69.2914601185354),c(13.763900773099152,27.455106839419685),c(35.624213765668394,31.377264959336785),c(41.29170231929746,44.45112535906044),c(51.00739698266156,77.1357763583696),c(27.527801546198305,33.99203703928151))
targetgene="mKatescrm"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(80.96412219470089,115.04997151756821),c(237.22487803047363,224.87039887524696))
targetgene="TCF20"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(157.07039705771973,194.80051995588255),c(55.05560309239661,69.2914601185354),c(30.766366433986338,32.68465099930915))
targetgene="hSMARCA4"
collabel=c("P1_F1.fq","P1_F2.fq")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("P1_01_day3_summary.Rnw");
library(tools);

texi2dvi("P1_01_day3_summary.tex",pdf=TRUE);

