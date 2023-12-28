
#基因集==================================================================================================


#稀释曲线===
gg = c('08.geneSet/pro_i50_c80_cluster.tsv','08.geneSet/pro_i90_c80_cluster.tsv','08.geneSet/pro_i95_c80_cluster.tsv')

library(parallel)
Threads = 50
func = function(y){
  dat = sample(1:nrow(p),size = y,replace = F)
  dat = data.frame(Step=y,idn=x,num=length(unique(p$V1[dat])),stringsAsFactors = F)
  return(dat)
}

Tab = c()
for(x in gg){
  p = read.csv(x,sep='\t',stringsAsFactors = F,header = F)
  p$V1 = as.numeric(as.factor(p$V1))
  
  Step = rep(seq(1,nrow(p),round(nrow(p)/50,0)),each=5)
  cl <- makeCluster(Threads) # 初始化四核心集群
  clusterExport(cl, c("p","x")) #设定全局变量
  results <- parLapply(cl,Step,func) # lapply的并行版本
  dat <- do.call('rbind',results) # 整合结果
  Tab <- rbind(Tab,dat)
  stopCluster(cl) # 关闭集群
}

#Raref = Tab

aa = aggregate(Tab$num,Tab[,1:2],mean)

ggplot(aa)+
  geom_line(aes(x=Step,y=x,group=idn))+
  theme_bw()


#与virgo交集===

fungi_info = read.csv('22.fungi/fungi.info',sep='\t',stringsAsFactors = F)

p = read.csv('08.geneSet/pro_i90_c80_cluster.tsv',sep='\t',stringsAsFactors = F,header = F)
p$V2 = sub('\\|.*','',p$V2)
p$V2 = sub('_\\d+_g$','',p$V2)
p$V3 = NA
p$V3[p$V2%in%mag_info$MAG] = 'ba'
p$V3[p$V2%in%fungi_info$ID] = 'fu'
p$V3[grepl('^V\\d+$',p$V2)] = 'virgo'
p$V3[is.na(p$V3)] = 'vir'
#p$V4 = as.numeric(as.factor(p$V1))

pf = unique(p[,c('V1','V3')])
FIG = venn.diagram(split(pf$V1,pf$V3),
                   fill=c("#729ECE","#FF9E4A","#67BF5C"), 
                   alpha=c(0.5,0.5,0.5), 
                   col = "#414451",
                   cex=1,
                   cat.pos = 6,
                   cat.dist = 0.05,
                   cat.fontface=4, 
                   #fontfamily=1,
                   force.unique = T,
                   filename = NULL,
                   height = 300, 
                   width = 300)
grid.draw(FIG)




p = read.csv('08.geneSet/pro_VIRGO_cluster.tsv',sep='\t',stringsAsFactors = F,header = F)

dat = p
dat$g = 'VMGC'
dat$g[grepl('^V',dat$V2)] = 'VIRGO'

dat = unique(dat[,c('V1','g')])
FIG = venn.diagram(split(dat$V1,dat$g),
                   fill=c("#729ECE","#FF9E4A"), 
                   alpha=c(0.5,0.5), 
                   col = "#414451",
                   cex=1,
                   cat.pos = 6,
                   cat.dist = 0.05,
                   cat.fontface=4, 
                   #fontfamily=1,
                   force.unique = T,
                   filename = NULL,
                   height = 300, 
                   width = 300)
grid.draw(FIG)



dat = merge(p,pf,by.x='V2',by.y='V1',all=T)
dat$V3[is.na(dat$V3)] = 'virgo'
Dat = dcast(V1~V3,data=dat[,c('V1','V3')])
dat$V4 = as.numeric(as.factor(dat$V1))
dat = unique(dat[,c('V4','V3')])
dat = split(dat$V4,dat$V3)

library(VennDiagram)
FIG = venn.diagram(dat,
                   fill=c("#729ECE","#FF9E4A","#67BF5C","#c49611"), 
                   alpha=c(0.5,0.5,0.5,0.5), 
                   col = "#414451",
                   cex=1,
                   cat.pos = 6,
                   cat.dist = 0.05,
                   cat.fontface=4, 
                   #fontfamily=1,
                   force.unique = T,
                   filename = NULL,
                   height = 300, 
                   width = 300)
grid.draw(FIG)


sum(Dat[rowSums(Dat[,2:4])!=0 & Dat$virgo!=0,2:4])
sum(Dat[rowSums(Dat[,2:4])!=0 & Dat$virgo!=0,5])
sum(Dat[rowSums(Dat[,2:4])!=0 & Dat$virgo==0,2:4])
sum(Dat[rowSums(Dat[,2:4])==0 & Dat$virgo!=0,5])






#功能===

mkegg = read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t',stringsAsFactors = F)
mkegg = mkegg[mkegg$LevelA!='Organismal Systems',]

mkeggf = unique(data.frame(KO=mkegg$KO,Group=mkegg$LevelB,stringsAsFactors = F))
aa = plyr::count(mkeggf$KO)
bb = unique(mkeggf$Group)
bb = data.frame(Group=bb,Col=mycolor[1:length(bb)])
bb$g  = mkegg$LevelA[match(bb$Group,mkegg$LevelB)]

Fun = read.csv('08.geneSet/pro_virgo_i90_c80_rep_seq.kegg.m8',header = F,sep='\t',stringsAsFactors = F)
Fun$V2 = sub('.*\\|','',Fun$V2)

n=c()
fig = c()
fig1 = c()
gg = c('ba','fu','vir','virgo')
for(x in gg){
  dat = Dat
  dat$KO = Fun$V2[match(dat$V1,Fun$V1)]
  dat$KO[is.na(dat$KO)] = 'unk'
  dat$val = 1/aa$freq[match(dat$KO,aa$x)]
  dat$val[is.na(dat$val)] = 1
  
  if(x!='virgo'){
    dat = dat[dat[,x]!=0 & as.numeric(apply(dat[,setdiff(gg,c(x,'virgo'))],1,function(x){sum(x!=0)})) == 0,]
  }else{
    dat = dat[as.numeric(apply(dat[,setdiff(gg,'virgo')],1,function(x){sum(x!=0)})) > 0,]
  }
  dat = data.frame(KO=dat$KO,val=dat$val,stringsAsFactors = F)
  
  ff = data.frame(Group=c('k','u'),val=c(sum(dat$KO!='unk'),sum(dat$KO=='unk')))
  ff$x = ff$val/sum(ff$val)
  fig1[[x]]=ggplot(ff, aes(x = "", y = x, fill = Group)) +
    geom_bar(stat = "identity",position = 'fill', width = 1,color='black') +    ## width >= 1 时中心的杂点将消失
    #facet_grid(~tax)+
    #geom_text(aes(y = POS, x = "", label = Lab), size = 2)+
    geom_text(aes(y=1-cumsum(x)+x/2, 
                  label=round(x*100,2)),size=3,color='black')+
    scale_fill_manual(values=c('#bbbbbb','white'))+
    guides(fill=F)+
    #scale_fill_manual(values=pf$Col)+
    coord_polar(theta = "y")+
    theme_test()+
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  
  dat = merge(dat,mkeggf,by='KO',all.x = T)
  
  dat = aggregate(dat$val,list(dat$Group),sum)
  dat = dat[order(-dat$x),]
  Or = c(setdiff(dat$Group.1,'Others')[1:15],'Others')
  dat$Group.1[!dat$Group.1 %in% Or] = 'Others'
  dat = aggregate(dat$x,list(dat$Group.1),sum)
  dat$Col = bb$Col[match(dat$Group.1,bb$Group)]
  dat$Group.1 = factor(dat$Group.1,Or) 
  dat = dat[order(dat$Group.1),]
  
  fig[[x]]=ggplot(dat, aes(x = "", y = x, fill = Group.1)) +
    geom_bar(stat = "identity",position = 'fill', width = 1,color='white') +    ## width >= 1 时中心的杂点将消失
    #facet_grid(~tax)+
    #geom_text(aes(y = POS, x = "", label = Lab), size = 2)+
    scale_fill_manual(values=dat$Col)+
    guides(fill=F)+
    #scale_fill_manual(values=pf$Col)+
    coord_polar(theta = "y")+
    theme_test()+
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  n = c(n,Or)
}
ggarrange(plotlist = fig1,align = 'hv')
ggarrange(plotlist = fig,align = 'hv')

bbf = bb[bb$Group%in%n,]
bbf = bbf[order(bbf$g,bbf$Group),]
bbf$Group = factor(bbf$Group,bbf$Group)
ggplot(bbf, aes(x = 1, y = Group, fill = Group)) +
  geom_bar(stat = "identity",position = 'fill', width = 1,color='white') +    ## width >= 1 时中心的杂点将消失
  scale_fill_manual(values=bbf$Col)




#map rate=======================================================================================================
cleanread = read.csv('13.bac_prof.list.cleanRead_gt_400k',header = F,stringsAsFactors = F,sep='\t')
cleanread = cleanread[cleanread$V2>1e6,]

mapClean = read.csv('08.geneSet_map_stat.mapRead',header = F,stringsAsFactors = F,sep='\t')
mapClean = mapClean[mapClean$V1%in%cleanread$V1,]

map = read.csv('../mapping.file/mapping.file',stringsAsFactors = F,sep='\t')
map$Country[map$Country=='USA'] = 'United States'

p = read.csv('08.geneSet_map_stat.all',header = F,stringsAsFactors = F,sep='\t')
p = p[p$V1%in%mapClean$V1,]
p[is.na(p)] = 0 
p$rate = p$V3/mapClean$V2[match(p$V1,mapClean$V1)]
p$country = map$Country[match(p$V1,map$Sample_alias)]
p$proj = map$Project_NCBI_BioProject[match(p$V1,map$Sample_alias)]
p$all = mapClean$V2[match(p$V1,mapClean$V1)]

dat = aggregate(p$rate,p[,c('country','V2')],mean)


ggplot(p)+
  geom_boxplot(aes(x=country,y=rate,fill=V2))+
  scale_y_sqrt()

aa = unique(p[,c('V1','country')])
plyr::count(aa$country)

ggplot(dat)+
  geom_bar(aes(x=country,y=x,fill=V2),position = 'stack',stat='identity',width = .7)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

aa = aggregate(p$rate,p[,-c(2,3,4)],sum)
mean(aa$x[aa$proj=='PRJNA695070'])

aa = aggregate(p$rate,p[,c('V1','country')],sum)
aa = plyr::count(aa$country)


bb=dcast(country~V2,data=dat)
bb$sample_num = aa$freq[match(bb$country,aa$x)]
