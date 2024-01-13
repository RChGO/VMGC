#基因集==================================================================================================



#稀释曲线===
gg = c('16.geneSet/pro_i95_c80_cluster.tsv','16.geneSet/pro_i90_c80_cluster.tsv','16.geneSet/pro_i50_c80_cluster.tsv')

library(parallel)
Threads = 30
func = function(y){
  dat = sample(1:nrow(p),size = y,replace = F)
  dat = data.frame(Step=y,idn=x,num=length(unique(p$V1[dat])),nosingle=sum(plyr::count(p$V1[dat])$freq>1),stringsAsFactors = F)
  return(dat)
}

Tab = c()
for(x in gg){
  p = read.csv(x,sep='\t',stringsAsFactors = F,header = F)
  #p = tidyr::separate_rows(p[,c(1,3)],V3,sep=',')
  #colnames(p) = c('V1','V2')
  p$V1 = as.numeric(as.factor(p$V1))
  
  Step = rep(seq(1,nrow(p),round(nrow(p)/50,0)),each=10)
  for(y in Step){
    dat = sample(1:nrow(p),size = y,replace = F)
    Tab = rbind(Tab,data.frame(Step=y,idn=x,num=length(unique(p$V1[dat])),nosingle=sum(plyr::count(p$V1[dat])$freq>1),stringsAsFactors = F))
    print(y)
  }
}

#Raref = Tab

aa = aggregate(Tab$num,Tab[,1:2],mean)
aa$sd = aggregate(Tab$num,Tab[,1:2],sd)$x
aa$nos = aggregate(Tab$nosingle,Tab[,1:2],mean)$x
aa$nos_sd = aggregate(Tab$nosingle,Tab[,1:2],sd)$x

f1=ggplot(aa)+
  geom_line(aes(x=Step,y=x,group=idn))+
  geom_errorbar(aes(x=Step,y=x,group=idn,ymin=x-sd,ymax=x+sd))+
  #geom_line(aes(x=Step,y=nos,group=idn))+
  #geom_errorbar(aes(x=Step,y=nos,group=idn,ymin=nos-nos_sd,ymax=nos+nos_sd))+
  #scale_x_continuous(breaks = c(0,4000,8000,12000))+
  theme_bw()

pdf(file = 'aa.pdf')
f1
dev.off()

#与virgo交集===
p = read.csv('16.geneSet/pro_virgo_i90_c80_cluster.tsv',sep='\t',stringsAsFactors = F,header = F)
#p = read.csv('17.VIRGO_ncbi_check/03.clrt/vmgc_virgo_i90_c80_cluster.tsv',sep='\t',stringsAsFactors = F,header = F)
#p = read.csv('17.VIRGO_ncbi_check/03.clrt/vmgc_virgo_66iso_i90_c80_cluster.tsv',sep='\t',stringsAsFactors = F,header = F)

dat = p
dat$g = 'VMGC'
dat$g[grepl('^V\\d+',dat$V2)] = 'VIRGO'
plyr::count(dat$g)

# dat = p
# dat$g = 'VIRGO-90'
# dat$g[grepl('^vmgc',dat$V2)] = 'VMGC-90'
# dat$g[grepl('^virgo',dat$V2)] = '66iso-90'


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

#交集与非交集丰度分布===
aa = split(dat$V1,dat$g)
Inter = intersect(aa$VIRGO,aa$VMGC)
setd = setdiff(aa$VIRGO,aa$VMGC)

dat = p
dat$g = 'VMGC'
dat$g[grepl('^V\\d+',dat$V2)] = 'VIRGO'

Inter = dat$V2[dat$V1%in%Inter & dat$g=='VIRGO']
setd = dat$V2[dat$V1%in%setd & dat$g=='VIRGO']

clr = read.csv('../data_generation/08.geneSet/VIRGO_i90_c80_cluster.tsv',sep='\t',stringsAsFactors = F,header = F)
Inter = clr$V2[clr$V1%in%Inter]
setd = clr$V2[clr$V1%in%setd]

map = read.csv('../mapping.file/mapping.file20231215',stringsAsFactors = F,sep='\t')
map$Country[map$Country=='USA'] = 'United States'
mapf = split(map,cut(1:nrow(map), breaks = 30, labels = FALSE))

func = function(y){
  map = mapf[[y]]
  #tab = c()
  for(x in 1:nrow(map)){
    p = read.csv(paste('22.virgo_map/rcv/',map$Sample_alias[x],'.rc',sep=''),sep='\t',header = F,stringsAsFactors = F)
    colnames(p)[2] = map$Sample_alias[x]
    if(x==1){
      pf = p
    }else{
      pf = merge(pf,p,by = 'V1',all = T)
    }
  }
  return(pf)
}
Threads=30
library(parallel)
cl <- makeCluster(Threads) # 初始化四核心集群
clusterExport(cl, varlist = list("mapf"),envir=environment()) #设定全局变量, 多个直接加进list，如list("p","d","a")
results <- parLapply(cl,1:length(mapf),func) # lapply的并行版本
#Tab <- do.call('rbind',results) # 整合结果
stopCluster(cl) # 关闭集群

# Tab = results[[1]]
# for(x in 2:length(results)){
#   Tab = merge(Tab,results[[x]],by = 'V1',all = T)
# }

Tab = c()
n = 0
for(x in results){
  n = n+ncol(x)-1
  x[is.na(x)] = 0 
  x[,-1] = sweep(x[,-1],2,colSums(x[,-1]),'/')
  Tab = rbind(Tab,data.frame(id=x[,1],val=as.numeric(rowSums(x[,-1]))))
}
dat = aggregate(Tab$val,list(Tab$id),sum)
dat$x = dat$x/n

dat$g[dat$Group.1%in%Inter] = 'overlap'
dat$g[dat$Group.1%in%setd] = 'nooverlap'

Inter = data.frame(id=Inter,g='overlap',val=dat$x[match(Inter,dat$Group.1)])
setd = data.frame(id=setd,g='nooverlap',val=dat$x[match(setd,dat$Group.1)])

dat = rbind(Inter,setd)
dat[is.na(dat)] = 1e-11
dat = dat[order(-dat$val),]
dat$id = factor(dat$id,dat$id)

#len = read.csv('22.virgo_map/VIRGO.ffn.len',sep='\t',header = F)
#len = len[len$V2>=1000,]
#dat = dat[dat$id%in%len$V1,]

aggregate(dat$val,list(dat$g),summary)
aggregate(dat$val,list(dat$g),sum)

ggplot(dat)+
  #geom_line(aes(x=val,group=g,color=g),alpha=.5,stat = 'count')+
  geom_histogram(aes(x=val,fill=g),position = 'dodge',bins = 50)+
  #geom_boxplot(aes(x=g,y=val))+
  scale_x_log10()

#病毒、细菌、真菌基因overlap==
p = read.csv('16.geneSet/vmgc.list',sep='\t',stringsAsFactors = F,header = F)
clr = read.csv('16.geneSet/pro_i90_c80_cluster.tsv',sep='\t',stringsAsFactors = F,header = F)

clr$g = p$V2[match(clr$V2,p$V1)]

dat = unique(clr[,c('V1','g')])

FIG = venn.diagram(split(dat$V1,dat$g),
                   fill=c("#729ECE","#FF9E4A","#ABC13E"), 
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


#map rate=======================================================================================================

#人群总map率比较###

map = read.csv('sample.info',stringsAsFactors = F,sep='\t')
map$Country[map$Country=='USA'] = 'United States'

gg1 = c('22.virgo_map.mapRate','18.geneSet_map.mapRate','21.vmgc_map.mapRate')
gg2 = c('VIRGO','VMGC_GeneSet','VMGC_GenomeSet')

dat = c()
for(x in 1:3){
  p = read.csv(gg1[x],sep='\t',header = F)
  p$Group = gg2[x]
  dat = rbind(dat,p)
}
dat = dat[dat$V1%in%map$Sample_alias,]

datf = dat
datf$Group2 = map$Country[match(datf$V1,map$Sample_alias)] 
dat$Group2 = 'All'
dat = rbind(dat,datf)
aa = plyr::count(dat$Group2)
aa$lab = paste(aa$x,' (n=',aa$freq/3,')',sep='')
dat$Group2 = aa$lab[match(dat$Group2,aa$x)]

ggplot(dat,aes(x=Group,y=V2,fill=Group))+
  geom_boxplot(outlier.shape = 21,outlier.fill = NA,outlier.size = .5)+
  #stat_compare_means(comparisons = list(c('VIRGO','VMGC_GeneSet'),c('VMGC_GeneSet','VMGC_GenomeSet'),c('VIRGO','VMGC_GenomeSet')),
  #                   method = 'wilcox')+
  facet_grid(~Group2)+
  theme_bw()

aa = aggregate(dat$V2,dat[,c('Group','Group2')],summary)


mapClean = read.csv('18.geneSet_map.mapRead',header = F,stringsAsFactors = F,sep='\t')
mapClean = mapClean[mapClean$V1%in%map$Sample_alias,]
gg = c('22.virgo_map/rcv/','18.geneSet_map/rcv/','21.vmgc_map/rcv/')

n=3
for(x in gg){
  for(y in 1:nrow(mapClean)){
    p  = read.csv(paste(x,mapClean$V1[y],'.rc',sep=''),sep='\t',header = F)
    mapClean[y,n] = sum(p$V2)/mapClean$V2[y]*100
  }
  n=n+1
}

colnames(mapClean)[3:5] = c('VIRGO','VMGC_GeneSet','VMGC_GenomeSet')

dat = melt(mapClean[,c(1,3:5)])
ggplot(dat)+
  geom_boxplot(aes(x=,y=,fill=),outlier.shape = NA)+
  geom_jitter(aes(x=,y=),width = .5,color='#d5d5d5')


####
genmap = read.csv('18.geneSet_map/all_i100_c100_cluster.map',sep='\t',stringsAsFactors = F,header = F)
genmap = unique(genmap[,c(1,3)])
aa = plyr::count(genmap$V1)
genmap$val = 1/aa$freq[match(genmap$V1,aa$x)]

mapf = split(map,cut(1:nrow(map), breaks = 30, labels = FALSE))

func = function(y){
  map = mapf[[y]]
  tab = c()
  for(x in 1:nrow(map)){
    p = read.csv(paste('18.geneSet_map/rcv/',map$Sample_alias[x],'.rc',sep=''),sep='\t',header = F,stringsAsFactors = F)
    pf = merge(p,genmap,by = 'V1')
    pf = aggregate(pf$V2*pf$val,list(pf$V3),sum)
    pf[nrow(pf)+1,2]=mapClean$V2[mapClean$V1==map$Sample_alias[x]]-sum(pf$x)
    pf[nrow(pf),1]= 'Unmap'
    pf$ID=map$Sample_alias[x]
    tab = rbind(tab,pf)
    #print(x)
  }
  return(tab)
}
Threads=30
library(parallel)
cl <- makeCluster(Threads) # 初始化四核心集群
clusterExport(cl, varlist = list("genmap","mapf","mapClean"),envir=environment()) #设定全局变量, 多个直接加进list，如list("p","d","a")
results <- parLapply(cl,1:length(mapf),func) # lapply的并行版本
Tab <- do.call('rbind',results) # 整合结果
stopCluster(cl) # 关闭集群

Tab$country = map$Country[match(Tab$ID,map$Sample_alias)]
Tab$rate = Tab$x/mapClean$V2[match(Tab$ID,mapClean$V1)]

aa = unique(Tab[,c('ID','country')])
aa = plyr::count(aa$country)
aa$lab = paste(aa$x,' (n=',aa$freq,')',sep='')

dat = aggregate(Tab$rate,Tab[,c('country','Group.1')],mean)
dat = dat[dat$Group.1!='Unmap',]
dat$lab = aa$lab[match(dat$country,aa$x)]

ggplot(dat)+
  geom_bar(aes(x=lab,y=x,fill=Group.1),position = 'stack',stat='identity',width=.7)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


Tab$Group.1 = factor(Tab$Group.1,c('Unmap','virus','fungi','bacteria'))
aa = Tab[Tab$Group.1=='bacteria',]
Tab$ID = factor(Tab$ID,aa$ID[order(aa$rate)])
ggplot(Tab)+
  geom_bar(aes(x=ID,y=rate,fill=Group.1),position = 'stack',stat='identity')+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+
  theme_test()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



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






#基因集功能===================================================================================


mkegg = read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t',stringsAsFactors = F)
mkegg = mkegg[mkegg$LevelA!='Organismal Systems',]

mkeggf = unique(data.frame(KO=mkegg$KO,Group=mkegg$LevelB,stringsAsFactors = F))
aa = plyr::count(mkeggf$KO)
bb = unique(mkeggf$Group)
bb = data.frame(Group=bb,Col=mycolor[1:length(bb)])
bb$g  = mkegg$LevelA[match(bb$Group,mkegg$LevelB)]

Fun = read.csv('16.geneSet/pro_i90_c80_rep_seq.kegg.m8',header = F,sep='\t',stringsAsFactors = F)
Fun$V2 = sub('.*\\|','',Fun$V2)

p = read.csv('16.geneSet/pro_i90_c80_cluster.tsv',header = F,sep='\t',stringsAsFactors = F)
m = read.csv('16.geneSet/vmgc.list',header = F,sep='\t',stringsAsFactors = F)
p$g = m$V2[match(p$V2,m$V1)]
Dat = dcast(V1~g,data=p[,c('V1','g')])

n=c()
fig = c()
fig1 = c()
gg = c("all","Bac","Vir","Fung")
for(x in gg){
  dat = Dat
  dat$KO = Fun$V2[match(dat$V1,Fun$V1)]
  dat$KO[is.na(dat$KO)] = 'unk'
  dat$val = 1/aa$freq[match(dat$KO,aa$x)]
  dat$val[is.na(dat$val)] = 1
  
  if(x!='all'){
    dat = dat[dat[,x]!=0 & as.numeric(apply(dat[,setdiff(gg,c(x,'all'))],1,function(x){sum(x!=0)})) == 0,]
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
