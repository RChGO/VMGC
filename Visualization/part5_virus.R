
#病毒物种==================================================================================================================
                
ckv = read.csv('votu.ckv',sep='\t')
ckv$id = sprintf("v%04x", as.hexmode(1:nrow(ckv)))
summary(ckv$contig_length)
sum(ckv$V2>200000)

vtax = read.csv('votu.tax_family',sep='\t',header = F)
ckv$tax = vtax$V4[match(ckv$contig_id,vtax$V1)]
ckv$tax[is.na(ckv$tax)] = 'Unclassified'
#write.table(ckv,'aa',row.names = F,sep='\t',quote = F)

dat = plyr::count(ckv$tax)
dat = dat[order(-dat$freq),]
dat$x[7:nrow(dat)] = 'Others'
dat = aggregate(dat$freq,list(dat$x),sum)
dat = dat[order(-dat$x),]
dat$Group.1 = factor(dat$Group.1,c(setdiff(dat$Group.1,'Others'),'Others'))

ggplot(dat)+
  geom_bar(aes(x=1,y=x,fill=Group.1),stat = 'identity')+
  coord_polar(theta = "y")+
  scale_fill_manual(values = c(mycolor,mycolor))+
  theme_test()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )



#votu稀释曲线=====================================================

library(parallel)
Threads = 50
func = function(y){
  dat = sample(1:nrow(pf),size = y,replace = F)
  dat = data.frame(Step=y,idn=x,num=length(unique(pf$V1[dat])),
                   stringsAsFactors = F)
  return(dat)
}


p = read.csv('virus.i95_c85.cluster',sep='\t',stringsAsFactors = F,header = F)
p = tidyr::separate_rows(p[,c(1,3)],V3,sep=',')
colnames(p) = c('V1','V2')
p$V1 = as.numeric(as.factor(p$V1))

aa = plyr::count(p$V1)
Tab = c()
for(x in 1:2){
  pf = p[p$V1%in%aa$x[aa$freq>=x],]
  Step = rep(seq(1,nrow(pf),round(nrow(pf)/500,0)),each=20)
  cl <- makeCluster(Threads) # 初始化四核心集群
  clusterExport(cl, c("pf","x","prof")) #设定全局变量
  results <- parLapply(cl,Step,func) # lapply的并行版本
  dat <- do.call('rbind',results) # 整合结果
  Tab <- rbind(Tab,dat)
  stopCluster(cl) # 关闭集群
}


#Raref = Tab

aa = aggregate(Tab$num,Tab[,1:2],mean)
aa$sd = aggregate(Tab$num,Tab[,1:2],sd)$x

ggplot(aa)+
  geom_line(aes(x=Step,y=x,group=idn))+
  geom_errorbar(aes(x=Step,y=x,group=idn,ymin=x-sd,ymax=x+sd))+
  scale_x_continuous(breaks = c(0,4000,8000,12000))+
  theme_bw()


#votu与数据库比较======================================================================================================================

p = read.csv('dbcomp.i95_c85.uniq.list',sep='\t',header = F,stringsAsFactors = F)
m = read.csv('dbcomp.map',sep='\t',header=F,stringsAsFactors = F)
p$group = m$V2[match(p$V2,m$V1)]
p$val = 1
#p = p[!(p$V2=='MG1145.k121_3894' & p$group=='votu'),]
dat = dcast(p[,c(1,3,4)],V1~group)

sum(dat$votu!=0)
sum(dat$votu!=0 & rowSums(dat[,2:6])==0)

sum(dat$votu!=0)
sum(dat$votu!=0 & dat$refseq!=0)


gut=dat$V1[rowSums(dat[,2:4])!=0]
oral=dat$V1[dat$OVD!=0]
Ref=dat$V1[dat$RefSeq!=0]
vmgc=dat$V1[dat$votu!=0]



library(VennDiagram)
FIG = venn.diagram(list(vmgc=vmgc,gut=gut,oral=oral,Ref=Ref), 
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


##votu进化树信息========================================================================================================
multip_vir_host = function(host){
  for(x in 3:8){
    dat = host[,c(1,x)]
    #dat[,2] = sub('_[A-Z] ',' ',dat[,2])
    #dat[,2] = sub('_[A-Z]$','',dat[,2])
    aa = unique(dat)
    aa = plyr::count(aa$V1)
    dat[dat[,1]%in%aa$x[aa$freq>1],2] = 'multiple'
    host[,x] = dat[,2]
  }
  return(unique(host[,c(1,3:8)]))
}

ckv = read.csv('votu.ckv',sep='\t')
votu_clstr = read.csv('virus.i95_c85.cluster', sep='\t', header = F, stringsAsFactors = F)
vtax = read.csv('votu.tax_family',sep='\t', header = F, stringsAsFactors = F)
prof = read.csv('votu.host.tax',sep='\t',header = F)
hostf = multip_vir_host(prof)


ckv$host = hostf$V4[match(ckv$contig_id,hostf$V1)]
ckv$host[is.na(ckv$host)] = 'unk'
ckv$tax = vtax$V4[match(ckv$contig_id,vtax$V1)]
ckv$tax[is.na(ckv$tax)] = 'Unclassified'
ckv$host[ckv$tax%in%c('Papillomaviridae','Virgaviridae','Metaviridae',
                      'Iflaviridae','Totiviridae','Polyomaviridae','Retroviridae')] = 'Eukaryote'
ckv$num = votu_clstr$V2[match(ckv$contig_id,votu_clstr$V1)]

#write.table(ckv,'aa',sep='\t',row.names = F,quote = F)

aa = aggregate(ckv$num,list(ckv$tax),sum)
bb = plyr::count(ckv$tax)
aa$votu = bb$freq[match(aa$Group.1,bb$x)]


aa$Group.1[aa$x<2] = 'Others'
aa = aggregate(aa[,2:3],list(aa$Group.1),sum)
aa = aa[order(-aa$votu),]

aa$Group.1 = factor(aa$Group.1,c(setdiff(aa$Group.1,'Others'),'Others'))
aa = aa[order(aa$Group.1),]
ggplot(aa)+
  geom_point(aes(x=Group.1,y=x),color='#e83823')+
  geom_point(aes(x=Group.1,y=votu),color='#3e90c9')+
  scale_y_log10()


ggplot(aa)+
  geom_point(aes(x=Group.1,y=x),color='#e83823')+
  geom_point(aes(x=Group.1,y=votu),color='#3e90c9')+
  #geom_xspline(aes(y=avg_rb,x=as.numeric(Name)),spline_shape = 0.5)+
  geom_bar(aes(y=max(x),x=Group.1,fill=rep(c('a','b'),100)[1:length(unique(Group.1))]),
           position = 'dodge',stat='identity',width = 1,alpha=0.05)+
  scale_y_continuous(trans = 'log10')+
  guides(fill=F)+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust =.5,
      hjust =1 
    )
  )


#ckv$host2 = hostf$V7[match(ckv$V1,hostf$V1)]
dat = ckv[,c('tax','host')]
dat$tax[!dat$tax%in%aa$Group.1]='Others'
dat = plyr::count(dat)

dat$col=mycolor[as.numeric(as.factor(dat$host))+20]
dat$Fill = paste(dat$tax,dat$host,sep=':')
dat = dat[order(dat$freq),]
dat$Fill = factor(dat$Fill,dat$Fill)
dat$tax = factor(dat$tax,aa$Group.1)

ggplot(dat)+
  geom_bar(aes(x=tax,y=freq,fill=Fill),stat = 'identity',position = 'stack',width=.65)+
  scale_fill_manual(values = dat$col)+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle=90,
      hjust=1,
      vjust=.5
    )
  )

aa = aggregate(dat$freq,list(dat$host),sum)
aa = aa[order(aa$x),]
aa$Group.1 = factor(aa$Group.1,aa$Group.1)
ggplot(aa)+
  geom_bar(aes(x=1,y=x,fill=Group.1),stat = 'identity')+
  coord_polar(theta = "y")+
  scale_fill_manual(values = c(mycolor,mycolor))+
  theme_test()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )


#votu和宿主对应关系========================================================================================
#magClstr = read.csv('00.data/mag.clrt',sep='\t',header = F)
#aa = tidyr::separate_rows(magClstr[,c(1,4)],'V4',sep=',')
dat = unique(prof[,c('V1','V9')])
#dat$sgb = aa$SGB[match(dat$V2,aa$MAG)]
#dat = unique(dat[,c('V1','V3','V8','sgb')])

dat = aggregate(dat[,-1],list(dat$V1),function(x){paste(sort(unique(x)),collapse = ', ')})
dat$x = gsub("([a-z]|[A-Z])_([a-z])", "\\1 \\2", dat$x)
#write.table(dat,'aa',sep = '\t',quote = F,row.names = F)


##votu病毒基因的kegg分布===================================================================================

mkegg = read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t',stringsAsFactors = F)
mkegg = mkegg[mkegg$LevelA!='Organismal Systems',]


Fun = read.csv('votu.kegg.m8',header = F,sep='\t',stringsAsFactors = F)
Fun$V2 = sub('.*\\|','',Fun$V2)

nrow(Fun)/sum(ckv$gene_count)

mkeggf = unique(data.frame(kegg=mkegg$KO,Group=mkegg$LevelA,stringsAsFactors = F))
mkeggf = unique(mkeggf)
Fun = merge(Fun,mkeggf,by.x='V2',by.y='kegg',all.x = T)
Fun$Group[is.na(Fun$Group)] = 'Others'
Fun$Group[Fun$Group=='Viral proteins'] = 'Others'

aa =  plyr::count(Fun$V1)
aa$freq =1/aa$freq 
Fun$val = aa$freq[match(Fun$V1,aa$x)]
dat = aggregate(Fun$val,list(Fun$Group),sum)
dat = dat[order(-dat$x),]
dat$Group.1 = factor(dat$Group.1,c(setdiff(dat$Group.1,'Others'),'Others'))

ggplot(dat, aes(x = "", y = x, fill = Group.1)) +
  geom_bar(stat = "identity",position = 'fill', width = 1,color='white') +    ## width >= 1 时中心的杂点将消失
  #facet_grid(~tax)+
  #geom_text(aes(y = POS, x = "", label = Lab), size = 2)+
  #scale_fill_manual(values=dd$Col)+
  #guides(fill=F)+
  #scale_fill_manual(values=pf$Col)+
  coord_polar(theta = "y")+
  theme_test()+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

Fun_m = Fun[Fun$Group=='Metabolism',]
dat= aggregate(Fun_m$val,list(Fun_m$V2),sum)
dat$anno = sub(' \\[EC.*','',mkegg$Anno[match(dat$Group.1,mkegg$KO)])
dat = dat[head(order(-dat$x),n=50),]

aa = merge(dat,unique(mkegg[,c('KO','LevelC')]),by.x='Group.1',by.y='KO')
dat$L = aa$LevelC[match(dat$Group.1,aa$Group.1)]
dat$Group.1 = factor(dat$Group.1,rev(dat$Group.1))


ggplot(dat, aes(x = x, y = Group.1, fill = L)) +
  geom_bar(stat = "identity",position = 'dodge', width = .8,color='white') +    ## width >= 1 时中心的杂点将消失
  #facet_grid(~tax)+
  #geom_text(aes(y = POS, x = "", label = Lab), size = 2)+
  #scale_fill_manual(values=dd$Col)+
  #guides(fill=F)+
  #scale_fill_manual(values=pf$Col)+
  scale_x_sqrt(breaks=c(50,100,200,300,400,500))+
  theme_classic()

dat$lab=paste(dat$Group.1,dat$anno,sep=': ')
ggplot(dat, aes(x = 1, y = Group.1))+
  geom_text(aes(label=lab))


p1 = read.csv('aa.prof',sep='\t',check.names = F)
p2 = read.csv('bb.prof',sep='\t',check.names = F)

dat = merge(p1,p2,all = T)
#write.table(dat,'aaf',sep='\t',row.names = F,quote = F)



##votu病毒基因的cazy分布========================================================================================================
vtax = read.csv('votu.tax_family',sep='\t',header = F,stringsAsFactors = F)

prof = read.csv('votu.cazy.m8',header = F,sep='\t',stringsAsFactors = F)
prof$V2 = sub('\\|$','',prof$V2)
prof$V2 = sub('\\|\\d+.*','',prof$V2)
prof$V2 = sub('(.*?)\\|','',prof$V2,perl = T)

tab = tidyr::separate_rows(prof[,1:2],"V2",sep='[|]')
aa =  plyr::count(tab$V1)
tab$val = 1/aa$freq[match(tab$V1,aa$x)]
tab$votu = sub('_\\d+$','',tab$V1)
tab$tax = vtax$V4[match(tab$votu,vtax$V1)]
tab$tax[is.na(tab$tax)] = 'Unclassified'

tab = aggregate(tab$val,tab[,c('tax','V2')],sum)

aa = aggregate(tab$x,list(tab$V2),sum)
aa = aa[order(-aa$x),]
Or = aa$Group.1[1:30]
tabf = tab[tab$V2%in%Or,]

aa = aggregate(tabf$x,list(tabf$tax),sum)
tabf$tax[!tabf$tax%in%aa$Group.1[aa$x>1]] = 'Others'
tabf = aggregate(tabf$x,tabf[,c('tax','V2')],sum)
tabf = merge(tabf,data.frame(tax=unique(tabf$tax),Col=brewer.pal(12,'Set1')[1:10]))
tabf = tabf[order(-tabf$x),]
tabf$or = paste(tabf$tax,tabf$V2,sep=':')
tabf$or = factor(tabf$or,tabf$or)
tabf$V2 = factor(tabf$V2,Or)

ggplot(tabf)+
  geom_bar(aes(x=x,y=V2,fill=or),color='black',position='stack',stat='identity',width = .75)+
  scale_fill_manual(values = tabf$Col)+
  guides(fill=F)+
  theme_test()

aa = unique(tabf[,c(1,4)])
aa = aa[order(aa$tax),]
ggplot(aa)+
  geom_bar(aes(x=1,y=tax,fill=tax),color='black',position='stack',stat='identity',width = .75)+
  scale_fill_manual(values = aa$Col)+
  guides(fill=F)+
  theme_test()


#HPV===================================================================================================================
p = read.csv('votu.hpv.tax',sep = '\t',stringsAsFactors = F)

dat = p[,c('vOTU','HPVtype','L1')]
dat$Group = 'L1'

dat = rbind(dat,data.frame(p[,c('vOTU','HPVtype')],L1=p$genome-p$L1,Group='no'))
dat  = dat[dat$L1!=0,]
dat$vOTU = factor(dat$vOTU,rev(p$vOTU))
dat$Group = factor(dat$Group,c('no','L1'))
datf = unique(dat[,c('vOTU','HPVtype')])

ggplot(dat)+
  geom_bar(aes(y=vOTU,x=L1,fill=Group),position = 'stack',stat='identity')+
  geom_text(data=datf,aes(y=vOTU,x=60,label=HPVtype))+
  scale_x_continuous(breaks = c(0,15,30,45,60))







p = read.csv('17.hpv/votu.hpv.tax',sep = '\t',header = F,stringsAsFactors = F)
m = read.csv('../mapping.file/mapping.file',sep='\t',stringsAsFactors = F)

p$Sample = sub('\\..*','',p$V1)
p$a = m$Country[match(p$Sample,m$Sample_alias)]

dat = plyr::count(p[,c('a','V16')])
dat = dat[dat$a%in%c('USA','China'),]
aa = aggregate(dat$freq,list(dat$a),sum)
dat$rate = dat$freq/aa$x[match(dat$a,aa$Group.1)]
aa = aggregate(dat$rate,list(dat$V16),mean)
aa = aa[order(-aa$x),]
dat$V16[!dat$V16%in%aa$Group.1[1:30]] = 'others'
dat = aggregate(dat$rate,dat[,1:2],sum)

ggplot(dat)+
  geom_bar(aes(x=V16,y=x),position = 'dodge',stat = 'identity')+
  facet_grid(a~.)+
  scale_fill_manual(values = c(mycolor,mycolor))+
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )


dat = plyr::count(p[,c('V16')])
dat = dat[order(-dat$freq),]
dat$x = factor(dat$x,dat$x)
ggplot(dat)+
  geom_bar(aes(x=x,y=freq),position = 'dodge',stat = 'identity',width=.7)+
  #scale_fill_manual(values = c(mycolor,mycolor))+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )


dat = p
dat$contigid = sub('_\\d+$','',dat$V1)
dd = votu_clstr
dd$votu = ckv$id[match(dd$V1,ckv$V1)]
dd = tidyr::separate_rows(dd[,c('votu','V3')],'V3',sep=',')
dat$votu = dd$votu[match(dat$contigid,dd$V3)]
dat = plyr::count(dat[,c('votu')])
#write.table(dat,'aa',sep='\t',row.names = F)

