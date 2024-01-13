#sgb功能分布============================================================================================================

#可注释率====

gg = c('sgb.kegg.m8','sgb.cazy.m8','sgb.vfdb.m8','sgb.antibio.m8')

tab = data.frame(ID=gg,val=NA,stringsAsFactors = F)
n=1
for(x in gg){
  p = read.csv(paste('10.sgb_fun/',x,sep=''),stringsAsFactors = F,header = F,sep='\t')
  tab$val[n] = nrow(p)
  n=n+1
}

# ##
# pf = p[,1:2]
# pf$V1 = sub('[|].*','',pf$V1)
# pf$V2 = sub('[|][798].*','',pf$V2)
# pf = pf[pf$V1%in%c('GCF_000164135.1','MG139.bin.2','MG482.bin.4','MG74.bin.3'),]

# dat = plyr::count(pf[,1:2])
# dat = dcast(V1~V2,data=dat)
# row.names(dat)= dat$V1
# dat = dat[,-1]
# dat = t(dat)
# dat[is.na(dat)] = 0
# ##

all = read.csv('sgb.gene_num',stringsAsFactors = F,header = F,sep='\t')
tab$rate = tab$val/sum(all$V2)
tab$ID = factor(tab$ID,gg)
ggplot(tab)+
  geom_bar(aes(x=ID,y=rate),position = 'stack',stat='identity',width = .7)+
  scale_fill_manual(values = c(brewer.pal(12,'Paired'),brewer.pal(8,'Set3')))+
  #guides(fill=F,color=F)+
  scale_y_sqrt()+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust =.5
    )
  )



#sgb功能相对贡献====
p = read.csv('prokaryote_SGB.species.profile',sep='\t',row.names = 1,stringsAsFactors = F)
magClstr = read.csv('sgb.info',sep='\t',stringsAsFactors = F)
row.names(p) = sub(':.*','',row.names(p))
rb = as.data.frame(rowMeans(p))
rb['SGB447',] = 0
colnames(rb) = 'rb'


TAB = c()


##kegg===
kegg = read.csv('sgb.kegg.m8',sep='\t',stringsAsFactors = F,header = F)
kegg$ko = sub('.*\\|','',kegg$V2)
aa = sub('\\|.*','',kegg$V1)
kegg$sgb = magClstr$SGB_ID[match(aa,magClstr$Representative_genome_ID)]
kegg = kegg[,c('V1','sgb','ko')]
kegg$val = rb$rb[match(kegg$sgb,row.names(rb))]

#可注释饼图==
dat = kegg
map = read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t')
map = unique(data.frame(group=map$LevelA,ko=map$KO,stringsAsFactors = F))
aa = plyr::count(map$ko)
dat$val = 1/aa$freq[match(dat$ko,aa$x)]
dat= merge(dat,map,by='ko')

dat = aggregate(dat$val,list(dat$group),sum)
dat = dat[order(-dat$x),]
dat[1+nrow(dat),] = NA 
dat$Group.1[nrow(dat)] = 'unknown'
dat$x[nrow(dat)] = sum(all$V2)-sum(dat$x,na.rm = T)
dat$Group.1 = factor(dat$Group.1,dat$Group.1)
dat$rate = dat$x/sum(dat$x)
dat$Alpha = round(1-(1:nrow(dat))*(0.8/nrow(dat)),2)

f1.1=ggplot(dat,aes(x=1,y=rate,fill=Group.1))+
  geom_bar(aes(alpha=Alpha),stat="identity",position="fill",color='black')+
  geom_text(aes(y=1-cumsum(rate)+rate/2, 
                label=paste(Group.1,round(rate*100,2),sep=': ')),size=3,color='black')+
  #coord_polar(theta="y")+
  scale_alpha_continuous(range = c(0,1))+
  scale_fill_manual(values = rep(c('#8DD3C7'),nrow(dat)))+
  #scale_color_manual(values = Col)+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()

dat$Group.1 = as.character(dat$Group.1)
dat$Group.1[dat$Group.1!='unknown'] = 'assig'
dat = aggregate(dat$rate,list(dat$Group.1),sum)
f1.2=ggplot(dat,aes(x=1,y=x,fill=Group.1))+
  geom_bar(stat="identity",position="fill",color='black')+
  scale_fill_manual(values = c('#8DD3C7','white'))+
  scale_alpha_continuous(range = c(0,1))+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()


#生殖道相关有害物质=
map = read.csv('focus.ko.mapping',sep='\t',stringsAsFactors = F)
keggf = kegg[kegg$ko %in% map$KO,c('sgb','ko','val')]
keggf = merge(keggf,map,by.x='ko',by.y='KO')
keggf$tax = magClstr$Order[match(keggf$sgb,magClstr$SGB_ID)]
keggf$sp = magClstr$Species[match(keggf$sgb,magClstr$SGB_ID)]

dat = aggregate(keggf$val,keggf[,c('tax','Group')],sum)
dat = dat[order(dat$Group,-dat$x),]
dat$Type = 'BT'

TAB = rbind(TAB,dat)

#BT特定功能的种水平分布===
keggf$tax = magClstr$Species[match(keggf$sgb,magClstr$SGB_ID)]
keggf$tax = paste(keggf$sgb,keggf$tax,sep=': ')
keggf$order = magClstr$Order[match(keggf$sgb,magClstr$SGB)]

dat = aggregate(keggf$val,keggf[,c('tax','Group','order')],sum)
dat = dat[order(-dat$x),]


Top=15
gg = c()
for(x in unique(dat$Group)){
  datf = dat[dat$Group==x,]
  datf = datf[1:Top,]
  gg = c(gg,datf$order)
}
gg = plyr::count(gg)
gg = gg[order(-gg$freq),]
gg$col = c(brewer.pal(8,'Set3')[-2],brewer.pal(12,'Paired'),'grey')[1:nrow(gg)]

fig = c()
for(x in sort(unique(dat$Group))){
  datf = dat[dat$Group==x,]
  datf = datf[1:Top,]
  datf$Col = gg$col[match(datf$order,gg$x)]
  datf$tax = factor(datf$tax,datf$tax)
  datf = datf[order(datf$order),]
  datf$order = factor(datf$order,unique(datf$order))
  fig[[x]]=ggplot(datf)+
    geom_bar(aes(x=tax,y=x,fill=order),position = 'dodge',stat = 'identity',width = .7,color='black')+
    scale_fill_manual(values = unique(datf$Col))+
    theme_bw()+
    guides(fill=F)+
    labs(title=x,x='',y='Weighted abundance of genes')+
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    )
}

ggarrange(plotlist = fig,align = 'hv')


gg = gg[order(gg$x),]
ggplot(gg)+
  geom_bar(aes(x=x,y=1,fill=x),position = 'dodge',stat = 'identity',width = .7,color='black')+
  scale_fill_manual(values = unique(gg$col))


#acid===

map = read.csv('GMM.mapping',sep='\t',stringsAsFactors = F)
map = map[map$HL1=='organic acid metabolism' & map$Position=='output' & map$Name!='fucose degradation',]
map$Name[map$GMMs == 'MF090'] = paste('butyrate: ',map$Name,sep='')[map$GMMs == 'MF090']
map$Name[map$GMMs == 'MF089'] = paste('acetate: ',map$Name,sep='')[map$GMMs == 'MF089']
map$Group = map$Name

keggf = kegg[kegg$ko %in% map$KO,c('sgb','ko','val')]
keggf = merge(keggf,map[,c('KO','Group')],by.x='ko',by.y='KO')
keggf$tax = magClstr$Order[match(keggf$sgb,magClstr$SGB_ID)]
#Phy = magClstr$phylum[match(keggf$sgb,magClstr$SGB)]
keggf = keggf[!(grepl('butyrate',keggf$Group) & keggf$tax=='Bacteroidales'),]

dat = aggregate(keggf$val,keggf[,c('tax','Group')],sum)
dat = dat[order(dat$Group,-dat$x),]
dat$Type = 'acid'

TAB = rbind(TAB,dat)



##cazy===

cazy = read.csv('10.sgb_fun/sgb.cazy.m8',sep='\t',stringsAsFactors = F,header = F)
cazy$ca = sub('(.*?)\\|','',cazy$V2)
cazy$ca = sub('\\|\\d.*','',cazy$ca)
cazy$ca = sub('\\|$','',cazy$ca)



#可注释饼图==
dat = cazy
dat$val = 1/(as.numeric(lapply(dat$ca,function(x){str_count(x,'\\|')}))+1)
dat = tidyr::separate_rows(dat[,c('ca','val')],'ca',sep='\\|')
dat$Group = sub('\\d.*|_NC$','',dat$ca)

dat = aggregate(dat$val,list(dat$Group),sum)
dat = dat[order(-dat$x),]
dat[1+nrow(dat),] = NA 
dat$Group.1[nrow(dat)] = 'unknown'
dat$x[nrow(dat)] = sum(all$V2)-sum(dat$x,na.rm = T)
dat$Group.1 = factor(dat$Group.1,dat$Group.1)
dat$rate = dat$x/sum(dat$x)
dat$Alpha = round(1-(1:nrow(dat))*(1/nrow(dat)),2)

f2.1=ggplot(dat,aes(x=1,y=rate,fill=Group.1))+
  geom_bar(aes(alpha=Alpha),stat="identity",position="fill",color='black')+
  geom_text(aes(y=1-cumsum(rate)+rate/2, 
                label=paste(Group.1,round(rate*100,2),sep=': ')),size=3,color='black')+
  #coord_polar(theta="y")+
  scale_fill_manual(values = rep(c('#8DD3C7'),nrow(dat)))+
  #scale_color_manual(values = Col)+
  scale_alpha_continuous(range = c(0,1))+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()

dat$Group.1 = as.character(dat$Group.1)
dat$Group.1[dat$Group.1!='unknown'] = 'assig'
dat = aggregate(dat$rate,list(dat$Group.1),sum)
f2.2=ggplot(dat,aes(x=1,y=x,fill=Group.1))+
  geom_bar(stat="identity",position="fill",color='black')+
  scale_fill_manual(values = c('#8DD3C7','white'))+
  scale_alpha_continuous(range = c(0,1))+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()

#热图
aa = sub('\\|.*','',cazy$V1)
cazy$sgb = magClstr$SGB_ID[match(aa,magClstr$Representative_genome_ID)]
cazy$val = rb$rb[match(cazy$sgb,row.names(rb))]
cazy$val = cazy$val/as.numeric(lapply(cazy$ca, function(x){str_count(x,'\\|')+1}))
cazy = tidyr::separate_rows(cazy,'ca', sep = "\\|")
cazy$Group = sub('\\d.*|_NC$','',cazy$ca)
cazy$tax = magClstr$Order[match(cazy$sgb,magClstr$SGB_ID)]

dat = aggregate(cazy$val,cazy[,c('tax','Group')],sum)
dat = dat[order(dat$Group,-dat$x),]
dat$Type = 'cazy'

TAB = rbind(TAB,dat)


##VFDB===
map = read.csv('/share/data1/Database/VFDB/VFDB_setB_pro.20220614.map',stringsAsFactors = F,header = F,sep='\t')
vf = read.csv('sgb.vfdb.m8',sep='\t',stringsAsFactors = F,header = F)
aa = sub('\\|.*','',vf$V1)
vf$sgb = magClstr$SGB_ID[match(aa,magClstr$Representative_genome_ID)]
vf$Group = sub(' \\(.*','',map$V5[match(vf$V2,map$V1)])
vf$Group[vf$Group=='Biofilm'] = 'Biofilm_vf'
vf$val = rb$rb[match(vf$sgb,row.names(rb))]
vf$tax = magClstr$Order[match(vf$sgb,magClstr$SGB_ID)]



#可注释饼图==
dat = vf
dat$val = 1

dat = aggregate(dat$val,list(dat$Group),sum)
dat = dat[order(-dat$x),]
dat[1+nrow(dat),] = NA 
dat$Group.1[nrow(dat)] = 'unknown'
dat$x[nrow(dat)] = sum(all$V2)-sum(dat$x,na.rm = T)
dat$Group.1 = factor(dat$Group.1,dat$Group.1)
dat$rate = dat$x/sum(dat$x)
dat$Alpha = round(1-(1:nrow(dat))*(1/nrow(dat)),2)

f3.1=ggplot(dat,aes(x=1,y=rate,fill=Group.1))+
  geom_bar(aes(alpha=Alpha),stat="identity",position="fill",color='black')+
  geom_text(aes(y=1-cumsum(rate)+rate/2, 
                label=paste(Group.1,round(rate*100,2),sep=': ')),size=3,color='black')+
  #coord_polar(theta="y")+
  scale_fill_manual(values = rep(c('#8DD3C7'),nrow(dat)))+
  #scale_color_manual(values = Col)+
  scale_alpha_continuous(range = c(0,1))+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()

dat$Group.1 = as.character(dat$Group.1)
dat$Group.1[dat$Group.1!='unknown'] = 'assig'
dat = aggregate(dat$rate,list(dat$Group.1),sum)
f3.2=ggplot(dat,aes(x=1,y=x,fill=Group.1))+
  geom_bar(stat="identity",position="fill",color='black')+
  scale_fill_manual(values = c('#8DD3C7','white'))+
  scale_alpha_continuous(range = c(0,1))+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()

#热图
dat = aggregate(vf$val,vf[,c('Group','tax')],sum)
dat = dat[order(dat$Group,-dat$x),]
dat$Type = 'vfdb'

TAB = rbind(TAB,dat)



#anti===

anti = read.csv('sgb.antibio.m8',sep='\t',stringsAsFactors = F,header = F)
anti$Group = sub('clrt(\\d+?)\\|','',anti$V2)
anti$Group = sub('\\|.*','',anti$Group)
aa = sub('\\|.*','',anti$V1)
anti$sgb = magClstr$SGB_ID[match(aa,magClstr$Representative_genome_ID)]
anti$val = rb$rb[match(anti$sgb,row.names(rb))]
anti$tax = magClstr$Order[match(anti$sgb,magClstr$SGB_ID)]


#可注释饼图==
dat = read.csv('sgb.antibio.m8',sep='\t',stringsAsFactors = F,header = F)
dat$val = 1
dat$Group = as.character(lapply(dat$V2, function(x){strsplit(x,'\\|')[[1]][4]}))

dat = aggregate(dat$val,list(dat$Group),sum)
dat = dat[order(-dat$x),]
dat[1+nrow(dat),] = NA 
dat$Group.1[nrow(dat)] = 'unknown'
dat$x[nrow(dat)] = sum(all$V2)-sum(dat$x,na.rm = T)
dat$Group.1 = factor(dat$Group.1,dat$Group.1)
dat$rate = dat$x/sum(dat$x)
dat$Alpha = round(1-(1:nrow(dat))*(1/nrow(dat)),2)

f4.1=ggplot(dat,aes(x=1,y=rate,fill=Group.1))+
  geom_bar(aes(alpha=Alpha),stat="identity",position="fill",color='black')+
  geom_text(aes(y=1-cumsum(rate)+rate/2, 
                label=paste(Group.1,round(rate*100,2),sep=': ')),size=3,color='black')+
  #coord_polar(theta="y")+
  scale_fill_manual(values = rep(c('#8DD3C7'),nrow(dat)))+
  #scale_color_manual(values = Col)+
  scale_alpha_continuous(range = c(0,1))+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()

dat$Group.1 = as.character(dat$Group.1)
dat$Group.1[dat$Group.1!='unknown'] = 'assig'
dat = aggregate(dat$rate,list(dat$Group.1),sum)
f4.2=ggplot(dat,aes(x=1,y=x,fill=Group.1))+
  geom_bar(stat="identity",position="fill",color='black')+
  scale_fill_manual(values = c('#8DD3C7','white'))+
  scale_alpha_continuous(range = c(0,1))+
  guides(color=F,fill=F,alpha=F)+
  theme_transparent()

ggarrange(f1.1,f1.2,f2.1,f2.2,f3.1,f3.2,f4.1,f4.2)

#热图
dat = aggregate(anti$val,anti[,c('tax','Group')],sum)
dat = dat[order(dat$Group,-dat$x),]
dat$Type = 'anti'

TAB = rbind(TAB,dat)



##绘图预处理scale===
TABf = TAB[!TAB$Group%in%c('Others','Biofilm_vf','Antimicrobial activity/Competitive advantage'),]
dat = dcast(Group~tax,data=TABf[,1:3])
#dat[is.na(dat)] = 0

dat[,-1] = t(apply(dat[,-1],1,function(x){scale(x,center = F,scale = T)}))

dat = melt(dat)
dat$Type = TAB$Type[match(dat$Group,TAB$Group)]

dat$Type = factor(dat$Type,c('BT','acid','cazy','anti','vfdb'))
aa = unique(magClstr[,c('Phylum','Order')])
aa = aa[aa$Order%in%dat$variable,]
aa = aa[order(aa$Phylum,aa$Order),]
dat$variable =factor(dat$variable,aa$Order)

ggplot(dat)+
  geom_tile(aes(x=variable,y=Group,fill=value))+
  facet_grid(Type~.,space = 'free',scales = 'free')+
  #scale_fill_gradientn(colors=brewer.pal(11,'Spectral'))+
  scale_fill_gradientn(trans='sqrt',colors=brewer.pal(9,'YlGnBu'))+
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = .5
    )
  )

#order水平丰度
dat = aggregate(p,list(magClstr$Order[match(row.names(p),magClstr$SGB_ID)]),sum)
dat = data.frame(id=factor(dat$Group.1,aa$Order),val=rowMeans(dat[,-1]))

ggplot(dat)+
  geom_bar(aes(x=id,y=val),stat = 'identity',width = .8)+
  #scale_fill_gradientn(colors=brewer.pal(11,'Spectral'))+
  scale_fill_gradientn(colors=brewer.pal(9,'YlGnBu'))+
  scale_y_continuous(trans = 'sqrt')+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = .5
    )
  )

