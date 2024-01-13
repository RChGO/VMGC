#SGB009, SGB034, and SGB080的功能特性(ipath pathway)==================================================================
komap = read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t',stringsAsFactors = F)
komap = komap[komap$LevelA!='Organismal Systems',]
sgb  =read.csv('sgb.info',sep='\t',stringsAsFactors = F)

Gen = read.csv('sgb.gene_num',sep='\t',header = F,stringsAsFactors = F)

p = read.csv('sgb.kegg.m8',sep='\t',header = F,stringsAsFactors = F)
p$id = sub('\\|.*','',p$V1)
p$ko = sub('.*\\|','',p$V2)
p$sgb = sgb$SGB_ID[match(p$id,sgb$Representative_genome_ID)]
p$sp = sgb$Species[match(p$id,sgb$Representative_genome_ID)]
p$or = sgb$Order[match(p$id,sgb$Representative_genome_ID)]
aa = plyr::count(p[,c('id','sgb')])
aa$all = Gen$V2[match(aa$id,Gen$V1)]
aa$rate = aa$freq/aa$all

pf= p[p$sgb%in%c('SGB009','SGB034','SGB080'),]
#write.table(pf,'aa',sep='\t',row.names = F,quote = F)

length(unique(pf$ko))

library(VennDiagram)
FIG = venn.diagram(list(SGB009=pf$ko[pf$sgb=='SGB009'],SGB034=pf$ko[pf$sgb=='SGB034'],SGB080=pf$ko[pf$sgb=='SGB080']), 
                   fill=c("#729ECE","#FF9E4A","#67BF5C"), 
                   alpha=c(0.5,0.5,0.5), 
                   col = "#414451",
                   cex=1,
                   cat.pos = 6,
                   cat.dist = 0.05,
                   cat.fontface=4, 
                   fontfamily=1,
                   force.unique = T,
                   filename = NULL,
                   height = 300, 
                   width = 300)
grid.draw(FIG)




ff = plyr::count(unique(pf[,c('ko','sgb')])$ko)
ff = pf[pf$ko%in%ff$x[ff$freq==3],]

#ff$anno = komap$Anno[match(ff$ko,komap$KO)]

bb=plyr::count(ff[,c('ko','sgb')])
cc=aggregate(bb$freq,list(bb$ko),sum)

pff = p#[!p$sgb%in%c('SGB001','SGB034','SGB080'),]
for(x in 1:nrow(cc)){
  cc$ra[x] = length(unique(pff$sgb[pff$ko==cc$Group.1[x]]))/632
}

sum(cc$ra<0.3)
plyr::count(unique(p[p$ko=='K00844',14:17])$or)

aa = unique(pf[,c('ko','sgb')])
aa$freq = 1
aa = dcast(aa,ko~sgb)
aa[is.na(aa)] = 0
aa$Group = 'all'
aa$Group[aa$SGB009==1 & aa$SGB034==0 & aa$SGB080 ==0] = 'SGB009'
aa$Group[aa$SGB009==1 & aa$SGB034==1 & aa$SGB080 ==0] = 'SGB009_SGB034'
aa$Group[aa$SGB009==1 & aa$SGB034==0 & aa$SGB080 ==1] = 'SGB009_SGB080'
aa$Group[aa$SGB009==0 & aa$SGB034==1 & aa$SGB080 ==0] = 'SGB034'
aa$Group[aa$SGB009==0 & aa$SGB034==1 & aa$SGB080 ==1] = 'SGB034_SGB080'
aa$Group[aa$SGB009==0 & aa$SGB034==0 & aa$SGB080 ==1] = 'SGB080'
#write.table(aa,'aa',sep='\t',row.names = F,quote = F)

ff = merge(ff,unique(komap[,c('KO','LevelB')]),by.x='ko',by.y='KO')
ff = plyr::count(ff$LevelB)

#AI-2信号系统基因分布==================================================================================
ff = c('K07173','K10558','K10555','K10556','K10557','K08321','K11530','K11216','K11531')
ff2 = c('LuxS', 'lsrA', 'lsrB', 'lsrC', 'lsrD', 'lsrF', 'lsrG', 'lsrK', 'lsrR')
pf = plyr::count(unique(p[p$ko%in%ff,c('ko','sgb')])$ko)
pf$rate=pf$freq/nrow(magClstr)

pf$x = factor(pf$x,ff)
pf = pf[order(pf$x),]
pf$x2 = ff2
pf$x2 = factor(pf$x2,ff2)
ggplot(pf)+
  geom_bar(aes(x=x2,y=rate),stat = 'identity',width = .7)+
  theme_light()


#Saccharofermentanales 完整性=====================================================================
p1 = read.csv('19.Saccharofermentanales/01.busco.stat',sep = '\t',stringsAsFactors = F,header = F)
p2 = read.csv('19.Saccharofermentanales/02.checkm1/checkm1.out',sep = '\t',stringsAsFactors = F)
p3 = read.csv('VMGC_prokaryote_MAG.info',sep = '\t',stringsAsFactors = F)
p3 = p3[p3$Genome_ID%in%p1$V1,]

aa = data.frame(mag=p3$Genome_ID,sgb=p3$Species.level_genomic_bin_.95._ANI.,ckm2=p3$X._Completeness,
           ckm1=p2$Completeness[match(p3$Genome_ID,p2$Bin.Id)],busco=p1$V2[match(p3$Genome_ID,p1$V1)])
aa = melt(aa[,2:5],id.vars = 'sgb')
aa$variable = factor(aa$variable,c('busco','ckm1','ckm2'))
ggplot(aa)+
  geom_boxplot(aes(x=sgb,y=value,fill=variable),outlier.size = .3)


dat = p3[p3$Species.level_genomic_bin_.95._ANI.%in%c('SGB009','SGB034','SGB080'),]

ggplot(dat)+
  geom_point(aes(y=X._Completeness,x=Genome_size_.bp.,color=Genome_quality),shape=21,size=.8)+
  facet_grid(~Species.level_genomic_bin_.95._ANI.,scales = 'free')+
  theme_bw()


#Saccharofermentanales USCO分布=====================================================================
USCOs = sub('.out','',list.files('19.Saccharofermentanales/01.busco/ERR10897634.mbin.2/run_bacteria_odb10/hmmer_output/','*.out'))

#gg = c('SGB009','SGB034','SGB080','Others')
p = read.csv('19.Saccharofermentanales/01.busco.missgene',sep = '\t',stringsAsFactors = F,header = F)
mag = read.csv('VMGC_prokaryote_MAG.info',sep='\t',stringsAsFactors = F)
mag = data.frame(mag = mag$Genome_ID,sgb=mag$Species.level_genomic_bin_.95._ANI.)
#mag = mag[mag$sgb%in%gg,]
mag = mag[mag$mag%in%p$V1,]
#mag$sgb[!mag$sgb%in%gg]='Others'

p = p[p$V1%in%mag$mag,]
p$sgb = mag$sgb[match(p$V1,mag$mag)]
p = plyr::count(p[,2:3])

gg = unique(p$sgb)
tab = c()
fig = c()
for(x in gg){
  pf = p[p$sgb==x,]
  pf = data.frame(id=USCOs,rate=pf$freq[match(USCOs,pf$V2)]/sum(mag$sgb==x)*100)
  pf$rate[is.na(pf$rate)] = 0 
  pf = pf[order(-pf$rate),]
  pf$sgb = x
  pf$x = 1:nrow(pf)
  tab = rbind(tab,pf)
  
  fig[[x]] = ggplot(pf)+
    geom_point(aes(x,rate),shape=21,size=.3)+
    theme_light()
}

ggarrange(plotlist = fig)

tab = dcast(id~sgb,data=tab[,c(1,3,2)])

#write.table(tab,'xx',row.names = F,sep='\t',quote = F)


dat = melt(tab,id.vars = 'id')
dat$lab = 'zz'
dat$lab[dat$value>=90] = 'yy'
dat=dat[order(dat$lab),]
dat$id = factor(dat$id,unique(dat$id))
ggplot(dat)+
  geom_tile(aes(x=id,y=variable,fill=lab))+
  theme_test()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = .5
    )
  )


#Sac进化树marker挑选======================================================================================
p = read.csv('19.Saccharofermentanales_pub/03.marker/all.clstr_cluster.tsv',header=F,sep='\t')
mp = read.csv('19.Saccharofermentanales_pub/03.marker/all.faa.len',header=F,sep='\t')
p$sample = sub('[|].*','',p$V2)
p$sample[grep('^SGB',p$sample)] = sub('[.].*','',p$sample[grep('^SGB',p$sample)])
p$len = mp$V2[match(p$V2,mp$V1)]


#基因出现率
aa = plyr::count(unique(p[,c(1,3)])$V1)
aa$prev = aa$freq/length(unique(p$sample))
#拷贝数
bb = plyr::count(p[,c(1,3)])
bb = aggregate(bb$freq,list(bb$V1),mean)
#序列长度CV
cc = aggregate(p$len,list(p$V1),function(x)(sd(x)/mean(x)))
#
p$prev = aa$prev[match(p$V1,aa$x)]
p$single_rate = bb$x[match(p$V1,bb$Group.1)]
p$len_cv = cc$x[match(p$V1,cc$Group.1)]
p$len_cv[is.na(p$len_cv)]=0

#length(unique(p$V1[p$prev>0.8 & p$single_rate<1.2 & p$len_cv<0.2]))
p = p[p$V1%in%unique(p$V1[p$prev>0.8 & p$single_rate<1.2 & p$len_cv<0.2]),]

#多拷贝，则保留序列长度最接近平均值的那条
tab = c()
for(x in split(p,p$V1)){
  pf = x
  aa = plyr::count(pf$sample)
  bb = mean(pf$len)
  dat = pf[pf$sample%in%aa$x[aa$freq==1],]
  pf = pf[pf$sample%in%aa$x[aa$freq==2],]
  pf = pf[order(abs(pf$len-bb)),]
  dat = rbind(dat,pf[!duplicated(pf$sample),])
  tab = rbind(tab,dat)
}

tab = tab[order(tab$sample,tab$V1),]

#write.table(tab,'19.Saccharofermentanales_pub/03.marker/marker.list',sep='\t',col.names=F,row.names = F,quote = F)

