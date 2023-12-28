
#Saccharofermentanales完整度评估=====================================================================
dat = tidyr::separate_rows(magClstr,'MAG',sep = ',')

ckm = mag_info[mag_info$MAG%in%unlist(strsplit(magClstr$MAG[magClstr$order=='Saccharofermentanales'],',')),]
#ckm = mag_info[mag_info$MAG%in%unlist(strsplit(magClstr$V3[magClstr$SGB=='SGB006'],',')),]
ckm$MAG[grepl('^GC',ckm$MAG)] = sub('.\\d$','',ckm$MAG)[grepl('^GC',ckm$MAG)]

aa = tidyr::separate_rows(magClstr,'MAG', sep = ",")
View(aa[aa$MAG%in%ckm$MAG[!is.na(ckm$X16s)],])

busco = read.csv('11.Saccharofermentanales_busco/Saccharofermentanales.busco.ass',sep='\t',header = F,stringsAsFactors = F)
busco$V1[grepl('^GC',busco$V1)] = sub('\\.\\d$','',busco$V1)[grepl('^GC',busco$V1)]
busco = busco[busco$V1%in%ckm$MAG,]
busco$ckm = ckm$Completeness[match(busco$V1,ckm$MAG)]

dat = tidyr::separate_rows(magClstr,MAG,sep=',')
dat$MAG[grepl('^GC',dat$MAG)] = sub('.\\d$','',dat$MAG)[grepl('^GC',dat$MAG)] 
busco$Group = dat$SGB[match(busco$V1,dat$MAG)]
busco$Group[busco$Group%in%c('SGB028','SGB130')] = 'BVAB3' 
busco$Group[busco$Group%in%c('SGB006')] = 'BVAB2' 
busco$Group[busco$Group=='SGB069'] = 'SGB069'
busco$Group[!busco$Group%in%c('BVAB2','BVAB3','SGB069')] = 'Others'
busco$Group[busco$V1=='GCF_000025225'] = 'GCF_000025225'


ggplot(busco)+
  geom_point(aes(V2,ckm,color=Group),shape=21,size=1)+
  #scale_x_continuous(limits = c(30,100))+
  #scale_y_continuous(limits = c(30,100))+
  #guides(color=F)+
  theme_light()


###124个bacterial USCO的存在情况===
LAB=c('SGB006','SGB028','SGB069')
Flist=c('11.Saccharofermentanales_busco/bvab2.missing_busco.stat',
        '11.Saccharofermentanales_busco/bvab3.missing_busco.stat',
        '11.Saccharofermentanales_busco/sgb069.missing_busco.stat')

Tab = c()
for(i in 1:3){
  p = read.csv(Flist[i],sep = '\t',header = F)
  p$OR = 1:nrow(p)
  p$rate = p$V2/max(p$V2)*100
  ggplot(p)+
    geom_point(aes(OR,rate),size=1,alpha=.5)
  
  ff = p$V1[1:14]
  gg = dat$MAG[dat$SGB%in%c(LAB[i])]
  #gg = dat$V3[dat$SGB%in%c('SGB028','SGB139')]
  tab = c()
  for(x in gg){
    if(x %in% c('GCF_002892425','GCF_000733145')){
      pf = read.csv(paste('11.Saccharofermentanales_busco/busco_ass/',x,'/run_bacteria_odb10/missing_busco_list.tsv',sep=''),
                    header = F,stringsAsFactors = F,comment.char = '#')
    }else{
      pf = read.csv(paste('00.old_result/11.Saccharofermentanales_busco/busco_ass/',x,'/run_bacteria_odb10/missing_busco_list.tsv',sep=''),
                    header = F,stringsAsFactors = F,comment.char = '#')
    }
    tab = rbind(tab,data.frame(mag=x,miss=length(pf$V1[!pf$V1%in%ff]),stringsAsFactors = F))
  }
  tab$newComp = (1-tab$miss/(124-length(ff)))*100
  tab$oldComp = busco$V2[match(tab$mag,busco$V1)]
  tab$sgb = LAB[i]
  Tab = rbind(Tab,tab)
}



ggplot(Tab)+
  geom_point(aes(x=oldComp,y=newComp))+
  scale_x_continuous(limits = c(30,100))+
  scale_y_continuous(limits = c(30,100))

dat = melt(Tab[,-2],id=c('mag','sgb'))
dat$variable = factor(dat$variable,c('oldComp','newComp'))
ggplot(dat)+
  geom_boxplot(aes(x=sgb,y=value,fill=variable))+
  geom_jitter(aes(x=sgb,y=value,fill=variable),position = position_jitterdodge(0.4),shape=21,size=.6)



Tab$contam = busco$V3[match(Tab$mag,busco$V1)]
#write.table(Tab,'xx',row.names = F,sep='\t',quote = F)


#BVAB基因组示意图====================================================================================
library(circlize)


#展示基因组 ERR2244509.bin.1


Contig = read.csv('07.BVAB/BVAB2/ERR2244509.bin.1.fasta.len',header = F,stringsAsFactors = F,sep='\t')
colnames(Contig) = c('chr','end')
Contig$start = 0
Contig = Contig[,c(1,3,2)]

GC = read.csv('07.BVAB/BVAB2/ERR2244509.bin.1.fasta.GC',header = F,stringsAsFactors = F,sep='\t')
GC$chr = sub('_sliding:.*','',GC$V1)
GC$start = as.numeric(str_extract(GC$V1,'(?<=:)(\\d+)(?=-)'))
GC$end = as.numeric(str_extract(GC$V1,'(?<=-)(\\d+)(?=$)'))
GC$GC = GC$V2
GC$GC_skew_g = '-'
GC$GC_skew_g[GC$V3 >0] = '+'
GC = GC[,c('chr','start','end','GC','GC_skew_g')]
GC$GC[GC$GC_skew_g=='-'] = -GC$GC[GC$GC_skew_g=='-']
#GC1 = GC; GC1$GC[GC1$GC_skew_g=='-'] = -GC1$GC[GC1$GC_skew_g=='-']; GC1$GC[GC1$GC_skew_g=='+']=0
#GC2 = GC; GC2$GC[GC2$GC_skew_g=='-']=0
#bed_list = list(GC1, GC2)

aa = plyr::count(GC[,c('chr','GC_skew_g')])
aa = dcast(data=aa,chr~GC_skew_g)
aa[is.na(aa)] = 0
aa$Or = aa$`+`/rowSums(aa[,2:3])
aa = aa[order(aa$Or),]

Or = aa$chr
#Or = c(setdiff(aa$chr,'k121_13705'),'k121_13705') #MG649.bin.21

Contig$chr = factor(Contig$chr,Or)




###rRNA/tRNA===


#rRNA/tRNA==

p = read.csv('04.rRNA/ERR2244509.bin.1.out',comment.char = '#',header = F)
aa = t(apply(p, 1, function(x){strsplit(x,'  *')[[1]]}))
anno = data.frame(chr='k121_16771',start=150,end=1682,value1=1,gene='16s rRNA')
anno = rbind(anno,data.frame(chr='k121_16771',start=2060,end=4955,value1=1,gene='23s rRNA'))
anno = rbind(anno,data.frame(chr='k121_16771',start=5064,end=5167,value1=1,gene='5s rRNA'))

p = read.csv('05.tRNA/ERR2244509.bin.1.out',comment.char = '#',header = F,skip = 3)
p = as.data.frame(t(apply(p, 1, function(x){strsplit(x,'\\s+')[[1]]})))
p[,'start'] = apply(p[,3:4],1,function(x){min(as.numeric(x))})
p[,'end'] = apply(p[,3:4],1,function(x){max(as.numeric(x))})
anno = rbind(anno,data.frame(chr=p$V1,start=p$start,end=p$end,value1=1,gene=paste('tRNA (',p[,5],')',sep='')))




#cds==
cds = read.csv('07.BVAB/BVAB2/ERR2244509.bin.1.cds.pos',sep='\t',stringsAsFactors = F,header = F)
cds$V1 = sub('_\\d+$','',cds$V1)
#==

circos.clear()
#genome
circos.par("gap.degree" = 0.5)
circos.genomicInitialize(Contig)
circos.genomicTrackPlotRegion(Contig,ylim=c(0,1),track.height=0.05,
                              #bg.col=rand_color(nrow(Contig)),
                              bg.col='grey90',bg.border = NA
)
#GC%
circos.genomicTrack(GC,track.height=0.08,bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, col = ifelse(value[[1]] > 0, "grey50", "grey50"),border = NA)
                      #circos.genomicLines(region, value,col='#3fafff',type='l')
                      circos.segments(x0=0, x1=max(Contig$end), y0=0, y1=0, lwd=0.6, lty="11", col="blue")
                    })
#cds
circos.genomicTrack(cds,track.height=0.1,bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, col = ifelse(value[[1]] > 0, "#f46d43", "#66c2a5"),border = NA)
                    })

# gene labels
circos.genomicLabels(anno, labels.column=5, cex=.3, col='yellow', line_lwd=0.5, line_col="grey80", 
                     side="inside", connection_height=0.03, labels_height=0.04)




#Clostridia类的功能差异================================================================================================
kw_test2=function(prof,map,Group){
  library(fdrtool)
  #prof=p;map=mp;Group='Group'
  p = prof[,map[,1]]
  tab = data.frame(ID=row.names(p),P_value=NA)
  for(x in 1:nrow(p)){
    m=map[!is.na(p[x,]),]
    ff=plyr::count(m[,Group])
    #if( all(ff$freq>=2) & nrow(ff)==length(unique(map[,Group])) ){
    dat = data.frame(ID=as.factor(m[,Group]),value=as.numeric(p[x,m[,1]]))
    tab$P_value[x]=kruskal.test(value~ID,dat)$p.value
    #}
  }
  tab$q_value=NA ; tab$FDR=NA
  tab$q_value[!is.na(tab$P_value)] = p.adjust(tab$P_value[!is.na(tab$P_value)],method = 'BH')
  tab$FDR[!is.na(tab$P_value)] = fdrtool(tab$P_value[!is.na(tab$P_value)],statistic = 'pvalue',plot = F,verbose = F)$qval
  Mean = aggregate(t(p),list(map[,Group]),function(x){mean(x,na.rm=T)})
  row.names(Mean)=Mean$Group.1; Mean=as.data.frame(t(Mean[,-1]))
  Mean$ID=row.names(Mean)
  tab = merge(Mean,tab,by = 'ID')
  return(tab)
}

dat = mag_info[mag_info$MAG%in%magClstr$V1[magClstr$class=='Clostridia'],]
dat$order = magClstr$order[match(dat$MAG,magClstr$V1)]
ff = plyr::count(dat$order)
ff = ff$x[ff$freq>=3]
dat = dat[dat$order%in%ff,]

p = read.csv('10.sgb_fun/sgb.ko.prof',sep='\t',row.names = 1,stringsAsFactors = F)
mp = dat[,c('MAG','order')]
tt = kw_test2(p,mp,'order')
ttf = tt[which(tt$q_value<0.01),]

#ff= apply(tt[,2:7],1,which.max)
#ttf = tt[ff==4,]

for(x in 1:nrow(ttf)){
  r1 = which.max(ttf[x,2:7])+1
  r2 = which.max(ttf[x,setdiff(2:7,r1)])
  ttf[x,'FD']=ttf[x,names(r1)]/ttf[x,names(r2)]
  ttf[x,'max']=names(r1)
}

write.table(ttf,'AA',sep='\t',row.names = F,quote = F)

#ttff=ttf[ttf$FD>=2,]
row.names(ttff) = ttff$ID
ttff =ttff[,2:7]
Or = matrix_or(ttff)
ttff = melt(as.matrix(ttff))
ggplot(ttff)+
  geom_tile(aes(x=Var1,y=Var2,fill=value))+
  scale_fill_gradientn(trans='sqrt',colors=brewer.pal(11,'Spectral'))


##===

foc_ko=c('K10555','K10556','K10557','K10558','K11216','K11531','K08321')
LuxS='K07173'
pf = p[c(foc_ko,LuxS),]

sum(pf[LuxS,]!=0)

#SGB006, SGB028, and SGB069的功能特性==================================================================
komap = read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t',stringsAsFactors = F)
komap = komap[komap$LevelA!='Organismal Systems',]

Gen = read.csv('10.sgb_fun/sgb.genes.count',sep='\t',header = F,stringsAsFactors = F)

p = read.csv('10.sgb_fun/sgb.kegg.m8',sep='\t',header = F,stringsAsFactors = F)
p$id = sub('\\|.*','',p$V1)
p$ko = sub('.*\\|','',p$V2)
p$sgb = magClstr$SGB[match(p$id,magClstr$V1)]
p$sp = magClstr$species[match(p$id,magClstr$V1)]
p$or = magClstr$order[match(p$id,magClstr$V1)]
aa = plyr::count(p[,c('id','sgb')])
aa$all = Gen$V2[match(aa$id,Gen$V1)]
aa$rate = aa$freq/aa$all

pf= p[p$sgb%in%c('SGB006','SGB028','SGB069'),]
#write.table(pf,'aa',sep='\t',row.names = F,quote = F)

length(unique(pf$ko))

library(VennDiagram)
FIG = venn.diagram(list(SGB006=pf$ko[pf$sgb=='SGB006'],SGB028=pf$ko[pf$sgb=='SGB028'],SGB069=pf$ko[pf$sgb=='SGB069']), 
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

pff = p#[!p$sgb%in%c('SGB001','SGB028','SGB069'),]
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
aa$Group[aa$SGB006==1 & aa$SGB028==0 & aa$SGB069 ==0] = 'SGB006'
aa$Group[aa$SGB006==1 & aa$SGB028==1 & aa$SGB069 ==0] = 'SGB006_SGB028'
aa$Group[aa$SGB006==1 & aa$SGB028==0 & aa$SGB069 ==1] = 'SGB006_SGB069'
aa$Group[aa$SGB006==0 & aa$SGB028==1 & aa$SGB069 ==0] = 'SGB028'
aa$Group[aa$SGB006==0 & aa$SGB028==1 & aa$SGB069 ==1] = 'SGB028_SGB069'
aa$Group[aa$SGB006==0 & aa$SGB028==0 & aa$SGB069 ==1] = 'SGB069'
#write.table(aa,'aa',sep='\t',row.names = F,quote = F)

ff = merge(ff,unique(komap[,c('KO','LevelB')]),by.x='ko',by.y='KO')
ff = plyr::count(ff$LevelB)

#AI-2信号系统基因分布==================================================================================
ff = c('K07173','K10558','K10555','K10556','K10557','K08321','K11530','K11216','K11531')
pf = plyr::count(unique(p[p$ko%in%ff,c('ko','sgb')])$ko)
pf$rate=pf$freq/658

pf$x = factor(pf$x,ff)
ggplot(pf)+
  geom_bar(aes(x=x,y=rate),stat = 'identity',width = .7)+
  theme_light()
