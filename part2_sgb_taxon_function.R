
#sgb进化树============================================================================================
library(stringr)
p = read.csv('03.gtdb/gtdbtk.ar122.summary.tsv',sep='\t',stringsAsFactors = F)
p = rbind(p,read.csv('03.gtdb/gtdbtk.bac120.summary.tsv',sep='\t',stringsAsFactors = F))


Clstr = read.csv('00.data/mag.clust',header = F,sep='\t',stringsAsFactors = F)
ff = Clstr[!grepl('GC',Clstr$V3),]

#Clstr = tidyr::separate_rows(dat,variable, sep = ",")

p$Group = 'Uncultured'
p$Group[!p$user_genome%in%ff$V1] = 'Cultured (human vagina)'
p$Group[p$user_genome%in%ff$V1 & !p$user_genome%in%p$user_genome[p$fastani_reference=='N/A']] = 'other sources'


p$tax = str_extract(p$classification,'(?<=f__)(.*?)(?=;)')
aa = plyr::count(p[,c('tax','Group')])
ggplot(aa)+
  geom_bar(aes(x=tax,y=freq,fill=Group),width = .7,
           position = 'stack',stat = 'identity')+
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

#PD多样性===
library(picante)
tree= read.tree('06.tree/sgb_best.contree.f')

p = magClstr
p$Reference_genome[grepl('^GC',p$Reference_genome)] = sub('\\..*','',p$Reference_genome[grepl('^GC',p$Reference_genome)])
d = data.frame(row.names = p$Reference_genome,all=rep(1,nrow(p)),stringsAsFactors = F)

d$new = 0
d$new[row.names(d)%in%p$Reference_genome[p$Cultured=='Uncultured']] = 1

tab = c()
gg = unique(p$phylum[p$Cultured=='Uncultured'])
for(x in gg){
  d_f = d[row.names(d)%in%p$Reference_genome[p$phylum==x],]
  aa = pd(t(d_f), tree, include.root=F)
  tab = rbind(tab,data.frame(id=x,all=aa$PD[1],new=aa$PD[2],new_num=sum(d_f$new),stringsAsFactors = F))
}

tab$rate = tab$new/tab$all*100
tab = tab[!is.na(tab$new),]
tab$x = paste(tab$id,' (n=',tab$new_num,')',sep='')
ggplot(tab)+
  geom_bar(aes(x=rate,y=x),stat='identity',width = .75)+
  theme_classic()


p$Cultured[p$Group=='female_reproductive_tract' & p$Cultured=='Cultured']= paste(p$Cultured,'rt',sep=':')[p$Group=='female_reproductive_tract' & p$Cultured=='Cultured']
aa = plyr::count(p[,c('Cultured','phylum')])
bb = plyr::count(p[,c('phylum')])
bb$id = paste(bb$x,' (n=',bb$freq,')',sep='')
aa$phylum = bb$id[match(aa$phylum,bb$x)]
aa = aa[!aa$phylum%in%c('Firmicutes_B (n=1)','Methanobacteriota (n=1)','Synergistota (n=2)','Verrucomicrobiota (n=1)'),]
ggplot(aa)+
  geom_bar(aes(x=phylum,y=freq,fill=Cultured),position ='fill',stat='identity',width = .75)+
  theme_classic()

View(dcast(phylum~Cultured,data=aa))

#细菌丰度表================================================================================
mag_info  = read.csv('00.data/mag.info',sep='\t',stringsAsFactors = F)
magClstr = read.csv('00.data/sgb.info',sep='\t',stringsAsFactors = F)
prof = read.csv('13.bac_prof.raw_prof',sep='\t',row.names = 1)

minR = 400000
fpath = read.csv('13.bac_prof.list.cleanRead_gt_400k',sep='\t',header = F,stringsAsFactors = F)

rc = prof[,colnames(prof)%in%fpath$V1]
rc = RandomMatrix(rc,minR)

aa = data.frame(ID=row.names(rc),rc)
aa$ID[aa$ID=='GCA_011082265.1'] = 'MG347.bin.1'
aa$ID = magClstr$SGB[match(aa$ID,magClstr$Reference.genome)]
#write.table(aa,'13.bac_prof.raw_prof.400k',sep='\t',row.names = F,quote = F)

Len = data.frame(MAG=c(mag_info$MAG,'GCA_011082265.1'),totalLength=c(mag_info$totalLength,1649642))
rb = sweep(rc,1,Len$totalLength[match(row.names(rc),Len$MAG)],'/')
rb = sweep(rb,2,colSums(rb),'/')
#row.names(rb)= paste(magClstr$SGB,magClstr$species,sep=':')[match(row.names(rb),magClstr$V1)]


aa = data.frame(ID=row.names(rb),rb)
aa$ID[aa$ID=='GCA_011082265.1'] = 'MG347.bin.1'
aa$ID = magClstr$SGB[match(aa$ID,magClstr$Reference.genome)]
aa = aa[order(aa$ID),]
#write.table(aa,'13.bac_prof.sgb.tpm',sep='\t',row.names = F,quote = F)


#各水平丰度表
p = read.csv('13.bac_prof.sgb.tpm',sep='\t',stringsAsFactors = F,row.names = 1)
tax = magClstr$genus[match(row.names(p),magClstr$SGB)]
p = aggregate(p,list(tax),sum)

Top=10

aa = data.frame(x=p$Group.1,val=apply(p[,-1],1,mean))
aa = aa[order(-aa$val),]
aa$x[(Top+1):nrow(aa)] = 'Others'
or = c(aa$x[1:Top],'Others')
aa = aggregate(aa$val,list(aa$x),sum)
aa$Group.1 = factor(aa$Group.1,or)
aa = aa[order(aa$Group.1),]
ggplot(aa,aes(x=1,y=x,fill=Group.1))+
  geom_bar(stat="identity",position="fill",alpha=1,color='black')+
  geom_text(aes(y=1-cumsum(x)+x/2, 
                label=paste(Group.1,round(x*100,2),sep=', ')),size=3,color='white')+
  coord_polar(theta="y")+
  #scale_fill_manual(values = brewer.pal(12,'Paired'))+
  #scale_color_manual(values =  brewer.pal(12,'Paired'))+
  guides(color=F,fill=F)+
  theme_transparent()


pf = p[p$Group.1%in%aa$Group.1,] 
pf[Top+1,] = NA 
pf$Group.1[nrow(pf)] = 'Others'
pf[nrow(pf),-1] = 1-colSums(pf[,-1],na.rm = T)

row.names(pf) = pf$Group.1
OR = matrix_or(pf[,-1])
pf = melt(pf)

pf$variable = factor(pf$variable,OR$c_order)
#pf$Group.1 = factor(pf$variable,OR$r_order)
pf$Group.1 = factor(pf$Group.1,aa$Group.1)

ggplot(pf,aes(x=variable,y=value,fill=Group.1))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values = brewer.pal(12,'Paired'))+
  scale_color_manual(values =  brewer.pal(12,'Paired'))+
  #guides(color=F,fill=F)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



#top75 SGB 丰度、binning数量排序=========================================================================================

p = read.csv('13.bac_prof.sgb.tpm',sep='\t',row.names = 1,stringsAsFactors = F)
#map = read.csv('../../data/mapping.file/mapping.file',sep='\t')

sum(magClstr$MAG_num[1:100])
#pf=p[paste(magClstr$SGB[1:75],magClstr$species[1:75],sep=':'),]
pf = p[magClstr$SGB[1:100],]
summary(colSums(pf))
summary(apply(pf, 1, function(x){sum(x!=0)/length(x)}))

tab = magClstr[1:100,c('SGB','species','MAG','phylum')]
tab = tidyr::separate_rows(tab,MAG, sep = ",")
tab$sample = sub('.bin.*','',tab$MAG)
tab$Name = sub('_',' ',paste(tab$SGB,tab$species,sep=':'))
tab$type = 'assem'
tab$type[grep('^GC[AF]',tab$MAG)] = 'cult'
tab$avg_rb = as.numeric(rowMeans(p))[match(tab$SGB,row.names(p))]
tab$Name = factor(tab$Name,rev(sort(unique(tab$Name))))

aa = plyr::count(tab[,c('Name','type')])
f1=ggplot(aa)+
  geom_point(aes(y=freq,x=Name,color=type))+
  guides(fill=F,color=F)+
  scale_y_continuous(trans = 'sqrt',breaks = c(0,50,500,1000,1500))+
  theme_classic()

library(ggalt)
aa = unique(tab[,c('Name','avg_rb')])
f2=ggplot(aa)+
  geom_point(aes(y=avg_rb,x=Name),shape=17)+
  geom_xspline(aes(y=avg_rb,x=as.numeric(Name)),spline_shape = 0.5)+
  geom_bar(aes(y=max(aa$avg_rb),x=as.numeric(Name),fill=rep(c('a','b'),100)[1:length(unique(aa$Name))]),
           position = 'dodge',stat='identity',width = 1,alpha=0.05)+
  scale_y_continuous(trans = 'log10')+
  guides(fill=F)+
  theme_classic()

ggarrange(f1,f2,nrow=2,align = 'hv')


#SGB分类旭日图======================================================================================

magClstr = read.csv('00.data/sgb.info',sep='\t',stringsAsFactors = F)
colnames(magClstr)
p = magClstr[,c('phylum','class','order','family','genus','species')]
p$val = 1

sunplot = function(Table,Count,label_show_cutoff=1){
  #Table，每个分类级别一列
  #Count,对应了扇形的大小，考虑有的时候不是从最低分类水平开始画，为保证个数一致
  #label_show_cutoff,用来展示图中字体，考虑到个数太少展示也看不清
  
  library(dplyr)
  library(plotly)
  
  Value = Table[,Count]
  Table = Table[,colnames(Table)!=Count] 
  
  #生成子级和父级的对应表格
  dat = c()
  for(x in 1:ncol(Table)){
    Table[,x] = paste(Table[,x],x,sep='_') #为了防止子级和父级存在相同名字，预处理
    if(x == 1){
      aa = aggregate(Value,list(Table[,x]),sum)
      colnames(aa) = c('ids','values')
      aa$parents = ""
      
    }else{
      aa = aggregate(Value,Table[,c(x-1,x)],sum)
      colnames(aa) = c('parents','ids','values')
    }
    aa = aa[order(aa$ids),]
    aa$labels = sub('_\\d+$','',aa$ids)
    dat = rbind(dat,aa)
  }
  
  dat$labels[dat$values<label_show_cutoff] = "" 
  
  fig <- plot_ly(branchvalues = 'total',spans = c(10, 20))
  fig <- fig %>%
    add_trace(
      ids = dat$ids,
      labels = dat$labels,
      parents = dat$parents,
      values = dat$values,
      colors = dat$Col,
      type = 'sunburst',
      maxdepth = 6,
      sort = F
    )
  fig <- fig %>%
    layout(
      colorway = c(brewer.pal(8,'Set3')[-2],brewer.pal(12,'Paired'))
    )
  return(fig)
  #return(dat)
}

aa = sunplot(p,'val',1)


#借用itol，添加mag数量和培养信息
p = magClstr[,c('phylum','class','order','family','genus','species','Reference_genome')]
p = p[order(p$phylum,p$class,p$order,p$family,p$genus,p$species),]

Or = read.csv('06.tree/IQK9FVdWDW2PDhaTsXZQhg_newick.label_order',header=F,sep='\t',stringsAsFactors = F)
p$Or = Or$V1
p$mag = magClstr$MAG_num[match(p$Reference_genome,magClstr$Reference_genome)]
p$Cult = paste(magClstr$Cultured,magClstr$Group,sep=':')[match(p$Reference_genome,magClstr$Reference_genome)]

write.table(p[,7:10],'aa',sep='\t',quote=F,row.names = F)


#sgb功能分布 ===========================================================================================================

#可注释率====

gg = c('sgb.kegg.m8','sgb.cazy.m8','sgb.vfdb.m8','sgb.antibio.m8')

tab = data.frame(ID=gg,val=NA,stringsAsFactors = F)
n=1
for(x in gg){
  p = read.csv(paste('10.sgb_fun/',x,sep=''),stringsAsFactors = F,header = F,sep='\t')
  tab$val[n] = nrow(p)
  n=n+1
}

all = read.csv('10.sgb_fun/sgb.gene_num',stringsAsFactors = F,header = F,sep='\t')
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
p = read.csv('13.bac_prof.sgb.tpm',sep='\t',row.names = 1,stringsAsFactors = F)
magClstr = read.csv('00.data/sgb.info',sep='\t',stringsAsFactors = F)
rb = as.data.frame(rowMeans(p))
colnames(rb) = 'rb'


TAB = c()


##kegg===
kegg = read.csv('10.sgb_fun/sgb.kegg.m8',sep='\t',stringsAsFactors = F,header = F)
kegg$ko = sub('.*\\|','',kegg$V2)
aa = sub('\\|.*','',kegg$V1)
kegg$sgb = magClstr$SGB[match(aa,magClstr$Reference_genome)]
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


#BT=
map = read.csv('../mapping.file/focus.ko.mapping',sep='\t',stringsAsFactors = F)
keggf = kegg[kegg$ko %in% map$KO,c('sgb','ko','val')]
keggf = merge(keggf,map,by.x='ko',by.y='KO')
keggf$tax = magClstr$order[match(keggf$sgb,magClstr$SGB)]

dat = aggregate(keggf$val,keggf[,c('tax','Group')],sum)
dat = dat[order(dat$Group,-dat$x),]
dat$Type = 'BT'

TAB = rbind(TAB,dat)

#BT特定功能的种水平分布===
keggf$tax = magClstr$species[match(keggf$sgb,magClstr$SGB)]
keggf$tax = paste(keggf$sgb,keggf$tax,sep=': ')
keggf$order = magClstr$order[match(keggf$sgb,magClstr$SGB)]

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

map = read.csv('../mapping.file/GMM.mapping.new',sep='\t',stringsAsFactors = F)
map = map[map$HL1=='organic acid metabolism' & map$Position=='output' & map$Name!='fucose degradation',]
map$Name[map$GMMs == 'MF090'] = paste('butyrate: ',map$Name,sep='')[map$GMMs == 'MF090']
map$Name[map$GMMs == 'MF089'] = paste('acetate: ',map$Name,sep='')[map$GMMs == 'MF089']
map$Group = map$Name

keggf = kegg[kegg$ko %in% map$KO,c('sgb','ko','val')]
keggf = merge(keggf,map[,c('KO','Group')],by.x='ko',by.y='KO')
keggf$tax = magClstr$order[match(keggf$sgb,magClstr$SGB)]
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
cazy$sgb = magClstr$SGB[match(aa,magClstr$Reference_genome)]
cazy$val = rb$rb[match(cazy$sgb,row.names(rb))]
cazy$val = cazy$val/as.numeric(lapply(cazy$ca, function(x){str_count(x,'\\|')+1}))
cazy = tidyr::separate_rows(cazy,'ca', sep = "\\|")
cazy$Group = sub('\\d.*|_NC$','',cazy$ca)
cazy$tax = magClstr$order[match(cazy$sgb,magClstr$SGB)]

dat = aggregate(cazy$val,cazy[,c('tax','Group')],sum)
dat = dat[order(dat$Group,-dat$x),]
dat$Type = 'cazy'

TAB = rbind(TAB,dat)


##VFDB===
map = read.csv('/share/data1/Database/VFDB/VFDB_setB_pro.20220614.map',stringsAsFactors = F,header = F,sep='\t')
vf = read.csv('10.sgb_fun/sgb.vfdb.m8',sep='\t',stringsAsFactors = F,header = F)
aa = sub('\\|.*','',vf$V1)
vf$sgb = magClstr$SGB[match(aa,magClstr$Reference_genome)]
vf$Group = sub(' \\(.*','',map$V5[match(vf$V2,map$V1)])
vf$Group[vf$Group=='Biofilm'] = 'Biofilm_vf'
vf$val = rb$rb[match(vf$sgb,row.names(rb))]
vf$tax = magClstr$order[match(vf$sgb,magClstr$SGB)]



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

anti = read.csv('10.sgb_fun/sgb.antibio.m8',sep='\t',stringsAsFactors = F,header = F)
anti$Group = sub('clrt(\\d+?)\\|','',anti$V2)
anti$Group = sub('\\|.*','',anti$Group)
aa = sub('\\|.*','',anti$V1)
anti$sgb = magClstr$SGB[match(aa,magClstr$Reference_genome)]
anti$val = rb$rb[match(anti$sgb,row.names(rb))]
anti$tax = magClstr$order[match(anti$sgb,magClstr$SGB)]


#可注释饼图==
dat = read.csv('10.sgb_fun/sgb.antibio.m8',sep='\t',stringsAsFactors = F,header = F)
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
aa = unique(magClstr[,c('phylum','order')])
aa = aa[aa$order%in%dat$variable,]
aa = aa[order(aa$phylum,aa$order),]
dat$variable =factor(dat$variable,aa$order)

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
dat = aggregate(p,list(magClstr$order[match(row.names(p),magClstr$SGB)]),sum)
dat = data.frame(id=factor(dat$Group.1,aa$order),val=rowMeans(dat[,-1]))

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


