#样品信息==================================================================================
library(maps)
library(sf)
library(proj)

p = read.csv('sample.info',sep='\t')
dd = plyr::count(p$Country)
dd$x[dd$x=='America'] = 'USA'
dd$x[dd$x=='Kenyan'] = 'Kenya'
#GG=c("China","Germany","France","Fiji","America","USA","Luxembourg","Italy","Australia","Spain")


world_map <- map_data("world")
world_map$Or = 1:nrow(world_map)
world_map <- merge(world_map,dd,by.x='region',by.y='x',all=T)
world_map$freq[is.na(world_map$freq)] <- 1
world_map <- world_map[order(world_map$Or),]

ggplot(world_map, aes(x = long, y = lat, group = group)) +
  #geom_sf()+
  #geom_polygon(fill='grey', colour =NA,size=0.1) +
  geom_polygon(aes(fill= freq),color=NA, size=0.1) +
  #scale_x_continuous(breaks = seq(-180, 210, 45), labels = function(x){paste0(x, "°")}) +
  #scale_y_continuous(breaks = seq(-60, 100, 30), labels = function(x){paste0(x, "°")}) +
  #scale_fill_gradientn(tran='log10',colours = c('#dddddd','#E6F5D0',brewer.pal(11,'RdYlBu')[c(9:11)]))+
  scale_fill_gradientn(tran='log10',colours = c('#dddddd',rev(brewer.pal(11,'Spectral')[c(1:4)])))+
  #scale_fill_gradientn(colors = c('grey'))+
  #scale_fill_gradient(tran='log10',low = "#fee0d2", high="#ef3b2c")+
  #theme_light()
  theme_test()
  #coord_map(projection = 'mollweide')
  #coord_map(projection = 'sinusoidal')
  #coord_map("polyconic")

#mag信息====================================================================================

#Col1 = c('#fdebe0','#fdfce9','#bbeae9')
#Col2 = c('#e85338','#fc9003','#28a5a8')

mag_info  = read.csv('mag.info',sep='\t',stringsAsFactors = F)

#完整性情况
mag_info$Genome_quality = factor(mag_info$Genome_quality,c("near-complete","high-quality","medium-quality"))
f3=ggplot(mag_info)+
  geom_point(aes(x=Completeness,y=Contamination,color=Genome_quality),fill=NA,alpha=.3,size=.5)+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_test()

f1=ggplot(mag_info)+
  geom_histogram(aes(x=Completeness,y=..count..,group=Genome_quality,fill=Genome_quality,color=Genome_quality),bins = 50)+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_classic()


f4=ggplot(mag_info)+
  geom_histogram(aes(y=Contamination,x=..count..,group=Genome_quality,fill=Genome_quality,color=Genome_quality),bins = 50)+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_classic()



aa = plyr::count(mag_info$Genome_quality)
aa$rate = aa$freq/sum(aa$freq)
aa$x = factor(aa$x,c("near-complete","high-quality","medium-quality"))
aa = aa[order(aa$x),] 
f2=ggplot(aa,aes(x=1,y=rate,fill=x,color=x))+
  geom_bar(stat="identity",position="fill",alpha=1)+
  geom_text(aes(y=1-cumsum(rate)+rate/2, 
                label=paste(x,freq,sep='\n')),size=3)+
  coord_polar(theta="y")+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_transparent()


ggarrange(f1,f2,f3,f4,align = 'hv')


##基因组组装情况
f1=ggplot(mag_info)+
  geom_boxplot(aes(x=Genome_quality,y=No._of_contigs,fill=Genome_quality),outlier.shape = 1,outlier.size = .5)+
  #scale_y_sqrt(breaks=c(0,200,400,600))+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_classic()

pf = mag_info[,c('Genome_quality','Genome_size')]
pf = rbind(pf,data.frame(Genome_quality='all',Genome_size=pf$Genome_size))
f2=ggplot(pf)+
  geom_boxplot(aes(x=Genome_quality,y=Genome_size,fill=Genome_quality),outlier.shape = 1,outlier.size = .5)+
  scale_y_log10()+
  scale_fill_manual(values = c(Col,'grey'))+
  scale_color_manual(values = c(Col,'grey'))+
  guides(color=F,fill=F)+
  theme_classic()

pf = mag_info[,c('Genome_quality','N50_length')]
pf = rbind(pf,data.frame(Genome_quality='Genome_quality',N50_length=pf$N50_length))

f3=ggplot(pf)+
  geom_boxplot(aes(x=Genome_quality,y=N50_length,fill=Genome_quality),outlier.shape = 1,outlier.size = .5)+
  scale_y_log10()+
  scale_fill_manual(values = c(Col,'grey'))+
  scale_color_manual(values = c(Col,'grey'))+
  guides(color=F,fill=F)+
  theme_classic()


ggarrange(f2,f3,f2,f3,align = 'hv')


summary(mag_info$N50_length)
summary(mag_info$Genome_size)


#真核基因组信息====================================================================================

fungi_info  = read.csv('fungi.info',sep='\t',stringsAsFactors = F)
aa = head(fungi_info,n=38)

ggplot(fungi_info)+
  geom_point(aes(x=X..Completeness,y=X..Contamination),color='#3E90C9',alpha=.8,size=1)+
  #scale_fill_manual(values = '3E90C9')+
  #scale_color_manual(values = '3E90C9')+
  guides(color=F,fill=F)+
  theme_test()


#病毒基因组信息====================================================================================
vs_info  = read.csv('allvirus.ckv',sep='\t',stringsAsFactors = F,header=F)

aa = plyr::count(vs_info$V8)
aa$rate = aa$freq/sum(aa$freq)
aa$x = factor(aa$x,c("Complete","High-quality","Medium-quality"))
aa = aa[order(aa$x),] 
f1=ggplot(aa,aes(x=1,y=rate,fill=x))+
  geom_bar(stat="identity",position="fill",alpha=1,color='black')+
  geom_text(aes(y=1-cumsum(rate)+rate/2, 
                label=paste(x,freq,sep='\n')),size=3,color='white')+
  coord_polar(theta="y")+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_transparent()

aa = vs_info
aa$g = 'High'
aa$g[aa$V12<20] = 'Med'
aa$g[aa$V12<10] = 'Low'
aa = plyr::count(aa$g)
aa$rate = aa$freq/sum(aa$freq)
aa$x = factor(aa$x,c("High","Med","Low"))
aa = aa[order(aa$x),] 
f2=ggplot(aa,aes(x=1,y=rate,fill=x))+
  geom_bar(stat="identity",position="fill",alpha=1,color='black')+
  geom_text(aes(y=1-cumsum(rate)+rate/2, 
                label=paste(x,freq,sep='\n')),size=3,color='white')+
  coord_polar(theta="y")+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_transparent()


pf=vs_info[,c('V8','V2')]
pf = rbind(pf,data.frame(V8='all',V2=pf$V2))
pf$V8 = factor(pf$V8,c("Complete","High-quality","Medium-quality",'all'))
f3=ggplot(pf)+
  geom_boxplot(aes(x=V8,y=V2,fill=V8),outlier.shape = 1,outlier.size = .5)+
  scale_y_log10()+
  scale_fill_manual(values = c(Col,'grey'))+
  scale_color_manual(values = c(Col,'grey'))+
  guides(color=F,fill=F)+
  theme_classic()

ggarrange(f1,f2,f3)

summary(vs_info$V2)

#SGB分类旭日图======================================================================================

magClstr = read.csv('sgb.info',sep='\t',stringsAsFactors = F)
colnames(magClstr)
p = magClstr[,c('Phylum','Class','Order','Family','Genus')]
p$val = 1

sunplot = function(Table,Count,label_show_cutoff=1,pieSort=F){
  #Table, 每个分类级别一列的表格
  #Count, 表格中，最小分类级别的计数列，将对应了扇形的大小（考虑有的时候不是从最低分类水平开始画，为保证个数一致）
  #label_show_cutoff,用来展示图中字体，考虑到个数太少展示也看不清
  #pieSort,F表示按每个分级内ids的字符顺序展示，T表示按每个分级内扇面的大小顺序展示,若要使用itol添加metadata，建议使用F模式
  
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
      sort = ifelse(pieSort==F,F,T)
    )
  fig <- fig %>%
    layout(
      colorway = c(brewer.pal(8,'Set3')[-2],brewer.pal(12,'Paired'))
    )
  return(fig)
  #return(dat)
}

aa = sunplot(p,'val',2)
aa


####添加mag数量和培养信息==

#按照sunplot的node顺序排序文本信息
p = magClstr[,c('Phylum','Class','Order','Family','Genus','Species','SGB_ID')]
p = p[order(p$Phylum,p$Class,p$Order,p$Family,p$Genus,p$Species),]

#借用gtdbtk和itol生成的进化树
Or = read.csv('sgb_best_newick.lab_order',header=F,sep='\t',stringsAsFactors = F)

#按sunplot的node顺序依次给进化树的node（进化树本身的node顺序不重要）添加metadata
p$Or = Or$V1
p$mag = magClstr$Number_of_genomes[match(p$SGB_ID,magClstr$SGB_ID)]
p$Cult = magClstr$Cultured[match(p$SGB_ID,magClstr$SGB_ID)]

#write.table(p[,7:10],'aa',sep='\t',quote=F,row.names = F)


#PD多样性===
library(picante)
tree= read.tree('prokaryote_SGB.newick')

p = magClstr
#p$Reference_genome[grepl('^GC',p$Reference_genome)] = sub('\\..*','',p$Reference_genome[grepl('^GC',p$Reference_genome)])
d = data.frame(row.names = p$SGB_ID,all=rep(1,nrow(p)),stringsAsFactors = F)

d$new = 0
d$new[row.names(d)%in%p$SGB_ID[p$Cultured=='Uncultured']] = 1

tab = c()
gg = unique(p$Phylum[p$Cultured=='Uncultured'])
for(x in gg){
  d_f = d[row.names(d)%in%p$SGB_ID[p$Phylum==x],]
  aa = pd(t(d_f), tree, include.root=F)
  tab = rbind(tab,data.frame(id=x,all=aa$PD[1],new=aa$PD[2],new_num=sum(d_f$new),stringsAsFactors = F))
}

dat=plyr::count(p$Phylum)
dat = dat[dat$freq>=5,]

tab$rate = tab$new/tab$all*100
#tab = tab[!is.na(tab$new),]
tab = tab[tab$id%in%dat$x,]
tab$x = paste(tab$id,' (n=',dat$freq[match(tab$id,dat$x)],')',sep='')
f1=ggplot(tab)+
  geom_bar(aes(x=rate,y=x),stat='identity',width = .7)+
  theme_classic()

dat=plyr::count(p[,c('Phylum','Cultured')])
dat = dat[dat$Phylum%in%tab$id,]
dat$Cultured = factor(dat$Cultured,c('Uncultured','Cultured (women vagina)','Cultured (others)'))
f2=ggplot(dat)+
  geom_bar(aes(y=Phylum,x=freq,fill=Cultured),stat='identity',position = 'fill',width = .7)+
  theme_classic()

ggarrange(f1,f2,widths = c(1,1.5))


#细菌丰度表=====================================================================================================================

#PRJNA881266项目为VLP测序，计算丰度时，注意踢掉相关样品
flist = list.files('11.bac_prof/rcv2','*.report')
flist = flist[!grepl('bracken',flist)]
map = read.csv('sample.info',sep='\t',stringsAsFactors = F)
flist = flist[flist%in%paste(map$Sample_alias,'.report',sep='')]


gg = c('d','p','c','o','f','g','s')
#gg = c('s')
prof = c()
n=1
for(x in flist){
  p = read.csv(paste('11.bac_prof/rcv2/',x,sep=''),stringsAsFactors = F,header = F,sep='\t')
  p$V6 = tolower(p$V6)
  p$V8 = sub('^ +','',p$V8)
  for(y in gg){
    pf = p[p$V6==y,]
    pf = data.frame(ID=pf$V8,val=pf$V1/sum(pf$V1))
    #pf = data.frame(ID=pf$V8,val=pf$V2)
    NANE = sub('.report$','',x)
    colnames(pf)[2] = NANE
    if(n == 1){
      prof[[y]] = pf
    }else{
      #prof[[y]] = merge(prof[[y]],pf,by='ID',all = T)
      ff = setdiff(pf$ID,prof[[y]]$ID)
      if(length(ff)!=0){
        aa = matrix(nrow = length(ff),ncol = ncol(prof[[y]]))
        aa[,1] = ff
        colnames(aa) = colnames(prof[[y]])
        prof[[y]] = rbind(prof[[y]],aa)
      }
      prof[[y]][,NANE] = pf[match(prof[[y]]$ID,pf$ID),2]
    }
  }
  n=n+1
  print(n)
}

for(x in names(prof)){
  #write.table(prof[[x]],paste('11.bac_prof/prof2/',x,'.prof',sep=''),row.names = F,sep='\t',na = '0',quote = F)
}



#各水平丰度表==================================================================================
p = read.csv('prokaryote_SGB.genus.profile',sep='\t',stringsAsFactors = F)
#p$ID = sub('.*:','',p$ID)
#p = aggregate(p[,-1],list(ID=p$ID),sum)
mp = read.csv('11.bac_prof/rcv.maprate',sep='\t',header = F)
#p = p[,colnames(p)%in%c('ID',mp$V1[mp$V3>=4000000])]
map = read.csv('sample.info',sep='\t',stringsAsFactors = F)
p = p[,colnames(p)%in%c('ID',map$Sample_alias)]

Top=15

aa = data.frame(x=p$ID,val=apply(p[,-1],1,mean))
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


pf = p[p$ID%in%aa$Group.1,] 
pf[Top+1,] = NA 
pf$ID[nrow(pf)] = 'Others'
pf[nrow(pf),-1] = 1-colSums(pf[,-1],na.rm = T)

row.names(pf) = pf$ID
OR = matrix_or(pf[,-1])
pf = melt(pf)
pf$value[pf$value<0] = 0
pf$variable = factor(pf$variable,OR$c_order)
#pf$Group.1 = factor(pf$variable,OR$r_order)
pf$ID = factor(pf$ID,aa$Group.1)

ggplot(pf,aes(x=variable,y=value,fill=ID))+
  geom_bar(stat="identity",position="fill")+
  #scale_fill_manual(values = c('#8dd3c7','#ffff99','#fb8072','#fdb462',brewer.pal(12,'Paired'),brewer.pal(9,'Set3')))+
  #scale_color_manual(values = c(brewer.pal(12,'Paired'),brewer.pal(9,'Set3')))+
  #guides(color=F,fill=F)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



#SGB稀释曲线============================================================================

library(parallel)
Threads = 50
func = function(y){
  dat = sample(1:nrow(pf),size = y,replace = F)
  dat = data.frame(Step=y,idn=x,num=length(unique(pf$V1[dat])),
                   rb=sum(prof$rb[prof$ID%in%unique(pf$V1[dat])]),
                   stringsAsFactors = F)
  return(dat)
}

prof = read.csv('prokaryote_SGB.species.profile',sep='\t',stringsAsFactors = F)
p = read.csv('mag.cluster',sep='\t',stringsAsFactors = F,header = F)
p = tidyr::separate_rows(p[,c(1,4)],V4,sep=',')
colnames(p) = c('V1','V2')
#p$V1 = as.numeric(as.factor(p$V1))
prof$ID = sub(':.*','',prof$ID)

prof = data.frame(ID=prof$ID,rb=as.numeric(apply(prof[,-1], 1, mean)),stringsAsFactors = F)

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

f1=ggplot(aa)+
  geom_line(aes(x=Step,y=x,group=idn))+
  geom_errorbar(aes(x=Step,y=x,group=idn,ymin=x-sd,ymax=x+sd))+
  labs(x='Number of MAGs',y='Number of SGBs')+
  theme_bw()

aa = aggregate(Tab$rb,Tab[,1:2],mean)
aa$sd = aggregate(Tab$rb,Tab[,1:2],sd)$x
aa = aa[aa$idn==1,]

f2=ggplot(aa)+
  geom_line(aes(x=Step,y=x,group=idn))+
  geom_errorbar(aes(x=Step,y=x,group=idn,ymin=x-sd,ymax=x+sd))+
  labs(x='Number of MAGs',y='Accumulation of relative abundance')+
  theme_bw()

ggarrange(f1,f2,align = 'hv')

#top100 SGB 丰度、binning数量排序=========================================================================================
p = read.csv('prokaryote_SGB.species.profile',sep='\t',row.names = 1,stringsAsFactors = F)
mag = read.csv('mag.cluster',sep='\t',header = F,stringsAsFactors = F)
mag = data.frame(SGB=mag$V1,MAG=mag$V4,phylum=magClstr$Phylum[match(mag$V1,magClstr$SGB_ID)],
                 species=magClstr$Species[match(mag$V1,magClstr$SGB_ID)])


tab = mag[1:100,c('SGB','species','MAG','phylum')]
tab = tidyr::separate_rows(tab,MAG, sep = ",")
tab$sample = sub('.[ms]bin.*','',tab$MAG)
tab$Name = sub('_',' ',paste(tab$SGB,tab$species,sep=':'))
#tab$type = magClstr$Cultured[match(tab$SGB,magClstr$SGB_ID)]
tab$type = 'MAG'
tab$type[grep('^GC[AF]',tab$MAG)] = 'Isolate'
tab$avg_rb = as.numeric(rowMeans(p))[match(tab$SGB,sub(':.*','',row.names(p)))]
tab$box = magClstr$Cultured[match(tab$SGB,magClstr$SGB_ID)]
tab$Name = factor(tab$Name,rev(sort(unique(tab$Name))))

aa = plyr::count(tab[,c('Name','type')])
f1=ggplot(aa)+
  geom_point(aes(y=freq,x=Name,color=type))+
  guides(fill=F,color=F)+
  scale_y_continuous(trans = 'sqrt',breaks = c(50,500,1000,1500))+
  theme_classic()

library(ggalt)
aa = unique(tab[,c('Name','avg_rb')])
f2=ggplot(aa)+
  geom_point(aes(y=avg_rb,x=Name),shape=17)+
  geom_xspline(aes(y=avg_rb,x=as.numeric(Name)),spline_shape = 0.5)+
  geom_bar(aes(y=0.01,x=as.numeric(Name),fill=rep(c('a','b'),100)[1:length(unique(aa$Name))]),
           position = 'dodge',stat='identity',width = 1,alpha=0.05)+
  scale_y_continuous(trans = 'log10')+
  guides(fill=F)+
  theme_classic()

ggarrange(f1,f2,nrow=2,align = 'hv')



aa = unique(tab[,c('Name','box')])
ggplot(aa)+
  geom_bar(aes(y=1,x=Name,fill=box),
           position = 'dodge',stat='identity',width = 1)+
  theme_classic()

p = read.csv('VMGC_prokaryote_MAG.info',sep='\t',stringsAsFactors = F)
p$country = map$Country[match(p$BioSample_ID,map$Sample_alias)]
p = p[!is.na(p$country) & p$Species.level_genomic_bin_.95._ANI.%in%tab$SGB,]
p = plyr::count(p[,20:21])
p$Name = tab$Name[match(p[,1],tab$SGB)]
p$Name = factor(p$Name,levels(tab$Name))

ggplot(p,aes(x=Name,y=country,fill=freq))+
  geom_tile()+
  scale_fill_gradientn(trans='sqrt',colours = rev(brewer.pal(12,'Spectral')))+
  theme_test()+
  theme(
    axis.text.x = element_text(
      angle=90
    )
  )
