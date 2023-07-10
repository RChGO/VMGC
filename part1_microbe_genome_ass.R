
#样品信息==================================================================================

p = read.csv('../../data/mapping.file/mapping.file',sep='\t')
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
  geom_polygon(aes(fill= freq), colour = "black",size=0.1) +
  scale_x_continuous(breaks = seq(-180, 210, 45), labels = function(x){paste0(x, "°")}) +
  scale_y_continuous(breaks = seq(-60, 100, 30), labels = function(x){paste0(x, "°")}) +
  scale_fill_gradientn(tran='log10',colours = c('white','#E6F5D0',brewer.pal(11,'RdYlBu')[c(9:11)]))+
  #scale_fill_gradient(tran='log10',low = "#fee0d2", high="#ef3b2c")+
  theme_light()



#原核基因组信息======================================================================================================
mag_info  = read.csv('00.data/mag.info',sep='\t',stringsAsFactors = F)

#完整性情况
mag_info$qg = factor(mag_info$qg,c("near-complete","high-quality","medium-quality"))
f3=ggplot(mag_info)+
  geom_point(aes(x=Completeness,y=Contamination,color=qg),alpha=.3,size=.5)+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_test()

f1=ggplot(mag_info)+
  geom_histogram(aes(x=Completeness,y=..count..,group=qg,fill=qg),color='black',bins = 50)+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_classic()


f4=ggplot(mag_info)+
  geom_histogram(aes(y=Contamination,x=..count..,group=qg,fill=qg),color='black',bins = 50)+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_classic()



aa = plyr::count(mag_info$qg)
aa$rate = aa$freq/sum(aa$freq)
aa$x = factor(aa$x,c("near-complete","high-quality","medium-quality"))
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


ggarrange(f1,f2,f3,f4,align = 'hv')


##基因组组装情况
f1=ggplot(mag_info)+
  geom_boxplot(aes(x=qg,y=contig_num,fill=qg),outlier.shape = 1,outlier.size = .5)+
  #scale_y_sqrt(breaks=c(0,200,400,600))+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_classic()

pf = mag_info[,c('qg','totalLength')]
pf = rbind(pf,data.frame(qg='all',totalLength=p$totalLength))
f2=ggplot(pf)+
  geom_boxplot(aes(x=qg,y=totalLength,fill=qg),outlier.shape = 1,outlier.size = .5)+
  scale_y_log10()+
  scale_fill_manual(values = c(Col,'grey'))+
  scale_color_manual(values = c(Col,'grey'))+
  guides(color=F,fill=F)+
  theme_classic()

pf = mag_info[,c('qg','N50')]
pf = rbind(pf,data.frame(qg='all',N50=p$N50))
f3=ggplot(pf)+
  geom_boxplot(aes(x=qg,y=N50,fill=qg),outlier.shape = 1,outlier.size = .5)+
  scale_y_log10()+
  scale_fill_manual(values = c(Col,'grey'))+
  scale_color_manual(values = c(Col,'grey'))+
  guides(color=F,fill=F)+
  theme_classic()

f4=ggplot(mag_info)+
  geom_boxplot(aes(x=qg,y=trna,fill=qg),outlier.shape = 1,outlier.size = .5)+
  scale_y_sqrt()+
  scale_fill_manual(values = Col)+
  scale_color_manual(values = Col)+
  guides(color=F,fill=F)+
  theme_classic()

ggarrange(f2,f3,f2,f3,align = 'hv')


#真核基因组信息====================================================================================

fungi_info  = read.csv('22.fungi/fungi.info',sep='\t',stringsAsFactors = F)

ggplot(fungi_info)+
  geom_point(aes(x=completeness.,y=contamination.),color='#3E90C9',alpha=.8,size=.5)+
  #scale_fill_manual(values = '3E90C9')+
  #scale_color_manual(values = '3E90C9')+
  guides(color=F,fill=F)+
  theme_test()



#病毒基因组信息====================================================================================
vs_info  = read.csv('00.data/final.virus.ckv',sep='\t',stringsAsFactors = F,header=F)

aa = plyr::count(vs_info$V8)
aa$rate = aa$freq/sum(aa$freq)
aa$x = factor(aa$x,c("Complete","High-quality","Medium-quality"))
aa = aa[order(aa$x),] 
ggplot(aa,aes(x=1,y=rate,fill=x))+
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
ggplot(pf)+
  geom_boxplot(aes(x=V8,y=V2,fill=V8),outlier.shape = 1,outlier.size = .5)+
  scale_y_log10()+
  scale_fill_manual(values = c(Col,'grey'))+
  scale_color_manual(values = c(Col,'grey'))+
  guides(color=F,fill=F)+
  theme_classic()
