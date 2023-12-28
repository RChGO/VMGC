

#fungi进化树marker挑选===============================================================================================

p = read.csv('22.fungi/09.tree/all.clstr_cluster.tsv',header=F,sep='\t')
mp = read.csv('22.fungi/09.tree/all.faa.len',header=F,sep='\t')

p$sample = sub('_\\d+_g','',p$V2)
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

length(unique(p$V1[p$prev>0.8 & p$single_rate==1 & p$len_cv<0.2]))
#write.table(p,'22.fungi/09.tree/all.clstr.stat',sep='\t',row.names = F,quote = F)




#fungi分布===============================================================================================
cleanread = read.csv('13.bac_prof.list.cleanRead_gt_400k',sep='\t',header = F,stringsAsFactors = F)
cleanread = cleanread[cleanread$V2>1e7,]


Flist = sub('.log','',list.files(path = '22.fungi/05.map/',pattern = '*.log'))
bb = c()
for(x in Flist){
  aa = readLines(paste('22.fungi/05.map/',x,'.log',sep=''))
  aa = aa[grep('reads; of these',aa)][1]
  bb = c(bb,as.numeric(sub(' .*','',aa)))
}
tab = data.frame(ID=Flist,all=bb)


fungi_info = read.csv('22.fungi/fungi.info',sep='\t',stringsAsFactors = F)
p = read.csv('22.fungi/05.map.raw_prof',sep='\t',stringsAsFactors = F)
p[,tab$ID[!tab$ID%in%colnames(p)]] = 0
p = p[,c('ID',tab$ID)]
p = p[,colnames(p)%in%c('ID',cleanread$V1)]

dat = sweep(p[,-1],2,tab$all,'/')
row.names(dat) = paste(fungi_info$ID,fungi_info$species,sep=':')[match(p$ID,fungi_info$ID)]
dat = melt(as.matrix(dat))
dd = aggregate(dat$value,list(dat$Var1),mean)
dd = dd[order(dd$x),]
dd$Group.1 = factor(dd$Group.1,dd$Group.1)
dd$lab = paste(round(dd$x*1000,2),'‰',sep='')
dd$g = 'Abun'

dat = p[,-1]
row.names(dat) = paste(fungi_info$ID,fungi_info$species,sep=':')[match(p$ID,fungi_info$ID)]
dat = melt(as.matrix(dat))
aa = aggregate(dat$value,list(dat$Var1),function(x){sum(x>2)})
aa$x = aa$x/nrow(cleanread)
aa$lab = paste(round(aa$x*100,2),'%',sep='')
aa$g = 'Pre'

dd = rbind(dd,aa)

ggplot(dd)+
  geom_bar(aes(y=Group.1,x=x),position = 'dodge',stat = 'identity',width = .7)+
  geom_text(aes(y=Group.1,x=x,label=lab),size=3)+
  facet_grid(~g,scales = 'free')+
  #scale_x_continuous(trans = 'sqrt')+
  theme_bw()+
  theme(
    #axis.text.x = element_text(
    #  angle = 45,
    #  hjust = 0,
    #  vjust = 0.5
    #)
  )



#p[p==1] = 0
dat = sweep(p[,-1],1,fungi_info$size[match(p$ID,fungi_info$ID)],'/')
dat = sweep(dat,2,colSums(dat),'/')
row.names(dat) = paste(fungi_info$ID,fungi_info$species,sep=':')[match(p$ID,fungi_info$ID)]
rowMeans(dat,na.rm = T)

dat = data.frame(ID=pf$ID,val=apply(pf[,-1],1,function(x){sum(x>3)}))
dat$rate =dat$val/length(Flist)

#pf[,-1] = sweep(pf[,-1],2,colSums(pf[,-1]),'/') 

dat$b=rowMeans(pf[,-1])

dat$ID = factor(dat$ID,dat$ID[order(-dat$rate)])
ggplot(dat)+
  geom_bar(aes(x=ID,y=rate),position = 'dodge',stat = 'identity',width = .7)+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )
