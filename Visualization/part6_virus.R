

#病毒物种==================================================================================================================

ckv = read.csv('00.data/votu.ckv',sep='\t',header = F)
ckv = ckv[ckv$V1!='MG1145.k121_3894',] #illu测序常用阳性对照，怀疑样品测序污染引入
ckv$id = paste('v',as.hexmode(1:nrow(ckv)),sep='')
#write.table(ckv,'aa',row.names = F,sep='\t',quote = F)
summary(ckv$V2)
sum(ckv$V2>200000)

vtax = read.csv('15.votu_tax/votu.tax_family',sep='\t',header = F)
ckv$tax = vtax$V4[match(ckv$V1,vtax$V1)]
ckv$tax[is.na(ckv$tax)] = 'unk'

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

votu_clstr = read.csv('00.data/final.virus.fa.bts.cvg.uniq.list', sep='\t', header = F, stringsAsFactors = F)
prof = read.csv('19.votu_host/votu.host.tax',sep='\t',header = F)
hostf = multip_vir_host(prof)



ckv$host = hostf$V3[match(ckv$V1,hostf$V1)]
ckv$host[is.na(ckv$host)] = 'unk'
ckv$host[ckv$tax%in%c('Papillomaviridae','Virgaviridae','Circoviridae','Metaviridae','Iflaviridae','Totiviridae')] = 'Eukaryote'
ckv$num = votu_clstr$V2[match(ckv$V1,votu_clstr$V1)]

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
  theme_classic()


#病毒宿主对应关系
aa = tidyr::separate_rows(magClstr[,c(1,4)],'MAG',sep=',')
dat = prof
dat$sgb = aa$SGB[match(dat$V2,aa$MAG)]
dat = unique(dat[,c('V1','V3','V8','sgb')])

dat = aggregate(dat[,-1],list(dat$V1),function(x){paste(sort(unique(x)),collapse = ', ')})
#write.table(dat,'aa',sep = '\t',quote = F,row.names = F)


#HPV菌分布===================================================================================================================

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
