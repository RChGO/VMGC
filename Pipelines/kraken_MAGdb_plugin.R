#!/usr/bin/env Rscript

##################################################
# Ruochun Guo (grchun@hostmail.com)
# Date: Oct., 2023 (version 1.0)
##################################################


argv<-commandArgs(trailingOnly = F)
sp <- sub('--file=','',argv[4])
args <- argv[-(1:5)]

if(length(args)<3){
  cat('\nUsage:\nRscript ',sp,' [db.tax] [db.seqid] [output dir] \n\n')
  quit(save = 'no')
}

intax=args[1]
inseq=args[2]
outd=args[3]

Level = c('superkingdom','phylum','class','order','family','genus','species')

#p = read.csv(intax,sep='\t',stringsAsFactors = F)
#p = data.frame(row.names = p[,1],p[,4:10],stringsAsFactors = F)

p = read.csv(intax,sep='\t',stringsAsFactors = F,row.names = 1)
db = read.csv(inseq,sep='\t',stringsAsFactors = F,header = F)

Matx = matrix(nrow = nrow(p),ncol = ncol(p),dimnames = list(row.names(p),Level))
Variable = c()
NodeNum = 1
for(x in 1:ncol(p)){
  Variable = paste(Variable,p[,x],sep=';')
  Num = as.numeric(as.factor(Variable)) + NodeNum
  Matx[,x] = Num
  NodeNum = max(Matx[,x])
}

Matx = data.frame('no rank' = '1',Matx,check.names = F)
Node = data.frame(Child='1',Parent='1',Level='no rank')
Name = data.frame(Node='1',Name='no rank',Anno='',Level='scientific name')
for(x in 2:ncol(Matx)){
  
  dat = data.frame(unique(Matx[,x:(x-1)]),Level=colnames(Matx)[x])
  colnames(dat)[1:2] = c('Child','Parent')
  dat$Child = as.character(dat$Child)
  dat$Parent = as.character(dat$Parent)
  Node = rbind(Node,dat)
  
  if(x!=ncol(Matx)){
    Anno=''
  }else{
    Anno=row.names(Matx)
  }
  dat = unique(data.frame(Node=as.character(Matx[,x]),Name=p[,x-1],
                          Anno=Anno,Level='scientific name'))
  Name = rbind(Name,dat)
}

Node = paste(apply(Node,1,function(x){paste(x,collapse = '\t|\t')}),'\t|',sep='')
Name = paste(apply(Name,1,function(x){paste(x,collapse = '\t|\t')}),'\t|',sep='')

write.table(Node,paste(outd,'nodes.dmp',sep='/'),row.names = F,col.names = F,quote = F)
write.table(Name,paste(outd,'names.dmp',sep='/'),row.names = F,col.names = F,quote = F)

db = data.frame(accession = sub('\\.\\d+$','',db[,2]), 
            accession.version = paste(db[,2],sep=''),
            taxid = Matx[match(db[,1],row.names(Matx)),ncol(Matx)],
            gi=1:nrow(db),
            check.names = F
  )
write.table(db,paste(outd,'db.accession2taxid',sep='/'),row.names = F,quote = F,sep='\t')


