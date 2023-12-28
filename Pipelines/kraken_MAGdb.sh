#!/usr/bin/bash

##################################################
# Ruochun Guo (grchun@hostmail.com)
# Date: Oct., 2023 (version 1.0)
##################################################



if [ "$#" -le "3" ]; then
  echo -e "\nusage: sh $0 <genome TaxInfo> <genome fasta list> <outdir> <Threads (default 10)> <brackendb readlen (default 150)>\n" 
  echo -e "\ngenome TaxInfo format (The colnames for taxonomy must exist in the NCBI's nodes.dmp file):"
  echo -e "ID      superkingdom    phylum          class           order             family             genus           species
           SGB001  Bacteria        Bacillota       Bacilli         Lactobacillales   Lactobacillaceae   Lactobacillus   Lactobacillus_iners
           SGB002  Bacteria        Actinomycetota  Coriobacteriia  Coriobacteriales  Atopobiaceae       Fannyhessea     Fannyhessea_vaginae" | column -t
  echo -e "\ngenome fasta list format:"
  echo -e "SGB001  00.data/sgb_fasta/SGB001.fa
           SGB002  00.data/sgb_fasta/SGB002.fa" | column -t
  exit 2

fi

intax=$1
inseqlist=$2
outdir=$3
Threads=${4:-10}
brackendb_readlen=${5:-150}

#intax='sgb.info'
#inseqlist='sgb.seq'
#outdir='ttest'
#Threads=10

mkdir -p $outdir/taxonomy

/share/data2/guorc/script/TableTreat.py -a $intax -b $inseqlist -x 0 -m table_merge -o $outdir/Temp.list
perl -i -ne 'if($. ne 1){$a=$.-1; s/(.*?)\t/KRA_$a\t/};print $_' $outdir/Temp.list 
awk '{if($NF!="NA"){print $1"\t"$NF}}'  $outdir/Temp.list | awk '{print "seqkit replace -p \"(.*)\" -r \""$1"CTG{nr}\" "$2}' | sh >  $outdir/db.fna
seqkit seq -n $outdir/db.fna | perl -ne '/(KRA_\d+)/;print "$1\t$_"' > $outdir/db.seqid
awk -F '\t' 'BEGIN{OFS=FS} {$NF=""; sub(/\t$/,"")} 1' $outdir/Temp.list > $outdir/db.tax

/share/data2/guorc/script/kraken_MAGdb_plugin.R $outdir/db.tax $outdir/db.seqid $outdir/taxonomy

kraken2-build --add-to-library $outdir/db.fna --db $outdir
kraken2-build --build --db $outdir  --threads $Threads
bracken-build -d $outdir -t $Threads -k 35 -l $brackendb_readlen

rm -rf $outdir/Temp.list $outdir/db.fna $outdir/db.seqid $outdir/db.tax


echo -e '\nFinish!!: kraken database for MAGs\n'
