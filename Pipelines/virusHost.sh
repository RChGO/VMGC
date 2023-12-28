#!/usr/bin/bash

##################################################
# Ruochun Guo (grchun@hostmail.com)
# Date: Apr, 2022 (v1.0)
##################################################

if [ "$#" -ne "5" ]; then
  echo -e "\nusage: sh $0 <genome fasta> <output file> <CRISPRs db> <Genome db> <Taxlist>\n" 
    exit 2
fi


fa=$1
outf=$2
db_crass=$3
db_genome=$4
db_tax=$5

#db_crass='/share/data1/Database/uhgg/csp.tot.fa'
#db_genome='/share/data1/Database/uhgg/genomes_rep.fasta'
#db_crass='/share/data1/Database/uhgg_20220401/csp.tot.fa'
#db_genome='/share/data1/Database/uhgg_20220401/genomes_rep.fasta'
#db_tax='/share/data1/Database/uhgg_20220401/uhgg.taxo.group'

blastn -query $fa -db $db_crass -evalue 1e-2 -out $outf.csp.bt -outfmt 6 -num_alignments 999999 -num_threads 50 -word_size 8
filter_blast -i $outf.csp.bt -o $outf.csp.bt.f --evalue 1e-5 --score 45 --tops 20
less $outf.csp.bt.f | perl -e 'while(<>){chomp;@s=split /\s+/;($a,$b)=($s[0],$s[1]);$b=~s/_\d+.sp\d+$//; push @{$g{$a}},$b unless exists $h{"$a $b"}; $h{"$a $b"}++;} for(sort keys %g){@a=@{$g{$_}};print "$_\t".(join ",",@a)."\n";}' > $outf.csp.list


blastn -query $fa -db $db_genome -evalue 1e-2 -out $outf.fa.bt -outfmt 6 -num_alignments 999999 -num_threads 50
connect_blast $outf.fa.bt $outf.fa.bt.conn 1
/share/data2/guorc/script/fasta_length -i $fa |sort -rnk2 > $outf.sort.len
filter_blast -i $outf.fa.bt.conn -o $outf.fa.bt.conn.f --evalue 1e-10 --qfile $outf.sort.len --qper 30 --identity 90
less $outf.fa.bt.conn.f | perl -e 'while(<>){chomp;@s=split /\s+/;($a,$b)=($s[0],$s[1]);$b=~s/_\d+$//; push @{$g{$a}},$b unless exists $h{"$a $b"}; $h{"$a $b"}++;} for(sort keys %g){@a=@{$g{$_}};print "$_\t".(join ",",@a)."\n";}' > $outf.fa.list

cat $outf.*.list |sort |perl -ne 'chomp;@a = split /[,\t]/; foreach $x (@a){print "$a[0]\t$x\n"}' |awk '$1!=$2'|sort -u > $outf.host
/share/data2/guorc/script/TableTreat.py  -a $outf.host -b $db_tax -o $outf.host.tax -m table_merge -x 1 -j 0
