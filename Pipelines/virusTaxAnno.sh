#!/bin/bash

##################################################
# Ruochun Guo (grchun@hostmail.com)
# Date: May, 2022 (v4.0)
##################################################


if [ "$#" -ne "3" ]; then
    echo -e "\nusage: sh $0 <AA faa> <Genome fasta> <output file>\n" 
    exit 2
fi


prot=$1
geno=$2
outf=$3
threads=30

db=/share/data2/guorc/Database/virus/20220331/tot.faa
db_tax=/share/data2/guorc/Database/virus/20220331/tot.tax
db_s=/share/data2/guorc/Database/virus/20220331/db_genome.fa
db_s_tax=/share/data2/guorc/Database/virus/20220331/db_genome.tax



#基于蛋白序列的科水平注释
diamond blastp --threads $threads --max-target-seqs 10 --db $db --query $prot --outfmt 6 --out $outf.pep.bt --quiet 
perl /share/data2/guorc/script/get_pep.length.pl $prot $outf.pep.len
filter_blast -i $outf.pep.bt -o $outf.pep.bt2 --qfile $outf.pep.len --qper 50 --identity 30 --tops 40 --score 50 
grep '^>' $prot | sed 's/^>//' |perl -pne 's/_(\d+) #.*//' |sort |uniq -c |awk '{print $2"\t"$1}' > $outf.ngenes
perl /share/data2/guorc/script/WGS/Virus/vctg_stat.v2.pl $outf.pep.bt2 $db_tax $outf.ngenes $outf.pep.bt2.f
awk '$2!="NA"' $outf.pep.bt2.f > $outf.pep.bt2.f2
msort -k 1,rn3,rn6 $outf.pep.bt2.f2 | perl -ne 'chomp;@s=split /\s+/;next if exists $h{$s[0]}; $pct=$s[2]/$s[4]; next unless ($pct>=0.2 or $s[2]>=10) and $s[5]>=30; printf "$s[0]\t$s[2]/$s[4]\t%.2f\t$s[1]\n",$s[5];$h{$s[0]}=1;' > $outf.tax_family


#基于基因组序列的种水平注释
blastn -db $db_s -query $geno -out $outf.nucl.bt -outfmt 6 -evalue 1e-5 -num_threads $threads -max_target_seqs 50 -word_size 28
/share/data2/guorc/script/fasta_length -i $geno -o $outf.nucl.len
/share/data1/zhangy2/scripts/blast_cvg.v1.py $outf.nucl.bt $outf.nucl.cvg 0 200
/share/data2/guorc/script/TableTreat.v4.py -a $outf.nucl.cvg -b $outf.nucl.len -m table_merge -o $outf.nucl.cvg
cat $outf.nucl.cvg | perl -ne '@a=split/\t/; $cov=$a[2]/$a[6]*100; if($cov > 75 && $a[3] > 95){$ord=$a[3]*$cov; print "$a[0]\t$a[1]\t$a[3]\t$cov\t$ord\n"}' | msort -t '\t' -k 1,rn5 | perl -ne '@s=split /\s+/;print "$s[0]\t$s[1]\t$s[2]\t$s[3]\n" unless $a eq $s[0];$a=$s[0];' > $outf.nucl.cvg2
/share/data2/guorc/script/TableTreat.v4.py -a $outf.nucl.cvg2 -b $db_s_tax -o $outf.nucl.cvg2 -x 1 -y 1  -m table_merge
cat $outf.nucl.cvg2 | perl -ne '@a = split(/\t/,$_); @b=/;f_(.*?);/ ; @c=/;s_(.*?)$/; print "$a[0]\t$b[0]\t$c[0]\n"' | perl -pne 's/\t\t/\tUnclassified\t/' | sed -e '1i ID\tfamily\tspecies' -e  's/ /_/g' > $outf.tax_species



