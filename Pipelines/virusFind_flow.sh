#!/usr/bin/bash


##################################################
# Ruochun Guo (grchun@hostmail.com)
# Date: May, 2022 (v2.0)
##################################################


if [ "$#" -ne "4" ]; then
  echo -e "\nusage: sh $0 <fasta> <minSize for virus> <outfile_prefix> <outdir>" 
  echo -e "version 2.0\n" 
    exit 2
fi

fa=$1
minLen=$2
id=$3
outdir=$4

#minLen=3000
Threads=20

mkdir -p $outdir

###### steq1
/share/data2/guorc/script/seqkit seq -m $minLen $fa | perl -ne 'if(/>/){$a++;print ">'${id}'_ctg$a\n"}else{print $_}' > $outdir/$id.fa
cd $outdir

###### step2 (found by checkv)
/share/data2/guorc/Software/conda/checkv/bin/checkv end_to_end -d /share/data2/guorc/Software/conda/checkv/checkv-db-v0.6/ -t $Threads $id.fa $id.ckv
awk '($3=="No"||$1=="contig_id") && !($7/($5+1e-16)>0.5 && ($6==0 || ($6>0 && $7/($6+1e-16)>5)))' $id.ckv/quality_summary.tsv > $id.vf
mv -f $id.vf $id.ckv/quality_summary.tsv.tmp
cut -f1 $id.ckv/quality_summary.tsv.tmp | seqkit grep -f - $id.fa > $id.fa2
mv -f $id.fa2 $id.fa
#perl -ne 'chomp;@b = split / /; $a=join(".",split(/_([^_]+)$/, $b[0]));print "$a\n"' $id.ckv/proviruses.fna >> $id.fa
sed 's/ .*//' $id.ckv/proviruses.fna >> $id.fa
/share/data2/guorc/Software/conda/checkv/bin/checkv end_to_end -d /share/data2/guorc/Software/conda/checkv/checkv-db-v0.6/ -t $Threads $id.ckv/proviruses.fna $id.ckv.tmp
sed '1d' $id.ckv.tmp/quality_summary.tsv |awk '$2>'${minLen}'' >> $id.ckv/quality_summary.tsv.tmp
cat $id.ckv.tmp/tmp/proteins.faa >> $id.ckv/tmp/proteins.faa
rm -rf $id.vf $id.ckv.tmp

###### step3 (found by vibrant)
source activate /share/data1/software/miniconda3/envs/vibrant
VIBRANT_run.py -i $id.fa -virome -l $minLen -f nucl -folder $id.vib -t $Threads -no_plot

###### step4 (found by virFinder)
source activate /share/data1/software/miniconda3/envs/dvf
python /share/data1/lvqb/software/DeepVirFinder/dvf.py -i $id.fa -o $id.vif -l $minLen -c $Threads



###### step5 (filter)
awk '$6>$7 && $8!="Not-determined"' $id.ckv/quality_summary.tsv.tmp  > $id.ckv/quality_summary.tsv.f 
perl -ne 's/_frag.*$//;print $_;' $id.vib/VIBRANT_$id/VIBRANT_phages_$id/$id.phages_combined.txt > $id.vib/$id.vib.f
awk '{if($4<0.01)print $1}' $id.vif/$id.fa_gt${minLen}bp_dvfpred.txt > $id.vif/$id.fa_gt${minLen}bp_dvfpred.txt.f

cat $id.ckv/quality_summary.tsv.f $id.vib/$id.vib.f $id.vif/$id.fa_gt${minLen}bp_dvfpred.txt.f | sort -u > $id.virus
/share/data2/guorc/script/diff_head $id.virus 1 $id.ckv/quality_summary.tsv.tmp 1 $id.virus.tsv
rm $id.ckv/quality_summary.tsv.f $id.vib/$id.vib.f $id.vif/$id.fa_gt${minLen}bp_dvfpred.txt.f $id.virus 



###### step6 (decontaminated by BUSCO)
cut -f1 $id.virus.tsv|sed '1d'|sed 's/$/_/' |seqkit grep -f - -r $id.ckv/tmp/proteins.faa > $id.virus.faa

hmmsearch --tblout $id.busco --cpu $Threads --noali /share/data1/lish/virome/2020_kideny/27.scaffolding/busco/bacteria_odb10/tot.hmm $id.virus.faa > $id.busco.log 
perl /share/data1/lish/virome/2020_kideny/27.scaffolding/busco/filter_hmm2.pl $id.busco /share/data1/lish/virome/2020_kideny/27.scaffolding/busco/bacteria_odb10/scores_cutoff $id.busco.f $id.virus.faa $id.busco.f2
awk '$3/$2<0.05' $id.busco.f2 > $id.busco.final

/share/data2/guorc/script/diff_head $id.busco.final 1 $id.ckv/quality_summary.tsv.tmp 1 $id.virus.checkv
#perl /share/data2/guorc/script/get_fasta.pl $id.virus.checkv $id.fa $id.virus.fa
sed '1d' $id.virus.checkv|cut -f1 | seqkit grep -f - $id.fa > $id.virus.fa
rm $id.fa $id.virus.faa $id.virus.tsv $id.busco.final $id.busco.f $id.busco.f2 $id.ckv/quality_summary.tsv.tmp


