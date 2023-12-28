## VMGC - Human Vaginal Microbiome Genome Collection
 The VMGC is a large-scale reference genome resource of the human vaginal microbiome, including over 33,000 genomes derived from 786 prokaryotes, 38 fungi, and 4,263 viruses associated with the human vagina. In terms of representation, the VMGC demonstrates high efficiency in capturing microbial sequences, with a median mapping rate of 91.7% across 4,472 vaginal metagenomic samples obtained from 14 countries.</font>
<br>
<br>
The microbial sequence data for VMGC is stored at http://puensum.tpddns.cn:8000/pub/VMGC/, including the data presented in the table below.
<br>
|Description|Size|Filename|
|  ----  | ---- | ---- |
|  Full prokaryotic genomes (n=19,541)  | 8.1 GB | VMGC_prokaryote_MAG.tar.gz |
| Annotations of prokaryotic genomes (Source/Quality/Clustering...)  | 2.3 MB | VMGC_prokaryote_MAG.info |
| Nonredundant prokaryotic genomes (n=786) | 440 MB | VMGC_prokaryote_SGB.tar.gz |
| Annotations of nonredundant prokaryotic genomes (Taxonomy/...) | 120 KB | VMGC_prokaryote_SGB.info |
| Full eukaryotic genomes (n=42, including 4 from parasites)  | 192 MB | VMGC_eukaryocyte.tar.gz |
| Annotations of eukaryotic genomes (Taxonomy/Quality/Clusting/...)  | 11 KB | VMGC_eukaryocyte.info |
| Full viral genomes (n=14,224) and vOTUs (n=4,263)| 197 MB | VMGC_virus.tar.gz |
| Annotations of eukaryotic genomes (Taxonomy/Quality/Clusting/Host...)  | 969 KB | VMGC_virus.info |
|  Kraken & Bracken database  | 2.1 GB | VMGC_prokaryote_SGB_KrakenDB.tar.gz |

<br>

---
### Single-sample and Mash-based multiple-sample binning
Metagenome-assembled genomes in the VMGC were obtained through a process that integrated both single-sample binning and Mash-based multiple-sample binning methods. The specific operational steps are outlined below.

* Required dependencies
    > Mash 2.3<br>
    > GNU parallel 20201122<br>
    > bwa-mem2 2.2.1<br>
    > MetaBAT v2<br>
    > perl 5.16.3<br>
    > [combine.pl](https://github.com/WatsonLab/single_and_multiple_binning/blob/main/scripts/combine.pl)

* Step 1: Calculate the Mash distances between the assembled files, and obtain the top 20 closest samples for each assembled sample.
    >mash sketch  -p 50 -s 100000 -k 32  -o mash.sketch contigs/*.fasta

    >find contigs/*.fasta | parallel -k -j 10 mash dist mash.sketch.msh {} \\| sort -nk3 \\| head -n 20 \\| sed -e "s/contigs\\\\\\///g" -e "s/.fasta//g" > mash.sketch.msh.top20

* Step 2: For a assembly file, clean reads from the top 20 closest samples (including its own reads) by Mash distance were used to calculate the sequencing depth of contigs.
    >find contigs/*.fasta | parallel --colsep '\t' -j 20 bwa-mem2 index -p {} {}

    >cat mash.sketch.msh.top20 | parallel --colsep '\t' -j 10 mkdir -p depth/{2} \\; bwa-mem2 mem -t 10 contigs/{2} clean_reads/{1}.1.fq.gz clean_reads/{1}.2.fq.gz \\| samtools view -bS - -@ 10 \\| samtools sort -@ 10 -o depth/{2}/{1}.sort.bam \\&\\& jgi_summarize_bam_contig_depths --outputDepth depth/{2}/{1}.sort.bam.depth depth/{2}/{1}.sort.bam

* Step 3: Integrate depth calculation files using the public script combine.pl, and perform  multiple-sample binning.
    >find depth/* -type d | parallel -j 10 combine.pl {}/*.sort.bam.depth \\> {}.depth

    >find depth/*.depth | parallel -j 10 mkdir -p bins/{/.}/ \\; metabat2 -i contigs/{/.}.fasta -a {} -o bins/{/.}/{/.}.mbin -m 2000 -s 200000 --saveCls --unbinned --seed 2020 

* Step 4: Perform single-sample binning.
    >find depth/* -type d | parallel -j 10 metabat2 -i contigs/{/.}.fasta -a depth/{/.}/{/.}.depth -o bins/{/.}/{/.}.sbin -m 2000 -s 200000 --saveCls --unbinned --seed 2020
<br>

---
### Taxonomic profiling
Based on the 786 species-level genome bins (SGBs) in the VMGC, we reconstructed the prokaryotic composition of the vagina using Kraken2 and Bracken tools.

* Required dependencies
    > Kraken 2.1.3<br>
    > Bracken 2.8<br>
    > Python 3.8.16<br>
    > R 4.2.3<br>
    > GNU parallel 20201122<br>
    > [kraken_MAGdb.sh](https://github.com/RChGO/VMGC/Pipelines/kraken_MAGdb.sh)<br>
    > [kraken_MAGdb_plugin.R](https://github.com/RChGO/VMGC/Pipelines/kraken_MAGdb_plugin.R)<br>

* Step 1: Create customized Kraken2 and Bracken databases.
    > kraken_MAGdb.sh [sgb.info](https://github.com/RChGO/VMGC/Documents/sgb.info) [sgb.seq](https://github.com/RChGO/VMGC/Documents/sgb.seq) KBdb<br>

        1. The information on "sgb.info" is derived from the annotation results of GTDB-tk.
        2. "sgb.seq" includes the file paths for the genomes of SGBs.
        3. "KBdb" is the output folder path, and the generated data includes the following files:
            >tree KBdb/
            KBdb/
            ├── database150mers.kmer_distrib
            ├── database150mers.kraken
            ├── database.kraken
            ├── hash.k2d
            ├── library
            │   └── added
            │       ├── nFPH00iyWZ.fna
            │       ├── nFPH00iyWZ.fna.masked
            │       ├── prelim_map.txt
            │       └── prelim_map_ZLPtIxoxoZ.txt
            ├── opts.k2d
            ├── seqid2taxid.map
            ├── taxo.k2d
            └── taxonomy
                ├── db.accession2taxid
                ├── names.dmp
                ├── nodes.dmp
                └── prelim_map.txt

* Step 2: The clean reads from each sample are mapped to the database, generating compositions at various taxonomic levels.
    > find clean_reads/*.1.fq.gz | sed 's/.1.fq.gz//'| parallel -j 5 kraken2 --threads 10 --db KBdb --report prof/{/}.report --report-minimizer-data --output prof/{/}.output {}.1.fq.gz {}.2.fq.gz

    > find prof/*.report | parallel -j 5 bracken -d KBdb -i {} -o {}.bracken -r 150 -l S -t 1

<br>

---
### Fungal genome identification
We utilized the aforementioned process, which integrates both single-sample binning and Mash-based multiple-sample binning methods, to generate raw bins. Subsequently, we employed EukRep and BUSCO tools to extract fungal genomes from raw bins.

* Required dependencies
    > EukRep 0.6.7<br>
    > BUSCO 5.4.2<br>
    > perl 5.16.3<br>
    > GNU parallel 20201122<br>

* Step 1: Remove the prokaryotic sequences from the raw bins.
    > find bins/\*/\*.[ms]bin.[0-9]\*.fa -size +3M | parallel -j 20 EukRep --min 2000 -i {} -o euk_bin/{/.}.fa

* Step 2: Assess the genome quality of candidate bins, and extract fungal genomes.
    > find euk_bin/*.fa -size +3M | parallel -j 20 busco -m genome -l fungi_odb10 -i {} -q -o busco/{/.}

    >grep "  C:" busco/*/*txt | perl -ne '/fungi_odb10.(.\*).txt.\*C:(.\*?)%.\*D:(.\*?)%/;print "ln -s euk_bin/$1.fa fungal_bins/\n" if $2>50 and $3<5' |sh

<br> 

---
### Viral genome identification

* Required dependencies
    > seqkit 2.1.0<br>
    > checkv 0.7.0<br>
    > VIBRANT 1.2.1<br>
    > DeepVirFinder 1.0<br>
    > BLAST 2.12.0+<br>
    > HMMER 3.3.2<br>
    > MinCED 0.4.2<br>
    > perl 5.16.3<br>
    > Python 3.8.16<br>
    > GNU Awk 4.0.2<br>
    > GNU parallel 20201122<br>
    > [virusFind_flow.sh](https://github.com/RChGO/VMGC/Pipelines/virusFind_flow.sh)<br>
    > [virusClstr.sh](https://github.com/RChGO/VMGC/Pipelines/virusClstr.sh)<br>
    > [vs_taxAnno.sh](https://github.com/RChGO/VMGC/Pipelines/vs_taxAnno.sh)<br>
    > [vs_host.sh](https://github.com/RChGO/VMGC/Pipelines/vs_host.sh)<br>
    > [filter_hmm2.pl](https://github.com/RChGO/VMGC/Pipelines/filter_hmm2.pl)<br>
    > [diff_head](https://github.com/RChGO/VMGC/Pipelines/diff_head)<br>
    > [blast_cvg.v1.py](https://github.com/RChGO/VMGC/Pipelines/blast_cvg.v1.py)<br>
    > [blast_cluster.v2.pl](https://github.com/RChGO/VMGC/Pipelines/blast_cluster.v2.pl)<br>
    > [fasta_length](https://github.com/RChGO/VMGC/Pipelines/fasta_length)<br>
    > [connect_blast](https://github.com/RChGO/VMGC/Pipelines/connect_blast)<br>
    > [TableTreat.py](https://github.com/RChGO/VMGC/Pipelines/TableTreat.py)<br>

* Step 1: Extract the viral sequences from the assembly files.
    > find contigs/*.fasta | parallel -j 20 virusFind_flow.sh {} 5000 {/.} virus/{/.}
* Step 2: Cluster viral sequences to generate Viral Operational Taxonomic Units (vOTUs).
    > virusClstr.sh virus/*.virus.fa votu clstr 50
* Step 3: Taxonomic annotation for vOTUs.
    > 
* Step 4: Virus-Host prediction based on the 19,541 prokaryotic MAGs in the VMGC.
    * Step 4.1: Find CRISPRs in all MAGs, and build BLAST databases.
        > find mag_fasta/\*.fa  | parallel -j 20 minced -minNR 2 {} csp/{/.}.csp csp/{/.}.gff
        
        > cat csp/\*.csp  |grep '^Sequence\\|\\[' | perl -ne 'if(/^Sequence.\*/){/'\''(\S+)'\''/; $a=$1;$b=1}else{/(\d+)(\s+)(\S+)(\s+)(\S+)/;print ">$a.sp$b\n$5\n";$b++}' > csp.tot.fa

        > makeblastdb -in csp.tot.fa -dbtype nucl  -out csp.tot.fa

        > cat mag_fasta/\*.fa > mag.fa && makeblastdb -in mag.fa -dbtype nucl -out mag.fa
    * Step 4.2:Predict hosts.
        > vs_host.sh clstr/votu.fa csp.tot.fa mag.fa [mag.tax](https://github.com/RChGO/VMGC/Documents/mag.tax)



<br>  

---

<br> 

##### See our paper for details:
<font size=2> Li S, Guo R, Zhang Y, et al. A catalogue of 48,425 nonredundant viruses from oral metagenomes expands the horizon of the human oral virome[J]. iScience, 2022: 104418. </font>

<font size=2> Correspondence and requests for materials should be addressed to grchun@hotmail.com </font>


