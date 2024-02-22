
The files in this folder contain custom scripts required for constructing VMGC. Example files related to their usage are stored at https://github.com/RChGO/VMGC/blob/main/Documents/. Refer to the following instructions for more details.</font>
<br>

---

#### TableTreat.py<br>
* The Python script was designed to match and merge two tables.</font>
    > \>TableTreat.py -h<br>
    ```
    optional arguments:
    -h, --help            show this help message and exit
    -a A                  input table1, required
    -b B                  input table2, required
    -o O                  input name of outfile, required
    -m {table_column,table_row,table_merge}
                            pick out the needed rows or columns according to table2, default table_row
    -r                    whether to output variables that aren't in table2, default closed
    -w W                  merge type [left/right/inner/outer], default left
    -x X                  column for table1, default 0
    -y Y                  column for table2, default 0
    -i I                  whether or not table1 has header, default None, e.g. 0,1,2 
    -j J                  whether or not table2 has header, default None, e.g. 0,1,2
    -f F                  replace this of vacant position, default NA
    ```
    > \>TableTreat.py -a [sgb.info](https://github.com/RChGO/VMGC/blob/main/Documents/sgb.info) -b [sgb.seq](https://github.com/RChGO/VMGC/blob/main/Documents/sgb.seq) -x 0 -m table_merge -o Temp.list <br>


#### blast_cvg.v1.py<br>
* In long sequence blastn alignments, it's common for the alignments between two sequences to be fragmented. This Python script was developed to merge these fragmented alignment results, providing an overall assessment of similarity and alignment length between the two sequences.</font>
    > \>blast_cvg.v1.py -h
    ```
    usage: blast_cvg.v1.py [-h] in_f out_f min_id min_len
    
    positional arguments:
    in_f        blast.m8
    out_f       output
    min_id      ignore the identity less than min_id[default 75]
    min_len     ignore the identity less than min_len[default 200]

    optional arguments:
    -h, --help  show this help message and exit
    ```
    > \>blast_cvg.v1.py [blast_cvg.input](https://github.com/RChGO/VMGC/blob/main/Documents/blast_cvg.input) [blast_cvg.output](https://github.com/RChGO/VMGC/blob/main/Documents/blast_cvg.output) 0 0 <br>

#### blast_cluster.v2.pl<br>
* This Perl script was developed to cluster sequences based on the output of blast_cvg.v1.py.
    > \>blast_cluster.v2.pl -h
    ```
    perl blast_cluster.v2.pl [in.len.sort] [in.blast.cvg] [out.fasta] [%cvg] [%iden]
    ```
    > \>blast_cluster.v2.pl [blast_cluster.len.sort](https://github.com/RChGO/VMGC/blob/main/Documents/blast_cluster.len.sort) [blast_cvg.output](https://github.com/RChGO/VMGC/blob/main/Documents/blast_cvg.output) [blast_cluster.output](https://github.com/RChGO/VMGC/blob/main/Documents/blast_cluster.output) 85 95<br>

#### connect_blast<br>
* This compiled C program performs a similar function to blast_cvg.v1.py, with the difference being that the output file retains the blastn m8 format.
    > \>connect_blast -h
    ```
    Usage:   connect_blast [*.blast] [*.blast.conn] [sorted (0:didnot)]
    Version: 2.0 (2011-6-24)...
    ```
    > \>connect_blast [blast_cvg.input](https://github.com/RChGO/VMGC/blob/main/Documents/blast_cvg.input) [connect_blast.output](https://github.com/RChGO/VMGC/blob/main/Documents/connect_blast.output) 1<br>

#### diff_head<br>
* Based on the information in the previous table, this compiled C program extracts or remove the specified rows from the subsequent table, which must have headers.
    > \>diff_head -h 
    ```
    Usage: diff_head [a.list] [column a] [b.list] [column b] [share.list] [1 or -1]
    ```
    > \>diff_head [sgb.seq](https://github.com/RChGO/VMGC/blob/main/Documents/sgb.seq) 1 [sgb.info](https://github.com/RChGO/VMGC/blob/main/Documents/sgb.info) 1 Temp.list -1 

#### filter_blast<br>
* This compiled C program is used to filter the results of blastn alignments.
    > \>filter_blast -h
    ```
    Version: 2.0 (lishenghui@genomics.org.cn, lishenghui1005@gmail.com)...
         3.0 (2011-6-3) Fixed some bugs...
         3.1 (2011-6-24) Combined with "connect_blast v2.0"...
         3.2 ...... Still some incompatible with "connect_blast v2.0"...
         4.0 (2012-11-20) Fixed some bugs...

    Usage: filter_blast <option> <value>...
    + Input the BLAST (m8) result...
        --input (-i) <file>
        --output (-o) <file>
    + Query ans Subject length filter (by percentage)...
        --qfile <file>
        --qper  <int> (80)
        --sfile <file>
        --sper  <int> (80)
    + Query and Subject length cutoff (for specified length of Q or S)...
        --qlength <int> (1)
        --slength <int> (1)
    + Identity, Evalue and Score cutoff...
        --identity <float> (9.0)
        --evalue   <float> (1e-1)
        --score    <float> (10.0)
    + Filter the Evalue or Score far from the best...
        --tope <int> 10^(99)  : top evalue times (10^1 is good)
        --tops <int> (100)    : top score percentage (10~20% is good)
        --wins <int> (9999)   : win score (the best NO. of scores)
    + Print Help...
        --help

    Example:
        (identity >= 95 and query length >= 90%): filter_blast -i f.blast -o f.blast.filter --identity 95 --qfile f.len --qper 90
        (evalue >= 1e-10 and top 10^1 times):     filter_blast -i f.blast -o f.blast.filter --evalue 1e-10 --tope 1
        (score >= 80 and top 20% of scores):      filter_blast -i f.blast -o f.blast.filter --score 80 --tops 20
    ```
#### filter_hmm2.pl
* Based on the alignment results from hmmsearch, this Perl script calculates the number of genes belonging to BUSCOs on each viral contig. It is used to filter viral contigs with high contamination of host genes (e.g., BUSCO rate ≥ 0.05).
    > \>filter_blast -h
    ```
    Usage: perl filter_hmm2.pl [*.hmm] [scores.cutoff] [*.hmm.out] [*.faa] [*.faa.busco]
    ```
    > \>hmmsearch --tblout vs.busco --cpu 10 --noali busco.bacteria_odb10.hmm vs.faa > $id.busco.log<br>
    > \>filter_hmm2.pl vs.busco [scores_cutoff](https://github.com/RChGO/VMGC/tree/main/Documents/scores_cutoff) vs.busco.f vs.faa vs.busco.f2<br>
    > \>awk '$3/$2<0.05' vs.busco.f2 > vs.busco.final

#### fasta_length & get_pep.length.pl<br>
* Two scripts are used respectively to calculate the lengths of nucleotide and protein sequences.
    > \>fasta_length -i uncl.fa -o uncl.fa.len<br>
    > \>get_pep.length.pl prot.fa prot.fa.len<br>

<br>

---

<font size=2> Correspondence and requests for materials should be addressed to grchun@hotmail.com </font>


