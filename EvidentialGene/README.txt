(base) bioinformatica@bioinformatica-MOOVE3-14:/media/bioinformatica/SEAGATE_GFG/8-TFM-2025/Bge_Transcriptome/Spa+Try_over1TPM_50c_perLib$ 
$ls
Tri_ctgs_TPMover1_subgr.txt  Spa_ctgs_TPMover1_subgr.txt

# Filtramos los ensamblajes originales
$seqtk subseq trinity_out_dir.Trinity.fasta Tri_ctgs_TPMover1_subgr.txt > over1sg_Trinity_tx.fasta
$seqtk subseq RNASpades_transcripts.fasta Spa_ctgs_TPMover1_subgr.txt > over1sg_RNASpades_tx.fasta

(base) gfg017@gfg017-A7:/media/gfg017/SEAGATE_GFG/8-TFM-2025/REPETICION_18.08/Bge_Transcriptome/Spa+Try_over1TPM_perSubGroup
$ wc -l *gr.txt
  62054 Spa_ctgs_TPMover1_subgr.txt
 153064 Tri_ctgs_TPMover1_subgr.txt
 215188 total

 
# Comprobamos tamaños de los ensamblajes filtrados
$ cat over1sg_RNASpades_tx.fasta | grep "^>" | wc -l
62054
$ cat over1sg_Trinity_tx.fasta | grep "^>" | wc -l
153064

# Fusionamos ambos ensamblajes y fusionamos
$ cat over1sg_* > over1_sg_Merged_tx.fasta
$ cat over1_sg_Merged_tx.fasta | grep "^>" | wc -l
215118


He visto en google forums que es un parametro a modificar,
si tenemos como referencia el HIFI de P.americana, de los
37240 cds, todos excepto 3 son >= 210 (70aa):

metaquast.py --threads 8 GCF_040183065.1_P.americana_PAMFEO1_priV1_cds_from_genomic.fna --max-ref-num 0 --contig-thresholds  0,90,120,150,210,240,300,320,400,500,600,700,800,900,1000,1500,2000,2500,5000,10000,25000,50000 


# Filtramos a 65 aa longitud
/home/gfg017/Documents/9-TFM-2026/Tools/arthropods.eugenes.org/EvidentialGene/evigene/scripts/prot/tr2aacds4.pl -MINAA 65 -strandedrna yes -cdnaseq over1_sg_Merged_tx.fasta -species=Bge
#t2ac: EvidentialGene tr2aacds.pl VERSION 2022.04.05
#t2ac: CMD: tr2aacds.pl  -MINAA 65 -strandedrna yes -cdnaseq over1_sg_Merged_tx.fasta -species=Bge
#t2ac: NOTE: evd_smallclassfilter not enough evidence nev=9718 for naa=109684
# Class Table for over1_sg_Merged_tx.trclass 
class    	%okay	%drop	okay	drop
althi    	4.7	0.5	5210	575
althi1   	17.3	2.7	19068	3004
althinc  	1.6	0	1812	0
altmfrag 	0.9	0.28	1075	312
altmid   	0.8	0.6	985	677
main     	12.3	0	13634	2
mainnc   	1.9	0	2142	0
noclass  	13	0	14349	1
noclassnc	10.6	0	11703	0
parthi   	0	4.6	0	5121
parthi1  	0	2.9	0	3284
perfdupl 	0	22	0	24244
perffrag 	0	2.5	0	2754
smallorf 	0	0	0	0
---------------------------------------------
total    	63.6	36.3	69978	39974
=============================================
# AA-quality for okay set of over1_sg_Merged_tx.aa.qual (no okalt): all and longest 1000 summary 
okay.top	 n=1000; average=1991; median=1701; min,max=1267,12509; nfull=978; sum=1991715; gaps=68,0
okay.all	 n=41828; average=265; median=107; min,max=65,12509; nfull=28980; sum=11105713; gaps=272,0
#ta2c: ERR:output skip dup id TRINITY_DN178764_c0_g1_i1 
#readPubidTab(publicset/over1_sg_Merged_tx.pubids)= 69978
# nin=220714, nok=144721, nfrag=383, nskipnotloc=64737, nskipdupfrag=3284, nskipdiffloc=7589
#insertUniqExons= 1658
#collectExonChains= 69199 of 69617 ids
#assignChainLoci
#n_class: ichain=54350 icalt=11140 icsub=3709 icdup=418


