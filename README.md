# LAMPgRNAtor

## Still Fixing Minor Bugs (Not Completed)

Installation
------------

Simple Installation (Advised to Use Conda):
```bash
git clone https://github.com/chewlabSB2/LAMPgRNAtor.git
cd LAMPgRNAtor 
python setup.py install
```

Installation via Conda
```bash
conda create -n LAMPgRNAtor python=3.8
conda activate LAMPgRNAtor

git clone https://github.com/chewlabSB2/LAMPgRNAtor.git
cd LAMPgRNAtor 
python setup.py install
```

Links to GLAPD & LAVA
---------------------
[GLAPD](https://github.com/jiqingxiaoxi/GLAPD)
[LAVA](https://github.com/pseudogene/lava-dna)

Usage
-----
###LAMPgRNAvalidate 
With reference to the test dataset (SARS-CoV-2)
```bash
mafft --add NonMSAEV71.fasta --keeplength --reorder --nomemsave --thread $thread referenceEV71.fasta 1> MSAEV71.fasta 2> error.log
LAMPgRNAvalidate -r referenceSARS2.fasta -bc NonMSASARS2.fasta -m MSASARS2.fasta -bo Offtarget.fasta -p test -l GLAPD_sample.txt LAVA_sample.txt 1>run1.o 2>run1.e
```

Expected Output
---------------
1.
2.
3. position_score.csv: Includes the nucleotide, individual entropy and conservation score of every nucleotide of the consensus sequence generated from the multiple sequence alignment as supplementary data 
4.


Add the flag --keep_tmp to retain the temporary files 

Citation 
--------
If you use LAMPgRNAtor please cite: [Will Be Updated]()