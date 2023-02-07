# LAMPgRNAtor

### LAMP and Cas12/13 gRNA generator

LAMPgRNAtor designs LAMP Primers and Cas12 gRNAs using GLAPD LAMP designing principles and Kim's CNN model respectively. Alternatively, it accepts GLAPD, LAVA and customized LAMP Primers as inputs to analyze its specificity and sensitivity. LAMPgRNAtor will search and predict for Cas12 gRNAs between F1c and B1c primers,  select the best LAMP Primer sets based on conservation and entropy scores, and plots its relative position and conservation score. It is available on all OS and designed with ease and user-friendliness in mind. 

Installation
------------

Simple Installation 
```bash
git clone https://github.com/chewlabSB2/LAMPgRNAtor.git
cd LAMPgRNAtor 
python setup.py build_ext install
```

Advised to use an environment as dependencies might clash with other software:

Installation via Anaconda
```bash
conda create -n LAMPgRNAtor python=3.8
conda activate LAMPgRNAtor

git clone https://github.com/chewlabSB2/LAMPgRNAtor.git
cd LAMPgRNAtor 
python setup.py build_ext install
```

Usage
-----
### LAMPgRNAtor-main 

With reference to the test dataset (MERS-CoV)

```bash
mafft --add NonMSAMERS.fasta --keeplength --reorder --nomemsave --thread $thread referenceMERS.fasta 1> MSAMERS.fasta 2> MAFFT_error.log
LAMPgRNAtor-main -r referenceMERS.fasta -bc NonMSAMERS.fasta -m MSAMERS.fasta --threads 8 -p MERS 1>run1.o 2>run1.e
```

### LAMPgRNAtor-LAMP

In order to design LAMP primers of any reference: 

```bash
LAMPgRNAtor-LAMP -r referenceMERS.fasta -p MERS
```

Expected Output
---------------
1. plot directory: Plots LAMP Primers and gRNAs positions relative to the conservation and entropy score of each position
2. csv directory: Contains the final LAMP Primers and gRNAs including thorough analysis with complete conservation and entropy score at each position. Includes alternative mapping of Primers/gRNAs to the consensus sequence  
3. position_score.csv: Includes the nucleotide, individual entropy and conservation score of every nucleotide of the consensus sequence generated from the multiple sequence alignment as supplementary data 

Add the flag --keep_tmp to retain the temporary files 

Links to GLAPD & LAVA
---------------------
[GLAPD](https://github.com/jiqingxiaoxi/GLAPD)

[LAVA](https://github.com/pseudogene/lava-dna)
