# Promoter Discovery in Bacteria using Statistical Gene Prediction
This aims at finding promoters within a given set of genomes of Escherichia coli bacteria using statistical gene prediction methods. This is performed as a partial requirement for the BM4321 - Genomic Signal Processing during the senior year at University of Moratuwa. The results from the analysis are presented in this [report](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/blob/main/Documentation/160005C_Assignment.pdf). 

![Pribnow Box Search](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/blob/main/Documentation/Promotor%20search.png)

The complete genomes of Escherichia coli bacteria in their FASTA format can be accessed in this [folder](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/tree/main/Genomes). The corresponding protein tables are provided in this [folder](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/tree/main/ProteinTables). You can access the genomes and protein tables by visiting [NCIB](https://www.ncbi.nlm.nih.gov/genome/) Genome Website. For further instruction on how to download the genomes and their corresponding protein tables, refer to [this](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#:~:text=To%20use%20the%20download%20service,button%20to%20start%20the%20download).

## Pseudo Algorithm for detecting Pribnow Box
...[1.] Read the protein tables to identify the *starting positions* of genes (proteins).
...[2.] Obtain 50 bases upstream and 3 bases downstream for all possible genes of the sense and antisense strands from their *starting positions*. 
...[3.] Perform a *standard local search (for an intact query)* to locate a Pribnow box promoter within upstream positions 5 to 30 of each sequence.
...[4.] Using the first 1000 sequences, obtain a position probability matrix (PPM) with 10 positions for the Pribnow box.
...[5.] Using a suitable entropy measure, eliminate the redundant positions of PPM obtained in [4.] to obtain *reduced PPM*. 
...[6.] Perform a statistical alignment for the remaining sequences of [2.] using the initial PPM of [4.] and the reduced PPM of [5.]. 
...[7.] For non-detected sequences, perform a *local search with a non-intact query* to locate any potential mutated promoters.


