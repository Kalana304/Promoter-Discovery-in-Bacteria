# Promoter Discovery in Bacteria using Statistical Gene Prediction
This aims at finding promoters within a given set of genomes of Escherichia coli bacteria using statistical gene prediction methods. This is performed as a partial requirement for the BM4321 - Genomic Signal Processing during the senior year at University of Moratuwa. The results from the analysis are presented in this [report](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/blob/main/Documentation/Final%20Report.pdf). 

![Pribnow Box Search](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/blob/main/Documentation/Promotor%20search.png)

The complete genomes of Escherichia coli bacteria in their FASTA format can be accessed in this [folder](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/tree/main/Genomes). The corresponding protein tables are provided in this [folder](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/tree/main/ProteinTables). You can access the genomes and protein tables by visiting [NCIB](https://www.ncbi.nlm.nih.gov/genome/) Genome Website. For further instruction on how to download the genomes and their corresponding protein tables, refer to [this](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#:~:text=To%20use%20the%20download%20service,button%20to%20start%20the%20download).

## Pseudo algorithm for detecting Pribnow Box

**Inputs :** Protein Table, Complete Genome (for the corresponding accession number) \
**Outputs :** Pribnow Box locations

1. Read the protein tables.
    > - The protein table corresponding to the accession number of the genome will contain the details on the starting and ending positions of the genes within the complete genome. This step identify the *starting positions* of genes (proteins). 
2. Obtain possible sequences where Pribnow Box is located.
    > - From the starting position of each protein in the protein table, obtain 50 bases upstream and 3 bases downstream for all possible genes of the sense and antisense strands. 
3. Obtain a *Position Probability Matrix (PPM)*. 
   > - Perform a *Standard Local Search - Intact Query* to locate a `Pribnow Box` within upstream positions 5 to 30 of each sequence.
   > - The Local Search will return the starting points of possible `Pribnow Box` in each sequence. 
   > - Using the first 1000 results, obtain a position probability matrix (PPM) with 10 positions for the `Pribnow Box`.
4. Obtain *Reduced PPM*. 
   > - Using a suitable entropy measure, eliminate the redundant positions of PPM obtained in [3.]. 
5. Statistical alignment for the remaining sequences of [2.] using the initial PPM of [4.] and the reduced PPM of [5.] to detect `Pribnow Box` positions. 
6. For non-detected sequences, perform a *local search with a non-intact query* to locate any potential mutated promoters.


**Same pseudo code is followed to detect the presence of `TTGACA Box`.**

## Executing the program

The scripts provided are developed using **SciLab - 6.1.0(x64)**. You can freely download it by visiting this [website](https://www.scilab.org/).

To run the program:
> execute `main.sce` in the [Scripts](https://github.com/Kalana304/Promoter-Discovery-in-Bacteria/tree/main/Scripts).

To change the protein table/genome:
> Change L34 and L35 of `main.sce`. 





