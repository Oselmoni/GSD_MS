reformatFastQ_v2.py

This script reformats a read file by: 1) Removing the first ten bases from each read and 2) reformatting the sequences header (remove space and add /1 or /2 subfix) needed for correct implementation of a fastQ in Trinity.
The input file must be a .fastq file, the read sequence must be on ONE LINE. More info about the standard Illumina output format for fastq files: http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm
To launch from the command line: $ python reformatFastQ_v2.py readfile.fastq

transcriptome_stats.py

This script computes count of gene and isoform per genes on the Trinity.fasta file, produced by the Trinity assembler.
The scripts uses the standard header of the Trinity assembler to sort isoforms by gene. The stats produced are: nb of transcripts, nb of genes, nb of gene with 1 to 15 or more isoforms.
To launch from the command line: $ python transcriptome_stats.py Trinity.fasta

parseannotation.py

The script filters out the transcripts from the Trinity.fasta assembly that do not have any match in similarity search against the swissport database.
The similarity search is previously done using the blastX command line application, with the following flags on:
$ blastx -db swiss -query Trinity.fasta -out annotation.txt -evalue 0.000001 -outfmt 6 -max_target_seqs 1
To launch from the command line: $ python parseannotation.py <annotation.txt> <Trinity.fasta>

parseDrerioannotation.py

The custom database can be downloaded directly from the uniprot website. Set filters of interest in the search and download the results.
Two files are required:
1) the sequences (download>fasta (canonical) > D.Rerio.fasta
2) the table of the attributes associated to the sequences containing: Entry,Entry name,Protein names,Gene ontology (biological process). Choose these four attributes in the results output formatting, then download with download>tab separated > Drerio_features.txt
Then blast filtered sequences in Assembly_filtered.fasta against the custom DB. 
$ blastx -subject <DRerio.fasta> -query <Assembly_Filtered.fasta> -out <annotationDR.txt> -evalue 0.000001 -outfmt 6 -max_target_seqs 1
Retrieve GO terms associated to the matches using:
$ python parseDrerioannotation.py <annotationDR.txt> <Drerio_features.txt>
This will produce GOs.R, ready to be used in DGEA in R.

tran2gene.py

Reformats the output of Kallisto by calculating the sum of the estimated counts per gene (Kallisto provides the estimated counts per each transcript).
The sorting is performed taking advantage of the Trinity header format (similarly to the transcriptome_stats.py script).
To launch from the command line: $ python tran2gene.py abundance.txt

customFunctions.R

Contains the voomWithQualityWeightsMOD and the voomMod function used during the normalization prior to the differential gene expression analysis.
The voomMod function was modified from the voom code by Julien Roux, and performs the normalization to cpm using the function from the edgeR package.
The voomWithQualityWeightsMOD calls the voomMod function.


