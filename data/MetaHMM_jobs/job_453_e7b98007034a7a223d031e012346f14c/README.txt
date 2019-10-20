The following files can be found in this archive:

input_unaligned.fasta
	The original, unaligned protein sequences comprising the model, as a FASTA file.
input_aligned.sto
	The aligned sequences of the model, as a STOCKHOLM file.
input_profile.hmm
	The HMM profile, built on the aligned sequences of the model.
output_all.csv
	All the matching domains found in the metagenomes, combined into one file.
	Matches are listed by full sequence E-value ascending (best first).
output_unaligned.fasta
	All the matching domains found in the metagenomes, combined into one file.
	Only those with E-value < 1e-9 are included.
output_***.txt
	Detailed search results, one file for each metagenome involved in the search.
	The files are in the native output format of hmmsearch.
	Those sequences which matched the model are listed after each other, ordered by E-value ascending (best first).
README.txt
	This readme file.
