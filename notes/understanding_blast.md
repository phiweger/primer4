## Understanding Blast

Q: Which information is contained in the backtrace string (BTOP)? Does it softclip?

Example construction:

Query  1   TTTATCGAGACGCGAGCAGCGAGCAGCAGAGCGACGAGCAGCGATTT  47
                || ||||||||||||||||||||||||||||||||||||   
Sbjct  47  ACGTCCGGGACGCGAGCAGCGAGCAGCAGAGCGACGAGCAGCGACGA  93

>qry
TTTATCGAGACGCGAGCAGCGAGCAGCAGAGCGACGAGCAGCGATTT 

>ref
ACGTCCGGGACGCGAGCAGCGAGCAGCAGAGCGACGAGCAGCGACGA

makeblastdb -in ref.fna -dbtype nucl

blastn -dust no -word_size 13 -evalue 100 -outfmt "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore nident btop sstrand" -query qry.fna -db ref.fna -out aln.m8

cat aln.m8
qry ref 97.436  39  1   0   6   44  6   44  7.89e-18    67.638  2AG36   plus

So the query is 47 nt long. The BTOP string however 2AG36 is shorter, ie it omits softclipping. Also not how the lengths only correspond to the aligned fraction.

blastn -task blastn-short -query qry.fna -subject ref.fna -outfmt '6 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore btop' -out aln2.m8 -word_size 7 -evalue 1000