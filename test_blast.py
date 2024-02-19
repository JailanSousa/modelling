from Bio.Blast.Applications import NcbiblastpCommandline

fasta = 'pepset.fasta'

outfmt = "10 qseqid sseqid evalue bitscore pident qcovs qseq sseq"
blastp = NcbiblastpCommandline(task='blastp',
                               query=fasta,
                               db="/home/jssousa/pdbdb/pdbaa",
                               outfmt=outfmt,        
                               out='blastp_tst_out.csv'
                               )

stdout, stderr = blastp()

#print(stdout, stderr)