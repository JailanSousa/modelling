import warnings
warnings.filterwarnings("ignore")

from Bio.Blast.Applications import NcbiblastpCommandline

def search(file):
    """blastp search."""
    
    # fields returned in output file
    outfmt = "10 qseqid sseqid evalue bitscore pident qcov qseq sseq"
                         
    blastp = NcbiblastpCommandline(task='blastp',
                                   query=file,
                                   db="/home/jssousa/pdbdb/pdbaa",
                                   outfmt=outfmt,        
                                   out=f'pep_blastp_search.csv'
                                   )
    print('\nRunning....\n')
    stdout, stderr = blastp()


search('pepset.fasta')

print('blastp_finished')
