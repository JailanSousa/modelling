import pandas as pd

def fasta(file):
    """generate a fasta file from a csv file."""
    
    pepdb = pd.read_csv(file)

    for entry, seq in zip(pepdb.Entry, pepdb.Sequence):
        
        lines = f">{entry}\n{seq}\n"
        
        with open('pepset.fasta', 'a') as writting:
            writting.write(lines)
        

file = '/home/jssousa/neuro_pep_dataset/pep_dataset_filted.csv'
fasta(file)
print('OK')