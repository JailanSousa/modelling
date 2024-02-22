import pandas as pd


def best_hits(file):

    """Select the best blast results based on blast output metrics."""
    
    blast = pd.read_csv(file)

    comboscore = blast['pident'] + blast['qcov'] + blast['bitscore'] - blast['evalue']

    blast['combo'] = comboscore

    hits = blast.loc[blast.groupby('qseqid')['combo'].idxmax(), :]

    hits.drop('combo', axis=1, inplace=True)
    return hits.reset_index(drop=True)
