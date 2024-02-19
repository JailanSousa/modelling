
def create_pir(self):
    """ 
    put the target sequence into the 
    PIR format readable by MODELLER.
    """
    with open(f'file.ali', 'a') as writing:
        writing.write(f'>P1;pep{self.index}\n //
                        sequence:pep{self.index}:::::::0.00: 0.00\n //
                        {self.peptide}*')
