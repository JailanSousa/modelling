import warnings
warnings.filterwarnings("ignore")

import re
import sys
import urllib
import os
import pandas as pd
import glob

from Bio.Blast.Applications import NcbiblastpCommandline
from tqdm import tqdm
from modeller import *
from modeller.automodel import *

tempalte_db_path = '/home/jssousa/mestrado/modeling_pep/pdbaa/pdbaa'
class Modelling:

    """This class take blast search and peforme a comparative with modeller."""
    
    def __init__(self, fasta='', pident=70, qcov=70,
                evalue=0.0001, models=5, db=tempalte_db_path) -> None:
        
        self.fasta = fasta
        self.pident = int(pident)
        self.qcov = int(qcov)
        self.evalue = float(evalue)
        self.models = int(models)
        self.db = db
        #self.tfile = tfile
        self.hits = None #best_hits(self.tfile)
        self.BASE_DIR = f'{os.getcwd()}/modeller_results'

        if not os.path.exists(f"{self.BASE_DIR}"):
            os.mkdir(f"{self.BASE_DIR}")
    
    def search_template(self):
        """blastp search."""
        
        # fields returned in output file
        outfmt = "10 qseqid sseqid evalue bitscore pident qcovs qseq sseq"
                         
        blastp = NcbiblastpCommandline(task='blastp',
                                   query=self.fasta,
                                   db=self.db,
                                   outfmt=outfmt,        
                                   out=f'blastp_search.csv'
                                   )
        print('\nRunning....\n')
        stdout, stderr = blastp()

    def best_hits(self):

        """Select the best blast results based on blast output metrics."""
        colnames = ['qseqid', 'sseqid','evalue', 
                    'bitscore', 'pident', 'qcovs', 'qseq', 'sseq']
        self.search_template()
    
        blast = pd.read_csv('blastp_search.csv', names=colnames, header=None)

        comboscore = blast['pident'] + blast['qcovs'] + blast['bitscore'] - blast['evalue']

        blast['combo'] = comboscore

        hits = blast.loc[blast.groupby('qseqid')['combo'].idxmax(), :]

        hits.drop('combo', axis=1, inplace=True)
        return hits.reset_index(drop=True)


    def get_template(self):
        self.hits = self.best_hits()

        """Select the best templates based on blast metrics (cuttoff)"""

        bad_templates = []
        templates = []
        for i in range(len(self.hits)):

            evalue = float(self.hits.evalue[i])
            pident = int(self.hits.pident[i])
            qcov = int(self.hits.qcovs[i])
            
            # Identinty greater than 25%
            # Cover equal or greater than 70%
            # evalue equal or lower than 0.0001
            cutoff = pident >= self.pident and qcov >= self.qcov and evalue < self.evalue
        
            if cutoff:
            
                pdb_id = self.hits.sseqid[i].split('_')
    
                templates_info = {'qseqid': self.hits.qseqid[i],
                                  'pdb_id': pdb_id[0],
                                  'chain': pdb_id[1]
                                  }
                
                templates.append(templates_info)
            
            else:
                
                templates_info = {'qseqid': self.hits.qseqid[i],
                                  'sseqid': self.hits.sseqid[i],
                                  'evalue': evalue,
                                  'pident': pident,
                                  'qcov': qcov
                                  }
                
                bad_templates.append(templates_info)
        
        df = pd.DataFrame(bad_templates)
        
        return pd.DataFrame(templates), df
    
    def download_templates(self):

        """Download templates in pdb format."""
        
        # Template directory
        templates_dir = f"{self.BASE_DIR}/templates"
        if not os.path.exists(templates_dir):
            os.mkdir(templates_dir)

        templates = self.get_template()[0]
        
        for i in range(len(templates)):
            
            pdb_id = templates.pdb_id[i].split('_')
        
            try:
                     
                URL = f'https://files.rcsb.org/download/{pdb_id[0]}.pdb'
                os.system(f'wget -P {templates_dir} {URL}')
                    
            except (urllib.error.HTTPError, urllib.error.URLError):
                
                log_msg = 'pdb file not found (download: error 404), unmodeled'
                print('urllib.error.HTTPError: HTTP Error 404: Not Found.')
                self.__writelog(log_msg, cod=pdb_id[0])
        
    def create_pir(self):

        """Create file which is modeller input file."""
        
        templates = self.get_template()[0]
        hits = self.hits
        
        # Pir directory
        pir_dir = f'{self.BASE_DIR}/ali'
        if not os.path.exists(pir_dir):
            os.mkdir(pir_dir)
        
        for i in range(len(hits)):
            
            for j in range(len(templates)):
                
                if hits.qseqid[i] == templates.qseqid[j]:
                    
                    f_name = hits.qseqid[i]
                    seq = re.sub('_', '', hits.qseq[i])
                    
                    # Pir file format
                    ali = f">P1;{f_name}\n" + \
                        f"sequence:{f_name}:::::::0.00: 0.00\n" + \
                        f"{seq}*"
                    
                    with open(f'{pir_dir}/{f_name}.ali', 'w') as writing:
                        writing.write(ali)
        
    def modelling(self):
        
        """Comparative modelling."""

        templates = self.get_template()[0]
        bad_template = self.get_template()[1]
        self.download_templates()
        self.create_pir()
        
        template_dir = f"{self.BASE_DIR}/templates"
        pir_dir = f'{self.BASE_DIR}/ali'

        # Aligment directory
        aligments_dir = f'{self.BASE_DIR}/aligments'
        if not os.path.exists(aligments_dir):
            os.mkdir(aligments_dir)

        # Models directory
        models_dir = f'{self.BASE_DIR}/models'
        if not os.path.exists(models_dir):
            os.mkdir(models_dir)
        
        outputs = []
        for i in tqdm(range(len(templates)), desc='Modelling'):
            
            pdb_id = templates.pdb_id[i]
            qseqid = templates.qseqid[i]
            chain = templates.chain[i]
            
            if os.path.exists(f'{template_dir}/{pdb_id}.pdb'):

                ali = f'{pir_dir}/{qseqid}.ali'
                template = f'{template_dir}/{pdb_id}.pdb'
                    
                print('template found')

                # Alignments 
                env = Environ()
                aln = Alignment(env)
                mdl = Model(env, file=template, 
                            model_segment=(f'FIRST:{chain}',f'LAST:{chain}'))
                aln.append_model(mdl, align_codes=f'{pdb_id}{chain}',
                                 atom_files=template)
                aln.append(file=ali, align_codes=qseqid)
                aln.salign(max_gap_length=20,
                            gap_function=True,
                            similarity_flag=True) 
                aln.write(file=f'{aligments_dir}/{qseqid}_{pdb_id}.ali', 
                          alignment_format='PIR')
                aln.write(file=f'{aligments_dir}/{qseqid}_{pdb_id}.pap', 
                          alignment_format='PAP')
                
                # Modelling
                # Query directory inside models directory
                models = f'{models_dir}/{qseqid}'
                if not os.path.exists(models):
                    os.mkdir(models)

                env = Environ()
                a = AutoModel(env, alnfile=f'{aligments_dir}/{qseqid}_{pdb_id}.ali',
                            knowns=f'{pdb_id}{chain}', sequence=qseqid,
                            assess_methods=(assess.DOPE, assess.GA341))
                
                a.starting_model = 1
                a.ending_model = self.models
                a.make()
                outputs.extend(a.outputs)

                # Move outputs to models directory
                os.system(f'mv {qseqid}* /{models}')
            
        bad_template.to_csv(f'{self.BASE_DIR}/bad_templates.csv', index=False)
        return pd.DataFrame(outputs)
    
    def best_model(self):

        """Select the best model based on  
        Discrete Optimized Protein Energy (DOPE score)"""
        outputs = self.modelling()
        code = [''.join(re.findall(r'.+(?=\.B)', x)) for x in outputs.name]
        outputs['Code'] = code
        lowestdope = outputs.loc[outputs.groupby('Code')['DOPE score'].idxmin(), :]
        best = lowestdope.reset_index(drop=True)
        failure = outputs[outputs.failure == None]

        # Best models directory
        best_models_dir = f'{self.BASE_DIR}/best_models'
        if not os.path.exists(best_models_dir):
            os.mkdir(best_models_dir)

        # move the best models from models to best_model directory
        models_dir = f'{self.BASE_DIR}/models/*/*.pdb'
        filer_names = glob.glob(models_dir)

        for i in range(len(best)):
            for f in filer_names:
                file = ''.join(re.findall(r'.*/([^/]+)$', f))

                if best.name[i] == file:
                    os.system(f'mv {f} /{best_models_dir}')

        # Register models failures.
        failure.to_csv(f'{self.BASE_DIR}/models_failure.csv', index=False)

    def run(self):
        
        """Set flags and run modelling"""
        argument = ''
        try:
            for i, arg in enumerate(sys.argv):
                if i % 2 == 1:
                    argument = arg
                    if argument == '-h' or argument == '--help':
                        print("""HELP:
                              -i --input file name in fasta or multfasta format.
                              -id --identity identity cutoff (default 70%)
                              -c --cover coverence cutoff (default 70%)
                              -e --evalue e-valiue cutoff (default 0.0001)
                              -m --models models quantity to be generated
                              -db --database path to the templates structure
                              """)
                        break
                elif i % 2 == 0 and i != 0:
                    value = arg
                    if argument == '-i' or argument == '--input':
                        self.fasta = value
                    elif argument == '-id' or argument == '--identity':
                        self.pident = int(value)
                    elif argument == '-c' or argument == '--cover':
                        self.qcov = int(value)
                    elif argument == '-e' or argument == '--evalue':
                        self.evalue = float(value)
                    elif argument == '-m' or argument == '--models':
                        self.models = int(value)
                    elif argument == '-db' or argument == '--database':
                       self.db = value
        except: 
            print('Look the instructions')
        
        if argument != '-h' and argument != '--help':
            self.best_model()
            os.system(f"rm -r {self.BASE_DIR}/templates")
            os.system(f"rm -r {self.BASE_DIR}/ali")
            os.system(f"rm -r {self.BASE_DIR}/aligments")
            os.system(f"rm -r {self.BASE_DIR}/models")

            


modelling = Modelling()
modelling.run()