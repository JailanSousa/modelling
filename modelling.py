import re
import urllib
import os
import pandas as pd
import glob

from tqdm import tqdm
from modeller import *
from modeller.automodel import *

from best_template import best_hits


class Modelling:

    """This class take blast search and peforme a comparative with modeller."""
    
    def __init__(self, tfile) -> None:
        
        self.tfile = tfile
        self.hits = best_hits(self.tfile)
        self.BASE_DIR = f'{os.getcwd()}/modeller_results'

        if not os.path.exists(f"{self.BASE_DIR}"):
            os.mkdir(f"{self.BASE_DIR}")

    def get_template(self):

        """Select the best templates based on blast metrics (cuttoff)"""

        bad_templates = []
        templates = []
        for i in range(len(self.hits)):

            evalue = self.hits.evalue[i]
            pident = self.hits.pident[i]
            qcov = self.hits.qcov[i]
            
            # Identinty greater than 25%
            # Cover equal or greater than 70%
            # evalue equal or lower than 0.0001
            cutoff = pident >= 25 and qcov >= 70 and evalue < 0.0001
        
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
        
        """Comparative modelling. create 5 models."""

        templates = self.get_template()[0]
        bad_template = self.get_template()[1]
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
                a.ending_model = 5
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

def run(file):

    modelling = Modelling(file)
    modelling.download_templates()
    modelling.best_model()

tfile = 'blastp_tst_out.csv'         

run(tfile)