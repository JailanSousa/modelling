import re
import urllib
import os
import pandas as pd
from modeller import *
from modeller.automodel import *

from best_template import best_hits
from tqdm import tqdm

class Modelling:
    
    def __init__(self, tfile) -> None:
        
        self.tfile = tfile
        self.hits = best_hits(self.tfile)
        self.BASE_DIR = f'{os.getcwd()}/modeller_results'

        if not os.path.exists(f"{self.BASE_DIR}"):
            os.mkdir(f"{self.BASE_DIR}")

    def get_template(self):

        bad_templates = []
        templates = []
        for i in range(len(self.hits)):

            evalue = self.hits.evalue[i]
            pident = self.hits.pident[i]
            qcov = self.hits.qcov[i]
            
            cutoff = pident >= 25 and qcov >= 60 and evalue < 0.8
        
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
        #df.to_csv(f'{self.templates_dir}/bad_templates.csv', index=False)
        
        return pd.DataFrame(templates), df
    
    def download_templates(self):

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
        
        templates = self.get_template()[0]
        hits = self.hits
        pir_dir = f'{self.BASE_DIR}/ali'
        
        if not os.path.exists(pir_dir):
            os.mkdir(pir_dir)
        
        for i in range(len(hits)):
            
            for j in range(len(templates)):
                
                if hits.qseqid[i] == templates.qseqid[j]:
                    
                    f_name = hits.qseqid[i]
                    seq = re.sub('_', '', hits.qseq[i])
                    
                    ali = f">P1;{f_name}\n" + \
                        f"sequence:{f_name}:::::::0.00: 0.00\n" + \
                        f"{seq}*"
                    
                    with open(f'{pir_dir}/{f_name}.ali', 'w') as writing:
                        writing.write(ali)
        
    def modelling(self):
        
        templates = self.get_template()[0]
        tst = self.create_pir()
        
        template_dir = f"{self.BASE_DIR}/templates"
        pir_dir = f'{self.BASE_DIR}/ali'


        aligments_dir = f'{self.BASE_DIR}/aligments'
        if not os.path.exists(aligments_dir):
            os.mkdir(aligments_dir)

        models_dir = f'{self.BASE_DIR}/models'
        if not os.path.exists(models_dir):
            os.mkdir(models_dir)
        
        outputs = []
        for i in range(len(templates)):
            
            pdb_id = templates.pdb_id[i]
            qseqid = templates.qseqid[i]
            chain = templates.chain[i]
            
            if os.path.exists(f'{template_dir}/{pdb_id}.pdb'):

                ali = f'{pir_dir}/{qseqid}.ali'
                template = f'{template_dir}/{pdb_id}.pdb'
                    
                print('template found')

                # Alignment 
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

                # Move outputs for models directory
                os.system(f'mv {qseqid}* /{models}')
            

        return pd.DataFrame(outputs)
    
    def best_model(self):

        outputs = self.modelling()
        code = [''.join(re.findall(r'.+(?=\.B)', x)) for x in outputs.name]
        outputs['Code'] = code
        failure = outputs[outputs.failure == None]
        best = outputs.loc[outputs.groupby('Code')['DOPE score'].idxmin(), :]
        print(best)

        """best_models_dir = f'{self.BASE_DIR}/best_models'"""

    """ if not os.path.exists(best_models_dir):
            os.mkdir(best_models_dir)"""

        #models_dir = f'{self.BASE_DIR}/models/*'

        # Get a list of all successfully built models from a.outputs
    """ ok_models = [x for x in a.outputs if x['failure'] is None]
        failures = [x for x in a.outputs if x['failure'] is not None]"""

        # Rank the models by DOPE score
    """key = 'DOPE score'
        ok_models.sort(key=lambda a: a[key])"""

        # Get top model
    """m = ok_models[0]
        if os.path.exists(f'{models_dir}/*.pdb'):"""



tfile = 'blastp_tst_out.csv'         
modelling = Modelling(tfile)
modelling.download_templates()
modelling.best_model()