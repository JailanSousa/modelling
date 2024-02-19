from modeller import *

env = Environ()
aln = Alignment(env)
mdl = Model(env, file='2YBR', model_segment=('FIRST:C','LAST:C'))
aln.append_model(mdl, align_codes='2YBRA', atom_files='2YBR.pdb')
aln.append(file='1CN2.ali', align_codes='1CN2')
aln.salign(max_gap_length=20,
            gap_function=True,
            similarity_flag=True) 
aln.write(file='1cn2-2ybr.ali', alignment_format='PIR')
aln.write(file='1cn2-2ybr.pap', alignment_format='PAP')