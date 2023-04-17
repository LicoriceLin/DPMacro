from DPMacro import FvStructureProcesser as fv
from DPMacro import util as u
from typing import Optional,Dict,Tuple
from DPMacro.Data import CDR_annotations,CDR_maps
id=Tuple[str,int,str]
idmap=Dict[id,id]


def cdr_tag(resmap:idmap,scheme):
    '''
    chain already renumbered by the given scheme.
    '''
    s=CDR_annotations[scheme]
    cdr3b=resmap[(' ',s['CDR3'],' ')][1]
    cdr3e=resmap[(' ',s['FR4']-1,' ')][1]
    return cdr3b,cdr3e



import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-i',metavar='INPUT.pdb',required=True,type=str)
parser.add_argument('-a',metavar='Antigens',default='A,B',type=str)
parser.add_argument('-hc',metavar='Hchain',default='H',type=str)
parser.add_argument('-lc',metavar='Lchain',default='L',type=str)
parser.add_argument('-s',metavar='Scheme',default='i',type=str)
parser.add_argument('-o',metavar='OUT.pdb',default='out.pdb',type=str)
parser.add_argument('-keep-ag',action='store_true')
parser.add_argument('-keep-scheme',action='store_true')
args=parser.parse_args()

file:str=args.i
ags:str=args.a
hchain:str=args.hc
lchain:str=args.lc
scheme:str=args.s
out:str=args.o
keep_ag:bool=args.keep_ag
keep_scheme:bool=args.keep_scheme

s=u.read_in(file)
f=fv(s)
f.parse_structure(scheme)

m=u.Model(0)
if keep_ag:
    for ag in ags.split(','):
        m.add(f.object[0][ag].copy())
else:
    for ag in ags.split(','):
        m.add(f.splited_chains[ag]['Ag0'].copy())
        m['Ag0'].id=ag
    
# u.extract_hetatm(m['A'],inplace=True)
m.add(f.splited_chains[hchain]['Fv0'].copy())
m['Fv0'].id='H'
m.add(f.splited_chains[lchain]['Fv0'].copy())
m['Fv0'].id='L'

if not keep_ag:
    for ag in ags.split(','):
        a,_=u.renum(m[ag])
        u.add_child(m,a,ag)
        
if not keep_scheme:
    h,_h=u.renum(m['H'])
    u.add_child(m,h,'H')
    l,_l=u.renum(m['L'])
    u.add_child(m,l,'L')



u.write_out(m,out)
if not keep_scheme:
    h3anot=cdr_tag(_h,scheme)
    l3anot=cdr_tag(_l,scheme)
    with open(out,'a') as f:
        f.write(f'{h3anot[0]},{h3anot[1]}')
else:
    b=CDR_annotations[scheme]['CDR3']
    e=CDR_annotations[scheme]['FR4']-1
    with open(out,'a') as f:
        f.write(f'{b},{e}')