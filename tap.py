from DPMacro.TAPFeature import TAPExtractor
import os
from functools import reduce

def write_tap(infile:str):
    file_id=os.path.split(infile)[1].replace('.pdb','')
    tap=TAPExtractor('i')
    tap.transform(infile)
    o_dict={}
    o_dict['cdr_lenth']=int(tap.reduced_frame['CDR_lenth'].sum())
    o_dict['psh']=float(tap.reduced_frame['PSH'].sum())
    o_dict['ppc']=float(tap.reduced_frame['PPC'].sum())
    o_dict['pnc']=float(tap.reduced_frame['PNC'].sum())
    o_dict['SFvCSP']=reduce(lambda x,y:x*y,tap.reduced_frame['sc_SFvCSP'])

    o_dict['psh_r']=tap.frame['psh'].to_list()
    o_dict['ppc_r']=tap.frame['ppc'].to_list()
    o_dict['pnc_r']=tap.frame['pnc'].to_list()    
    # out=f'{file_id},{cdr_lenth},{psh},{ppc},{pnc},{SFvCSP}\n'
    # with open(outfile,'a') as f:
    #     f.write(out)
    return o_dict

if __name__=='__main__':
    import argparse
    import json
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    parser=argparse.ArgumentParser()
    parser.add_argument('-i',metavar='INPUT.pdb',required=True)
    parser.add_argument('-o',metavar='DIR/TO/OUTPUT',default='./tap')
    parser.add_argument('-s',metavar='STANDARD.csv',default='./sample/abfold.csv')
    args=parser.parse_args()

    if os.path.isdir(args.o):
        raise ValueError('OUTDIR already exists!')
    else:
        os.mkdir(args.o)
    outp=lambda x:os.path.join(args.o,x)

    o_dict=write_tap(args.i)


    with open(outp('tap.json'),'w') as f:
        json.dump(o_dict,f)

    bg=pd.read_csv(args.s)
    thresholds={}
    for colname in ['cdr_lenth', 'psh', 'ppc', 'pnc', 'SFvCSP']:
        c=bg[colname]
        thresholds[colname]=[]

        p0=np.min(c)
        p05=np.percentile(c,5)
        if colname not in ['ppc', 'pnc']:
            thresholds[colname].append(p0)
            thresholds[colname].append(p05)
        p95=np.percentile(c,95)
        thresholds[colname].append(p95)
        p100=np.max(c)
        thresholds[colname].append(p100)

        f,axs=plt.subplots(1,1)
        a=axs
        a.hist(c,bins=20,density=True,color=(26/255,134/255,163/255,0.6))

        ymin,ymax=a.get_ylim()
        if colname not in ['ppc', 'pnc']:
            a.vlines(p0,ymin=ymin,ymax=ymax,colors=(1,0,0,0.9),linestyles='--',label='100% threshold')
            a.vlines(p05,ymin=ymin,ymax=ymax,colors=(251/255,132/255,2/255,0.9),linestyles='--',label='95% threshold')
        a.vlines(p95,ymin=ymin,ymax=ymax,colors=(251/255,132/255,2/255,0.9),linestyles='--')
        a.vlines(p100,ymin=ymin,ymax=ymax,colors=(1,0,0,0.9),linestyles='--')
        a.vlines(o_dict[colname],ymin=ymin,ymax=ymax,colors=(14/255,91/255,118/255,0.9),label='your score')
        
        # plt.legend()
        a.set_title(colname)
        f.set_dpi(400)
        f.savefig(outp(f'{colname}.png'))
        plt.close()

    with open(outp('threshold.json'),'w') as f:
        json.dump(thresholds,f)

    
