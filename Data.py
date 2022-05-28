'''
library of constant,

'''
# from Bio.Data.IUPACData import protein_letters_3to1,protein_letters_1to3
from Bio.Data import SCOPData


#3-letter code <-> 1-letter code
amino1to3dict={
                  'A':'ALA',
                  'R':'ARG',
                  'N':'ASN',
                  'D':'ASP',
                  'C':'CYS',
                  'E':'GLU',
                  'Q':'GLN',
                  'G':'GLY',
                  'H':'HIS',
                  'I':'ILE',
                  'L':'LEU',
                  'K':'LYS',
                  'M':'MET',
                  'F':'PHE',
                  'P':'PRO',
                  'S':'SER',
                  'T':'THR',
                  'W':'TRP',
                  'Y':'TYR',
                  'V':'VAL',
                  'X':'XXX'}

nucleic_acid_dict = {
    'A':'A',
    'C':'C',
    'G':'G',
    'U':'U',
    'T':'T'
}
amino3to1dict = {
    **SCOPData.protein_letters_3to1,
    **{'ASH': 'A','CYX': 'C','HYP': 'P',
        'HID': 'H','HIE': 'H','HIP':'H','MSE': 'M'}
}


#scheme information
base_CDR_annotations = {'imgt':{
                        'FR1':1,
                        'CDR1':27,
                        'FR2':39,
                        'CDR2':56,
                        'FR3':66,
                        'CDR3':105,
                        'FR4':118,
                        'end':130
                    },
                    'aho':{
                        'FR1':1,
                        'CDR1':25,
                        'FR2':41,
                        'CDR2':58,
                        'FR3':78,
                        'CDR3':109,
                        'FR4':138,
                        'end':149 
                    },
                    'chothia|H':{
                        'FR1':1,
                        'CDR1':26,
                        'FR2':33,
                        'CDR2':52,
                        'FR3':57,
                        'CDR3':95,
                        'FR4':103,
                        'end':113 
                    },
                    'chothia|L':{
                        'FR1':1,
                        'CDR1':24,
                        'FR2':35,
                        'CDR2':50,
                        'FR3':57,
                        'CDR3':89,
                        'FR4':98,
                        'end':113 
                    },
                    'kabat|H':{
                        'FR1':1,
                        'CDR1':31,
                        'FR2':36,
                        'CDR2':50,
                        'FR3':66,
                        'CDR3':95,
                        'FR4':103,
                        'end':113 
                    },
                    'kabat|L':{
                        'FR1':1,
                        'CDR1':24,
                        'FR2':35,
                        'CDR2':50,
                        'FR3':77,
                        'CDR3':89,
                        'FR4':98,
                        'end':107 
                    },
                    'kabat|K':{
                        'FR1':1,
                        'CDR1':24,
                        'FR2':35,
                        'CDR2':50,
                        'FR3':77,
                        'CDR3':89,
                        'FR4':98,
                        'end':108 
                    },
             }

CDR_annotations = {
    **base_CDR_annotations,
    'a':base_CDR_annotations['aho'],
    'i':base_CDR_annotations['imgt'],
    'chothia|K':base_CDR_annotations['chothia|L'],
    'c|H':base_CDR_annotations['chothia|H'],
    'c|L':base_CDR_annotations['chothia|L'],
    'c|K':base_CDR_annotations['chothia|L'],
    'k|H':base_CDR_annotations['kabat|H'],
    'K|L':base_CDR_annotations['kabat|L'],
    'k|K':base_CDR_annotations['kabat|K']
             }


#static residue property
hydrophobicity_scale = {
   'ILE':4.5,
   'VAL':4.2,
   'LEU':3.8,
   'PHE':2.8,
   'CYS':2.5,'CYX':2.5,
   'MET':1.9,
   'ALA':1.8,
   'GLY':-0.4,
   'THR':-0.7,
   'TRP':-0.9,
   'SER':-0.8,
   'TYR':-1.3,
   'PRO':-1.6,'HYP':-1.6,
   'HIS':-3.2,'HIP':-3.2,'HID':-3.2,'HIE':-3.2,
   'GLU':-3.5,'GLH':-3.5,
   'GLN':-3.5,
   'ASP':-3.5,'ASH':-3.5,
   'ASN':-3.5,
   'LYS':-3.9,
   'ARG':-4.5,
}

charge_scale = {
   'ARG':1,
   'LYS':1,
   'HIS':0.1,
   'HIE':0.1,
   'HID':0.1,
   'GLU':-1,
   'ASP':-1
}

SASA_scale = {
   'ILE':1.850,
   'VAL':1.537,
   'LEU':1.831,
   'PHE':2.007,
   'CYS':1.404,'CYX':1.404,
   'MET':2.001,
   'ALA':1.102,
   'GLY':0.787,
   'THR':1.387,
   'TRP':2.405,
   'SER':1.172,
   'TYR':2.137,
   'PRO':1.419,'HYP':1.419,
   'HIS':1.819,'HIP':1.819,'HID':1.819,'HIE':1.819,
   'GLU':1.747,'GLH':1.747,
   'GLN':1.786,
   'ASP':1.441,'ASH':1.441,
   'ASN':1.464,
   'LYS':2.057,
   'ARG':2.290,
}

GeorgeDSASA_scale = {
   'ILE':1.850,
   'VAL':1.645,
   'LEU':1.931,
   'PHE':2.228,
   'CYS':1.461,'CYX':1.461,
   'MET':2.034,
   'ALA':1.118,
   'GLY':0.881,
   'THR':1.525,
   'TRP':2.663,
   'SER':1.298,
   'TYR':2.368,
   'PRO':1.468,'HYP':1.468,
   'HIS':2.025,'HIP':2.025,'HID':2.025,'HIE':2.025,
   'GLU':1.862,'GLH':1.862,
   'GLN':1.932,
   'ASP':1.587,'ASH':1.587,
   'ASN':1.655,
   'LYS':2.258,
   'ARG':2.560,
}


hydrophobic_set={'VAL','ILE','LEU','PHE','MET','TYR','TRP','CYS'}
Charg_set={'ARG','LYS','ASP','GLU'}
Aroma_set={'PHE','TYR','TRP'}