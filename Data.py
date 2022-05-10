'''
'''
from Bio.Data.IUPACData import protein_letters_3to1,protein_letters_1to3


amino1to3dict = protein_letters_1to3

amino3to1dict = {
    **protein_letters_3to1,
    **{'ASH': 'A','CYX': 'C','HYP': 'P',
        'HID': 'H','HIE': 'H','HIP':'H','MSE': 'M'}
}

CDR_annotations = {'imgt':{
                        'FR1':1,
                        'CDR1':27,
                        'FR2':39,
                        'CDR2':56,
                        'FR3':66,
                        'CDR3':105,
                        'FR4':118,
                        'end':130
                    },
                    'i':{
                        'FR1':1,
                        'CDR1':27,
                        'FR2':39,
                        'CDR2':56,
                        'FR3':66,
                        'CDR3':105,
                        'FR4':118,
                        'end':130
                    },
                    'a':{
                        'FR1':1,
                        'CDR1':25,
                        'FR2':41,
                        'CDR2':58,
                        'FR3':78,
                        'CDR3':109,
                        'FR4':138,
                        'end':149 
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
                    }
             }

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


# GeorgeDSASA = {
#    'I':1.850,
#    'V':1.645,
#    'L':1.931,
#    'F':2.228,
#    'C':1.461,
#    'M':2.034,
#    'A':1.118,
#    'G':0.881,
#    'T':1.525,
#    'W':2.663,
#    'S':1.298,
#    'Y':2.368,
#    'P':1.468,
#    'H':2.025,
#    'E':1.862,
#    'Q':1.932,
#    'D':1.587,
#    'N':1.655,
#    'K':2.258,
#    'R':2.560,
# }