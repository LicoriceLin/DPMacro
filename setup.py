from setuptools import setup,find_packages
#please use the environment.yaml to build env!
#key component is annotated inside the install_requires
install_requires = [
    "biopython>=1.79", "pandas", "biopandas", 
    # higher version of biopython doesn't contains ShrakeRupley module!
    
    ## optional, only work in Fv modules
    # "anarci==2021.02.04",
    
    ## optional, only work in AFill.
    # "pymol-open-source" #missing will lead to ImportError
    # "clustalw","dssp==3.0.0" # missing will lead to Error when evoke the binaries. 
    
    ## maybe useful?
    # "rdkit", 
]
setup(
    name='DPMacro',
    version='0.3',
    author='LicoriceLin',
    author_email='dengzf@dptech.net',
    description=('PDB structure manipulation and feature engineering.'),
    url='https://git.dp.tech/macromolecule/DPMacro',
    license=None,
    keywords='PDB processing & embedding',
    install_requires=install_requires,
    packages=find_packages(),
    zip_safe=False,
    #packages=packages,
    include_package_data=True)
