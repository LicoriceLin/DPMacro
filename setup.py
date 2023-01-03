from setuptools import setup,find_packages
#please use the environment.yaml to build env!
#key component is annotated inside the install_requires
install_requires = [
    # "biopython==1.79", "pandas==1.4.0", "biopandas==0.2.9",
    # "biopandas", "rdkit", "anarci==2021.02.04","rdkit","clustalw","dssp==3.0.0"

]
setup(
    name='DPMacro',
    version='0.1',
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
