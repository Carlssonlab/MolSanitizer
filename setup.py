from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='MolSanitizer',
    version='1.0',
    description='Code to clean SMILES for ML and Docking (protonation, tautomers,...)',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
    'Development Status :: Alpha',
    'License :: MIT License',
    'Programming Language :: Python :: 3.11',
    "Operating System :: OS Linux",
    'Topic :: Molecular Mechanics :: Docking :: AI :: ML',
    ],
    keywords='docking, drug design',
    url='https://github.com/Isra3l/MolSanitizer.git',
    author='Israel Cabeza de Vaca Lopez',
    author_email='israel.cabezadevaca@icm.uu.se',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'markdown', 'pandas', 'numpy'
    ],
    entry_points={
        'console_scripts': ['msani=MolSanitizer.molSanitizer:main'],
    },
    include_package_data=True,
    python_requires='>3.9',
    zip_safe=False)