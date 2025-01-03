from setuptools import setup, find_packages

setup(
    name='Proteon',
    version='0.1',
    description='A package for protein sequence and structure analysis and manipulation',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Mohammad Assadizadeh',
    author_email='mohammadassadizadeh@gmail.com',
    url='https://github.com/mohamadasadizade/Proteon', 
    packages=find_packages(),  
    install_requires=[
        'numpy',     
        'mdtraj',     
        'matplotlib',  
        'biopython',   
    ],
    python_requires='>=3.6'
)
