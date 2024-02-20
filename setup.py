from setuptools import setup, find_packages

with open("src/rrikindp/__init__.py") as f:
    for line in f.readlines():
        if line.startswith('__version__'):
            VERSION = line.strip().split()[-1][1:-1]

#with open('requirements.txt') as f:
#    requirements = f.read().splitlines()

setup(
    name='rrikindp',
    version=VERSION,
    packages=find_packages('src'),
    package_dir={'':'src'},
    description='Evaluation of thermodynamic and kinetic features of RNA-RNA interactions',
    author='Maria Waldl',
    author_email='code@waldl.org',
    url='https://github.com/mwaldl/RRIkinDP',
#    install_requires=requirements,
    python_requires='>=3.6'
)
