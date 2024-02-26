from setuptools import setup, find_packages, Extension
import sys
import os.path
import pybind11

with open("src/rrikindp/__init__.py") as f:
    for line in f.readlines():
        if line.startswith('__version__'):
            VERSION = line.strip().split()[-1][1:-1]

#with open('requirements.txt') as f:
#    requirements = f.read().splitlines()

print(sys.prefix)

setup(
    name='rrikindp',
    version=VERSION,
    ext_modules=[
        Extension(
            'libRRIkinDP', ['src/rrikindp/libRRIkinDP.cpp'],
            include_dirs=['src/rrikindp',  os.path.join(sys.prefix, 'include')],  # Include directories for header files
            libraries=['boost_program_options' , 'boost_filesystem' ,'boost_system', 'boost_regex','IntaRNA', 'RNA', 'easylogging', 'omp'],         # Libraries to link against
            library_dirs=['src/rrikindp', pybind11.get_include()],      # Library directories
            extra_compile_args=['-std=c++14', '-fopenmp', '-fpermissive'],  # Additional compiler arguments
            language='c++'
            )
        ],
    packages=find_packages('src'),
    package_dir={'':'src'},
    description='Evaluation of thermodynamic and kinetic features of RNA-RNA interactions',
    author='Maria Waldl',
    author_email='code@waldl.org',
    url='https://github.com/mwaldl/RRIkinDP',
#    install_requires=requirements,
    python_requires='>=3.6'
)
