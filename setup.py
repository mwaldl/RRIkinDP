import os.path
import sys

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Package metadata is declared in pyproject.toml.

# Prefixes searched for the headers and libraries of IntaRNA, the Vienna RNA
# package and boost. sys.base_prefix is included so that they are still found
# when the build runs inside a virtual environment, as it does for isolated
# builds.
PREFIXES = list(dict.fromkeys([sys.prefix, sys.base_prefix]))


def get_version():
    """Return __version__ of the package.

    The package can not be imported to get its version, as importing it requires
    the extension that is built by this script. A leading 'v' of the release tag
    is dropped to get a PEP 440 conform version.
    """
    version_file = os.path.join("src", "rrikindp", "__init__.py")
    with open(version_file) as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("=")[-1].strip().strip("\"'").lstrip("v")
    raise RuntimeError(f"no __version__ found in {version_file}")


ext_modules = [
    Pybind11Extension(
        "libRRIkinDP", ["src/rrikindp/libRRIkinDP.cpp"], # RRIkinDP.cpp is included by libRRIkinDP.cpp
        include_dirs=["src/rrikindp"] + [os.path.join(prefix, "include") for prefix in PREFIXES],
        library_dirs=[os.path.join(prefix, "lib") for prefix in PREFIXES],
        libraries=[
            "boost_program_options",
            "boost_filesystem",
            "boost_system",
            "boost_regex",
            "IntaRNA",
            "RNA",
            "easylogging",
        ],
        cxx_std=17,
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
    )
]

setup(
    version=get_version(),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
