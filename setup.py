from distutils.command.sdist import sdist as sdist_orig
from distutils.errors import DistutilsExecError
from setuptools import setup, find_packages
from os import path

PKG_NAME = MOD_NAME = "LAMPgRNAtor"

DESCRIPTION = """ 
LAMPgRNAtor
RfxCas13d guide RNA Scoring and Selection
Based on Guide Efficacy, Shannon Entropy and Conservation Scores of gRNA 
"""

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as fh, open(path.join(here, "requirements.txt")) as req:
    install_requires = [pkg.strip() for pkg in req]

__version__ = ""
exec(open("{}/_version.py".format(MOD_NAME)).read())

data = [
    "misc/DeepCpf1_weights.h5",
    "misc/Seq_deepCpf1_weights.h5",
]

setup(
    name=PKG_NAME,
    version=__version__,
    author="Muhammad Irfan",
    author_email="muhammad_irfan@gis.astar.edu.sg",
    description="LAMP and Cas12 guide RNA generator",
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/chewlabSB2",
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    package_dir={'LAMPgRNAtor': 'LAMPgRNAtor'}, 
    package_data={'LAMPgRNAtor': data},
    entry_points={
        "console_scripts": [
            #"LAMPgRNAtor={}.LAMPgRNAtor:main".format(MOD_NAME),
            "LAMPgRNAvalidate={}.LAMPgRNAvalidate:main".format(MOD_NAME),
        ],
    },
    
    install_requires=install_requires,
    include_package_data=True,
    python_requires=">=3.5",
)
