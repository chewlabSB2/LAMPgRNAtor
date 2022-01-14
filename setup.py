#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys

try:
    from setuptools import setup, Extension, find_packages
    from setuptools.command import install_lib, sdist, build_clib
except ImportError:
    from distutils.core import setup, Extension, find_packages
    from distutils.command import install_lib, sdist, build_clib

from distutils import log as setup_log

from os.path import join as pjoin
from os.path import relpath as rpath
from os.path import abspath 

PKG_NAME = MOD_NAME = "LAMPgRNAtor"

DESCRIPTION = """ 
LAMPgRNAtor
LAMP Primers and Cas12 guide RNA generator
Based on Cas12 Guide RNA Efficacy, Shannon Entropy and Conservation Scores of LAMP Primers and gRNAs 
"""

PACKAGE_PATH = abspath(os.path.dirname(__file__))

with open(pjoin(PACKAGE_PATH, "README.md")) as fh, open(pjoin(PACKAGE_PATH, "requirements.txt")) as req:
    install_requires = [pkg.strip() for pkg in req]

__version__ = ""
exec(open("{}/_version.py".format(MOD_NAME)).read())

data = [
    "misc/DeepCpf1_weights.h5",
    "misc/Seq_deepCpf1_weights.h5",
]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Platform-dependent binary for `make` commands
MAKE_BIN = 'mingw32-make' if sys.platform == 'win32' else 'make'
# Indicates whether the Primer3 library has been built during install
P3_BUILT = False

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Package paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
MODULE_PATH = pjoin(PACKAGE_PATH, 'LAMPgRNAtor')
LIBPRIMER3_PATH = pjoin(MODULE_PATH, 'src')
THERMO_PARAMS_PATH = pjoin(LIBPRIMER3_PATH, 'Primer3_Par')

LIBPRIMER3_FPS = [
    rpath(pjoin(LIBPRIMER3_PATH, 'common.c')),
]

LIBPRIMER3_THERMO_FPS = [
    rpath(pjoin(root, fp), MODULE_PATH) for root, _, fps in
    os.walk(THERMO_PARAMS_PATH) for fp in fps
]   

PACKAGE_FPS = (
    LIBPRIMER3_THERMO_FPS +
    ['LAMPgRNAtor/GLAPD.pyx']
)

# ~~~~~~~~~~~~~~~~~~~~~ Primer3 C library build helpers ~~~~~~~~~~~~~~~~~~~~~ #

def p3Clean():
    '''
    Run `make clean` for libprimer3 in a platform-dependent manner
    '''
    proc = subprocess.Popen(
        '{} clean'.format(MAKE_BIN),
        shell=True,
        cwd=LIBPRIMER3_PATH,
    )
    proc.wait()


def p3Build():
    '''
    Run `make clean && make` for libprimer3 in a platform-dependent manner
    '''
    proc = subprocess.Popen(
        '{} clean; {}'.format(MAKE_BIN, MAKE_BIN),
        shell=True,
        cwd=LIBPRIMER3_PATH,
    )
    proc.wait()


# Insure that the copied binaries are executable
def makeExecutable(fp):
    '''
    Adds the executable bit to the file at filepath `fp`
    '''
    mode = ((os.stat(fp).st_mode) | 0o555) & 0o7777
    setup_log.info("Adding executable bit to %s (mode is now %o)", fp, mode)
    os.chmod(fp, mode)

LIBPRIMER3_BINARIES = []
LIBPRIMER3_BINARY_FPS = []
class CustomInstallLib(install_lib.install_lib):
    '''
    Custom library installer to ensure that libprimer3 binaries are installed
    and made executable.
    Installed and invoked internally by `setuptools`/`distutils`
    '''
    def run(self):
        global P3_BUILT
        super().run()
        if not P3_BUILT:
            p3Clean()
            p3Build()
            P3_BUILT = True
        binary_dest_fps = [
            pjoin(
                self.install_dir,
                'src',
                fn,
            ) for fn in LIBPRIMER3_BINARIES
        ]
        for src_fp, dest_fp in zip(LIBPRIMER3_BINARY_FPS, binary_dest_fps):
            shutil.copyfile(src_fp, dest_fp)
            makeExecutable(dest_fp)


class CustomSdist(sdist.sdist):
    '''
    Custom sdist packager, ensures that libprimer3 build artifacts are removed
    prior to packaging.
    Installed and invoked internally by `setuptools`/`distutils`
    '''
    def run(self):
        global P3_BUILT
        # Clean up the primer3 build prior to sdist command to remove
        # binaries and object/library files
        p3Clean()
        P3_BUILT = False
        # Remove the build Cython
        os.remove(os.path.join(MODULE_PATH, 'GLAPD.c'))
        super().run()


class CustomBuildClib(build_clib.build_clib):
    '''
    Custom C library builder, ensures that libprimer3 is built prior to C
    library builds.
    Installed and invoked internally by `setuptools`/`distutils`
    '''
    def run(self):
        global P3_BUILT
        # Build primer3 prior to building the extension, if not already built
        if not P3_BUILT:
            p3Clean()
            p3Build()
            P3_BUILT = True
        super().run()


# ~~~~~~~~~~~~~~~~~~~~~~ Cython / C API extension setup ~~~~~~~~~~~~~~~~~~~~~ #

if sys.platform == 'win32':
    EXTRA_COMPILE_ARGS = ['']
else:
    EXTRA_COMPILE_ARGS = [
        '-Wno-error=declaration-after-statement',
        '-Wno-unused-function',
    ]

GLAPD_ext = Extension(
    name = 'GLAPD',
    sources=[pjoin('LAMPgRNAtor', 'GLAPD.pyx')] + LIBPRIMER3_FPS,
    include_dirs=[LIBPRIMER3_PATH],
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Ensure that libprimer3 is built prior to the `build_ext` or `install` cmds
if ('build_ext' in sys.argv or 'install' in sys.argv):
    if not P3_BUILT:
        p3Clean()
        p3Build()
        P3_BUILT = True

setup(
    name=PKG_NAME,
    version=__version__,
    license='GPLv2',
    author="Muhammad Irfan",
    author_email="muhammad_irfan@gis.astar.edu.sg",
    description="LAMP and Cas12 guide RNA generator",
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/chewlabSB2/LAMPgRNAtor",
    classifiers=[
        'Programming Language :: C',
        'Programming Language :: Cython',
        #'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
    ],
    packages= find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    ext_modules=[GLAPD_ext],
    package_dir={'LAMPgRNAtor': 'LAMPgRNAtor'}, 
    package_data={'LAMPgRNAtor': data + PACKAGE_FPS},
    
    cmdclass={
        'install_lib': CustomInstallLib,
        'sdist': CustomSdist,
        'build_clib': CustomBuildClib,
    },
    entry_points={
        "console_scripts": [
            "LAMPgRNAtor-main={}.LAMPgRNAtor:main".format(MOD_NAME),
            "LAMPgRNAtor-LAMP={}.LAMP:main".format(MOD_NAME),
            #"LAMPgRNAtor-Cas12={}.Cas12gRNAtor:main".format(MOD_NAME),
            #"LAMPgRNAtor-Cas13={}.Cas13gRNAtor:main".format(MOD_NAME),
        ],
    },
    setup_requires=['Cython', 'setuptools>=18.0'],
    zip_safe=False,
    install_requires=install_requires,
    include_package_data=True,
    python_requires=">=3.5",
)
