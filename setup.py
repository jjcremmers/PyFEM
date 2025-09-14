# setup.py
from pathlib import Path
from setuptools import setup, find_packages

README = Path(__file__).with_name("README.md")
long_description = README.read_text(encoding="utf-8") if README.exists() else ""

setup(
    name="pyfem",
    version="0.1.0",
    description="A Python finite element code (educational)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jjcremmers/PyFEM",
    author="J.J.C. Remmers",
    license="GPL-3.0-only",
    packages=find_packages(include=["pyfem", "pyfem.*"]),
    include_package_data=True,  # works with MANIFEST.in if you add one
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "meshio",
        "h5py",
        "PySide6",
        "vtk",
        #"pyyaml",   # if you load YAML props in api.py
        # add other runtime deps used by pyfem/*
    ],
    #extras_require={
    #    "dev": [
    #        "pytest",
    #        "black",
    #        "mypy",
    #    ],
    #},
    entry_points={
        "console_scripts": [
            "pyfem2=pyfem.cli:main",
            "pyfem-gui=pyfem.gui.app:main"
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    project_urls={
        "Homepage": "https://github.com/jjcremmers/PyFEM",
        "Issues": "https://github.com/jjcremmers/PyFEM/issues",
    },
    zip_safe=False,
)

