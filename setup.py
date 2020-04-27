import setuptools

#with open("README.md", "r") as fh:
#    long_description = fh.read()

setuptools.setup(
    name="PyFEM-TUE",
    version="3.0.2",
    author="J.J.C. Remmers",
    author_email="j.j.c.remmers@tue.nl",
    description="A Python Finite Element package",
    keywords = ['Finite Element Method', 'Mechanics', 'Wiley'],
    long_description="PyFEM is a python-based finite element code accompanies the book: 'Non-Linear Finite Element Analysis of Solids and Structures' by R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel. The code is open source and intended for educational and scientific purposes only. If you use PyFEM in your research, the  developers would be grateful if you could cite the book in your work.",
    long_description_content_type="text/markdown",
    url="https://github.com/jjcremmers/PyFEM",
    packages=setuptools.find_packages(),
    install_requires=['numpy','scipy','matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
