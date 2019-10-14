import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyFEM",
    version="0.0.1",
    author="J.J.C. Remmers",
    author_email="author@example.com",
    description="A python FEM package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jjcremmers/PyFEM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
