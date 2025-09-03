
Windows
-------

1. Install Python

In this instruction, we use the Anaconda Python distribution. 
Download the Python 3.6 version of Anaconda from: https://www.anaconda.com/download 
You are advised to use the default settings in the installation process, and to not select 
the option to also install Visual Studio. Besides the Python 3 interpreter the Anaconda 
distribution contains tools, such as the linear algebra package NumPy, the scientific 
computing package SciPy, the plotting library MatplotLib, and the Integrated Development 
Environment (IDE) Spyder. 

2. Install Git

* Go to the official Git website https://git-scm.com/ and download and install Git by following the download link and selecting the appropriate operating system (most likely 64-bit Windows). 
* Plenty of documentation is available through the Git site, which you will have to consult while working with Git. In particular the first two chapters of the online Git reference manual https://git-scm.com/doc are highly recommended to read.Git Bash and Git GUI are now installed on your computer. Try them both with an example project you make on your own Gitlab account. You may choose the interface you prefer, however, for the rest of this instruction the Git Bash is used. A site that may be usefull: https://www.atlassian.com/git/tutorials/git-bash In the Git Bash window, you can type Git Commands that can be found in documentation/internet. Tip: there are cheat sheets  available for the commands on the internet.

3. Install PyFEM using Git Bash

* Go to the directory on your computer where you would like to save/clone the PyFEM files, preferably this directory does not contain any spaces, e.g. C:\Users\...... Right click and select ‘Git Bash here’. 
* In the bash window you can now type   Git Commands. You can clone the PyFEM repository (or required branch) by typing ‘git clone’ followed by the https link of the repository you would like to clone and press enter. This https link can be found on the Gitlab repository by clicking on the blue ‘clone’ button (first select the required branch as this is different as the standard shown master branch). The files associated with PyFEM are automatically downloaded.  
* Open the Anaconda prompt: click the windows button and search & select it. Navigate to the directory install.py  (in the PyFEM folder you just cloned with Git Bash). You can change the current directory and navigate to your folders by using ‘cdfoldername’ to open a subfolder or ‘cd ..’ to move up again. Type now ‘python install.py’. This will check whether all required packages are available and gives you the path PyFEM is installed. An error message may occur stating 
that meshi is not installed. In this case, install meshio first bytyping conda install -c conda-forge meshio in the anaconda prompt. Then run install.py again. 

Linux
-----
