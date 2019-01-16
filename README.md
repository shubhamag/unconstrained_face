# unconstrained_face
dataset prep code in python for being able to use Roth et al 's photo-stereo approach for face reconstruction
From : http://cvlab.cse.msu.edu/project-face-recon.html

MATLAB Compiler

1. Prerequisites for Deployment 
. Python with SciPy and dlib.
  The easiest installation is to use Anaconda Python for windows: https://www.continuum.io/downloads
  To note that, we test the code under Python 2.7, Anaconda 4.1. Any version below Anaconda 4.1 may lead to a failure of the code.
  This includes all of the dependencies besides dlib.

  The easiest install of dlib is using the python package manager
  >>pip install dlib

. Verify the MATLAB Runtime is installed and ensure you    
  have installed version 9.0 (R2015b).   

. If the MATLAB Runtime is not installed, do the following:
  (1) enter
  
      >>mcrinstaller
      
      at MATLAB prompt. The MCRINSTALLER command displays the 
      location of the MATLAB Runtime installer.

  (2) run the MATLAB Runtime installer.

Or download the Windows 64-bit version of the MATLAB Runtime for R2015b 
from the MathWorks Web site by navigating to

   http://www.mathworks.com/products/compiler/mcr/index.html
   
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
Package and Distribute in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.    


NOTE: You will need administrator rights to run MCRInstaller. 

2. Preprocess photo collection

Your photo collection should be a folder containing *.jpg of the same person.  An example collection is provided.

To preprocess, run the python script from this folder with the path to the folder.

>>python Python\detect_landmarks.py Collections\LDiCaprio

The script will detect faces, and perform landmark alignment.  The processed faces will be placed inside a folder title 'Export' within the collection folder.
