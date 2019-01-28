# unconstrained_face

Python 3.5 port for the original implentation of Roth et al 's photo-stereo approach for face reconstruction ( http://cvlab.cse.msu.edu/project-face-recon.html )
Fixes and changes for py35, makes it much easier to run in windows w/ conda, since the code is provided as an executable

Readme, for running with Anaconda and python 3.5

MATLAB Compiler

1. Prerequisites for Deployment 
. Python with SciPy and dlib.
  The easiest installation is to use Anaconda Python for windows: https://www.continuum.io/downloads


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
