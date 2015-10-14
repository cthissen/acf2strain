# acf2strain #

What is it?
----------------- 
This program calculates the best-fit anisotropy parameters, and estimates the 
parameter errors, from a 2-D image or a set of 3-D images (e.g. tomography). 
For details of the method see: Thissen, C. J., & Brandon, M. T. (2015). 
An Autocorrelation Method for Three-Dimensional Strain Analysis. Journal of 
Structural Geology.

For comments, questions, or suggestions, please email cthissen@gmail.com or 
leave a comment under the issues tab at github.com/cthissen/acf2strain

Christopher J. Thissen, Yale University  
Mark T. Brandon, Yale University  

 
Motivation
------------------ 
Photomicrographs and X-ray tomography provide image data of fabrics. This 
program finds the parameters that best describe the anisotropy of these images. 
One example application is to provide quantitative estimates of deformation in
geologic samples. The figure below gives one example from the paper, where the 
anisotropy parameters are the maximum, intermediate, and minimum stretch directions, 
X, Y and Z, and their magnitudes, Sx, Sy and Sz.

![Example 92810-3](https://github.com/cthissen/acf2strain/blob/master/Example.png)


Requirements
------------------ 
acf2strain requires Matlab 2014 or later. 

This program requires the following toolboxes and additional codes: 
1. Curve Fitting Toolbox. 
2. Image Processing Toolbox. 
3. splinefit.m is included with the code, and is also available online at:
   http://www.mathworks.com/matlabcentral/fileexchange/13812-splinefit
4. sort_nat is included with the code, and is also available online at:
   http://www.mathworks.com/matlabcentral/fileexchange/10959-sort-nat--natural-
    order-sort
5. export_fig is a useful suite of programs for saving publication-quality 
   figures, and is included with the code. It is also available online at: 
   http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
6. coolwarm is included with the code. 
7. cmapscale is included with the code.

This program can be RAM intensive due to the need to store many
arrays of equivalent size to the three-dimensional tomogram. I have found that
16GB is sufficient for the examples included here. However, the code can be
modified to run on systems with less RAM. Please contact the author for
details.

Installation
------------------ 
No installation is necessary.

Usage
------------------ 
We have included three examples: 1) the undeformed ooids, 2) the 3D phantom, 3)
tomogram of sample 92810-3 from the Olympic Mountains. Each sample has a script
with the prefix "run_" that can be used to run the strain analysis. For example,
to run the strain analysis on the 3D phantom, open the "run_3Dphantom.m" file 
and push the "run" button, which looks like a green play button.

The required computation time for these examples will vary from computer to
computer, but should be on the order of about 10 minutes.

Output
------------------ 
Output is stored in the Results folder, including best-fit deviatoric strain 
results and goodness of fit parameters. A number of useful images are plotted
and saved, such as those found in Thissen and Brandon (2015). 
If Flag.saveMatFile is set to true, the program will also save a.mat file in the
Results folder. If present, this file will be loaded to avoid running the 
non-linear fitting routine again. 

The Latest Version
------------------ 
Details of the latest version can be found on the github project page under 
  server project page under https://github.com/cthissen/acf2strain
  
Contributors
------------------ 
Christopher Thissen, Yale University. christopher.thissen@yale.edu  
Mark Brandon, Yale University. mark.brandon@yale.edu

Feedback
------------------ 
Your comments are welcome! If you find any bugs or have feature requests report them to
Christopher Thissen, christopher.thissen@yale.edu. 

Issues can also be reported online: https://github.com/cthissen/Drex-MATLAB/issues


License
------------------ 
See License file.