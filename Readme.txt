README for acf2strain Matlab program. 

This program calculates the best-fit deformation parameters from a 2-D image or
a set of 3-D images. For given parameter values set in the Sample structure
(see runsamples.m), the program will load the set of images, calculate the
autocorrelation function (ACF) for the sample, and run a least squares fit for
the deformation parameters that return the ACF to isotropy. The program also
calculates estimated error parameters using the method of Hext (1963) and
outputs both the results and a number of useful figures. For details of the
method see: Thissen, C. J., & Brandon, M. T. (2015). An Autocorrelation Method
for Three-Dimensional Strain Analysis. Journal of Structural Geology.

For comments, questions, or suggestions, please email cthissen@gmail.com or 
leave a comment under the issues tab at github.com/cthissen/acf2strain

Christopher J. Thissen, Yale University 
Mark T. Brandon, Yale University
Time-stamp: <Oct 14 2015>

Version history: 
This version of the code is included as supplemental information associated with
the following publication: Thissen, C. J., & Brandon, M. T. (2015). An 
Autocorrelation Method for Three-Dimensional Strain Analysis. Journal of 
Structural Geology. For the most up to date version, 
visit github.com/cthissen/acf2strain or search the mathworks website.

v. 1.0 10/14/15: First release, Thissen and Brandon (Jour. of Structural 
				      Geology, 2015)
v. 0.4 08/07/15: Updated file structure to use external functions
v. 0.3 07/08/15: Added 2-component model for grain size estimate
v. 0.2 05/19/15: removed weighting by lag
v. 0.1 05/12/15: Added self-correlation grain size estimates
v. 0.0 12/20/14: alpha version

REQUIREMENTS: This program can be RAM intensive due to the need to store many
arrays of equivalent size to the three-dimensional tomogram. I have found that
16GB is sufficient for the examples included here. However, the code can be
modified to run on systems with less RAM. Please contact the author for
details.

This function requires the following toolboxes and additional codes: 
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

INSTALLATION: 
No installation is required.

RUNNING The CODE: 
We have included three examples: 1) the undeformed ooids, 2) the 3D phantom, 3)
tomogram of sample 92810-3 from the Olympic Mountains. Each sample has a script
with the prefix "run_" that can be used to run the strain analysis. For example,
to run the strain analysis on the 3D phantom, open the "run_3Dphantom.m" file 
and push the "run" button, which looks like a green play button.

The required computation time for these examples will vary from computer to
computer, but should be on the order of about 10 minutes.


OUTPUT: 
The program will make a Results folder to store the output from the code, 
including best-fit deviatoric strain results and goodness of fit parameters.
A number of useful images are plotted and saved, such as those found in 
Thissen and Brandon (2015). These can be found in the Results folder. 
If Flag.saveMatFile is set to true, the program will also save a.mat file in the
Results folder. If present, this file will be loaded to avoid running the 
non-linear fitting routine again. 



