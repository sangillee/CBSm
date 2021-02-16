# CBSm 1.0.2

MATLAB package for using Cubic Bezier Spline as a function approximator in modeling latent utility functions.  
While Cubic Bezier Splines (CBS) are heavily used in the graphics software industry, it can also be used as a flexible  
function approximation tool given the right constraints. The CBS package provides a method to calculate the y value  
from a x value given an appropriately constrained CBS curve. It then uses this method to approximate latent utility  
functions in intertemporal choice and risky choice data.  

The folder 'CBSm' contains the actual functions necessary to run. This is the only folder that's technically needed.  
Add this folder to the MATLAB path, and you should be good to go.

The folder 'examples' contains example MATLAB script and data to show how to use the functions in the folder 'CBSm'.  
This is for your aid only, and is not a necessary part of the package.

The folder 'java_src' contains the original java code for the function 'CBScalc.class', which is inside CBSm.  
This is here so that you can look at the source code if you want to. But this is not a necessary part of the package  
as the compiled code 'CBScalc.class' is already located inside the folder 'CBSm'.

The last stable build version can be found under the release tab.

Fun demo video: 
https://www.youtube.com/watch?v=obR1dpddYow&ab_channel=ArthurSangilLee

Citation:
Lee, S., Glaze, C. M., Bradlow, E. T., & Kable, J. W. (2020). Flexible Utility Function Approximation via Cubic Bezier Splines. psychometrika, 85(3), 716-737. https://doi.org/10.1007/s11336-020-09723-4