/**@author Sangil Arthur Lee
 * CBS class for static methods regarding calculation of y given x. 11/13/2018
 * modified: 1/27/2020: 1. moving all error checking to outer-shells and not here.
 *                      2. combining two methods into one for faster speed.
 *                      3. checking for monotonicity constraints and minimum length
 * Use this directly from matlab or R (or others) to calculate coordinates.
 */

import java.lang.Math;

public class CBScalc {
    public static double[] getyhat(double[] xpos, double[] ypos, double[] x){
        int numpoints = xpos.length;
        int numpiece = (numpoints-1)/3;
        int exitflag = 0; // an exitflag to make sure that the input data x is within the CBS range
        double[] y = new double[x.length]; // array for returning coordinates
        double[][] abcd = new double[3][numpiece]; // array for storing cubic equation coefficients for a given piece of CBS
        double d; // storage space for cubic equation coefficient
        int i,j,piecenum; // storage spaces for indice
        double t = 0, ft; // storage spaces for t and f(t)
        double precision = 0.0000001; // parameter for solution sensitivity
        double[] bound = {0,0}, boundval = {0,0}; // storage space for t min max bounds and ft at those bounds
        boolean check1, check2, check3; // monotonicity constraint checking
        
        
        // iterate through each piece and store cubic equation coefficients
        for(i=3;i<numpoints;i+=3){
            piecenum = i/3-1;
            abcd[0][piecenum] = -xpos[i-3]+3*xpos[i-2]-3*xpos[i-1]+xpos[i];
            abcd[1][piecenum] = 3*xpos[i-3]-6*xpos[i-2]+3*xpos[i-1];
            abcd[2][piecenum] = -3*xpos[i-3]+3*xpos[i-2];
            if((xpos[i]-xpos[i-3])< precision) {
                if((xpos[i]-xpos[i-3])< 0) {throw new RuntimeException("A CBS anchor point is placed before its previous anchor point. Possibly order of xpos in CBS is reversed.");}
                else {throw new RuntimeException("one or more CBS piece is too short for stable computation. If you repeatedly see this message, consider reducing the number of pieces");}
            }
            check1 = -Math.sqrt((xpos[i]-xpos[i-1])*(xpos[i-2]-xpos[i-3])) < (xpos[i-1]-xpos[i-2]);
            check2 = xpos[i-3] <= xpos[i-2];
            check3 = xpos[i-1] <= xpos[i];
            if(!check1 || !check2 || !check3){
                throw new RuntimeException("X coordinates not monotonic as a function of t. Multiple y-values may exist for x");
            }
        }
        
        // iterate through all data points and calculate their y-coordinates
        for(i=0;i<x.length;i++){
            exitflag = 1;
            for(j=3;j<numpoints;j+=3){ // looping through anchor positions skipping the first one
                if(xpos[j-3]<=x[i] && x[i]<=xpos[j]){
                    piecenum = j/3-1;
                    // root-finding algorithm
                    d = xpos[j-3]-x[i];
                    bound[0] = 0; bound[1] = 1; // min-max bound of t
                    boundval[0] = d; boundval[1] = abcd[0][piecenum]+abcd[1][piecenum]+abcd[2][piecenum]+d; // value of cubic equation at min-max bound
                    if(Math.abs(boundval[0]) < precision){t = 0;}
                    else if(Math.abs(boundval[1]) < precision){t = 1;}
                    else{
                        t = (-boundval[0])/(boundval[1]-boundval[0]); // linear initial approximation
                        ft = abcd[0][piecenum]*t*t*t + abcd[1][piecenum]*t*t + abcd[2][piecenum]*t + d; // initial function value at initial t
                        while (Math.abs(ft) > precision){ // loop until function value gets very close to zero
                            // update bounds
                            if((ft*boundval[0])>0){bound[0] = t;boundval[0] = ft;}
                            else{bound[1] = t;boundval[1] = ft;}
                            // newton-raphson update
                            t = t-ft/(3*abcd[0][piecenum]*t*t + 2*abcd[1][piecenum]*t + abcd[2][piecenum]);
                            if(bound[0]>t || t>bound[1]){t = (bound[0]+bound[1])/2;} // if newton-raphson is out-of-bound, use bisection
                            ft = abcd[0][piecenum]*t*t*t + abcd[1][piecenum]*t*t + abcd[2][piecenum]*t + d; // update function value at new t
                        }
                    }
                    exitflag = 0;
                    break;
                }
            }
            if(exitflag==1) {throw new RuntimeException("input to CBS function is out-of-bounds covered by CBS");}
            //conversion from t to y using Casteljau's theorem
            y[i] = (1-t)*((1-t)*((1-t)*ypos[j-3]+t*ypos[j-2])+t*((1-t)*ypos[j-2]+t*ypos[j-1]))+t*((1-t)*((1-t)*ypos[j-2]+t*ypos[j-1])+t*((1-t)*ypos[j-1]+t*ypos[j]));
        }
        return y;
    }
}