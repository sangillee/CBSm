/**@author Sangil Arthur Lee
 * CBS class for static methods regarding calculation of y given x. 11/13/2018
 * One should be able to use this directly from matlab to calculate coordinates,
 * or use it from another java class
 */

public class CBScalc {
	// methods for calculating y given x
	public static double[] getyhat(double[] xpos, double[] ypos, double[] x){
		int numpoints = xpos.length; // number of points
		if(ypos.length != numpoints) {throw new RuntimeException("The length of x and y coordinates should match");}
		if(numpoints<4) {throw new RuntimeException("Insufficient number of coorindates");}
		if((numpoints-1)%3 != 0) {throw new RuntimeException("The number of coordinates need to be 1+3n (n > 0)");}

		int numpiece = (numpoints-1)/3;
		double[] y = new double[x.length]; // array for returning coordinates
		double[][] abcd = new double[3][numpiece]; // array for storing cubic equation coefficients for a given piece of CBS
		int piecenum; // a storage space for current piece number
		double t = 0; // a storage space for intermediate variable t
		int exitflag = 0; // an exitflag to make sure that the input data x is within the CBS range

		// iterate through each piece and store cubic equation coefficients
		for(int i=3;i<numpoints;i+=3){
			piecenum = i/3-1;
			abcd[0][piecenum] = -xpos[i-3]+3*xpos[i-2]-3*xpos[i-1]+xpos[i];
			abcd[1][piecenum] = 3*xpos[i-3]-6*xpos[i-2]+3*xpos[i-1];
			abcd[2][piecenum] = -3*xpos[i-3]+3*xpos[i-2];
		}

		// iterate through all data points and calculate their y-coordinates
		for(int i=0;i<x.length;i++){
			int j;
			exitflag = 1;
			for(j=3;j<numpoints;j+=3){ // looping through anchor positions skipping the first one
				if(xpos[j-3]<=x[i] && x[i]<=xpos[j]){
					piecenum = j/3-1;
					t = NewtonT(abcd[0][piecenum],abcd[1][piecenum],abcd[2][piecenum],xpos[j-3]-x[i]);
					exitflag = 0;
					break;
				}
			}
			if(exitflag==1) {throw new RuntimeException("x is out-of-bounds covered by CBS");}
			//conversion from t to y using Casteljau's theorem
			y[i] = (1-t)*((1-t)*((1-t)*ypos[j-3]+t*ypos[j-2])+t*((1-t)*ypos[j-2]+t*ypos[j-1]))+t*((1-t)*((1-t)*ypos[j-2]+t*ypos[j-1])+t*((1-t)*ypos[j-1]+t*ypos[j]));
		}
		return y;
	}

	private static double NewtonT(double a, double b, double c, double d){
		double precision = 0.0000001; // parameter for solution sensitivity
		double[] bound = {0,1}; // min-max bound of t
		double[] boundval = {d,a+b+c+d}; // value of cubic equation at min-max bound

		// combination of Newton-Raphson and linear approximation algorithm to find root in [0,1]
		if(Math.abs(boundval[0]) < precision){return 0;}
		else if(Math.abs(boundval[1]) < precision){return 1;}
		else{
			double t = (-boundval[0])/(boundval[1]-boundval[0]); // linear initial approximation
			double ft = a*t*t*t + b*t*t + c*t + d; // initial function value at initial t
			while (Math.abs(ft) > precision){ // loop until function value gets very close to zero
				// update bounds
				if((ft*boundval[0])>0){
					bound[0] = t;
					boundval[0] = ft;
				}
				else{
					bound[1] = t;
					boundval[1] = ft;
				}
				// newton-raphson update
				t = t-ft/(3*a*t*t + 2*b*t + c);
				if(bound[0]>t || t>bound[1]){t = (bound[0]+bound[1])/2;} // if newton-raphson is out-of-bound, use bisection
				ft = a*t*t*t + b*t*t + c*t + d; // update function value at new t
			}
			return t;
		}
	}
}