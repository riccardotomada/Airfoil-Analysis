package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//This class defines a FreeStream type object.  					                              //
////////////////////////////////////////////////////////////////////////////////////////////////////

public class FreeStream {

	public double P;
	public double rho;
	public double v;
	public double alpha;
	public double[] vvec = new double[2];
	public double dyn_visc;
	public double kin_visc;

	public FreeStream(double P, double rho, double v, double alpha) {

		this.P = P;
		this.rho = rho;
		this.v = v;
		this.alpha = alpha;

		this.vvec[0] = this.v*Math.cos(this.alpha);
		this.vvec[1] = this.v*Math.sin(this.alpha);

		this.dyn_visc = 181.0/10000000.0;
		this.kin_visc = 181.0/10000000.0/1.225;
	}

}

