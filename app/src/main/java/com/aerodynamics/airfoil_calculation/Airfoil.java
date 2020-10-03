package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
// This class defines an Airfoil type object.                            					      //
////////////////////////////////////////////////////////////////////////////////////////////////////

public class Airfoil {

	public int id;
	public String airfoil_str;
	public double chord;
	public double theta;
	public double xcRefPoint;
	public double[] refPoint = new double[2];
	public int nChordPanels;

	public Airfoil(int id, String airfoil_str, double chord, double theta, double xcRefPoint,
				   double[] refPoint, int nChordPanels) {
		this.id = id;
		this.airfoil_str = airfoil_str;
		this.chord = chord;
		this.theta = theta;
		this.xcRefPoint = xcRefPoint;
		this.refPoint[0] = refPoint[0];
		this.refPoint[1] = refPoint[1];
		this.nChordPanels = nChordPanels;
	}
}
