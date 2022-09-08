package com.aerodynamics.airfoil_calculation;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.widget.TextView;

public class ThinAirfoil extends AppCompatActivity {

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_thin_airfoil);

		String nacaDW;
		double chLnW, staggW, aOA;
		int nPan;
		GlobalFunctions myGlobal = new GlobalFunctions();

		nacaDW   = "NACA";
		nacaDW   = nacaDW + getIntent().getStringExtra("NACA_DW");
		chLnW    = getIntent().getDoubleExtra("ChordW", 1.0);
		staggW   = getIntent().getDoubleExtra("StaggerW", 0.0);
		staggW   = Math.toRadians(staggW);
		nPan     = getIntent().getIntExtra("nPan", 20);
		aOA      = getIntent().getDoubleExtra("AngleOA", 0.0);
		aOA      = Math.toRadians(aOA);

		double[] refPoint = {0,0};

		Airfoil airfoilW = new Airfoil(1, nacaDW, chLnW, staggW,
				0.25,refPoint, nPan);

		//Thin walled theory --- Calculation of the Cl using the approximation of the thin walled
		//theory, camber problem.

		char m,p;
		char[] mm;
		double M,POS,alfa0,beta0,cl_thin,cm_thin,alfa_th;
		int coeff;

		if(airfoilW.airfoil_str.length() == 8){
			m = airfoilW.airfoil_str.charAt(4);
			p = airfoilW.airfoil_str.charAt(5);
			M = Double.parseDouble(String.valueOf(m))/100;
			POS = Double.parseDouble(String.valueOf(p));
			if(POS==0)
				POS = Math.pow(10,-7);
			POS = POS/10;

			double eta1 = Math.acos(1-2*POS);
			if(POS == Math.pow(10,-8))
				eta1 = 0.0;

			alfa0 = (M/Math.PI)*((1/Math.pow(POS,2))*((2*POS-1)*eta1+2*(1-POS)*Math.sin(eta1)-eta1/2-0.25*Math.sin(2*eta1))+(1/(Math.pow(1-POS,2)))
							*((2*POS-1)*Math.PI-Math.PI/2-(2*POS-1)*eta1-2*(1-POS)*Math.sin(eta1)+eta1/2+0.25*Math.sin(2*eta1)));
			beta0 = (2*M/Math.PI)*((1/Math.pow(POS,2))*((2*POS-1)*(eta1/2-Math.sin(2*eta1)/4)+1.0/3.0*Math.pow(Math.sin(eta1),3))+(1/Math.pow((1-POS),2))*((2*POS-1)*Math.PI/2-
					(2*POS-1)*(eta1/2-Math.sin(2*eta1)/4)-1.0/3.0*Math.pow(Math.sin(eta1),3)));
			alfa_th = (2*M/Math.PI)*(1/Math.pow(POS,2)*((POS-0.5)*eta1+0.5*Math.sin(eta1))+(1/Math.pow((1-POS),2))*((POS-0.5)*Math.PI-(POS-0.5)*eta1-0.5*Math.sin(eta1)));
			cl_thin = 2*Math.PI*(staggW+aOA-alfa0);
			cm_thin = -Math.PI/2*(beta0-alfa0);
		}
		else{
			mm = new char[]{airfoilW.airfoil_str.charAt(4), airfoilW.airfoil_str.charAt(5),
					airfoilW.airfoil_str.charAt(6)};
			coeff = Integer.parseInt(String.valueOf(mm));
			double q,K;

			switch(coeff){
				case 210:
					q = 0.0580;
					K = 361.4;
					break;
				case 220:
					q = 0.1260;
					K = 51.64;
					break;
				case 230:
					q = 0.2025;
					K = 15.957;
					break;
				case 240:
					q = 0.2900;
					K = 6.643;
					break;
				case 250:
					q = 0.3910;
					K = 3.230;
					break;
				default:
					q = 0;
					K = 0;
			}

			double a2 = 0.5*K;
			double a1 = K*(0.5-q);
			double a0 = K/6*(3.0/4.0-3*q+3*Math.pow(q,2)-Math.pow(q,3));
			double b0 = -K*Math.pow(q,3)/6;

			double[] I = new double[3];

			for(int check = 0; check<3; check++){
				I[check] = myGlobal.gauss_legendre_integration(-0.5,q-0.5,q-0.5,0.5,a2,a1,
						a0,b0,check);
			}

			alfa0 = 2/Math.PI*I[0];
			beta0 = 8/Math.PI*I[1];

			cl_thin = 2*Math.PI*(staggW+aOA-alfa0);
			cm_thin = Math.PI/2*(beta0-alfa0);
			alfa_th = 1/Math.PI*I[2];
		}

		// Print the results

		TextView textView = findViewById(R.id.textCLthinwalledValue);
		textView.setText(String.format("%.5f", cl_thin)+ " [-]");

		textView = findViewById(R.id.textCMthinwalledValue);
		textView.setText(String.format("%.5f", cm_thin)+ " [-]");

		textView = findViewById(R.id.textAlfaTHthinwalledValue);
		textView.setText(String.format("%.5f", alfa_th)+ " [-]");
	}
}