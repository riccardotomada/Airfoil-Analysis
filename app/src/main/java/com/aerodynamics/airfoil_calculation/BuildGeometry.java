package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//This class contains all the Java code about the "Geometery, CP, CL, CD" screen   			      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import androidx.appcompat.app.AppCompatActivity;
import android.graphics.Color;
import android.os.Bundle;
import android.widget.TableLayout;
import android.widget.TableRow;
import android.widget.TextView;
import android.widget.Toast;
import com.google.firebase.analytics.FirebaseAnalytics;
import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.stream.DoubleStream;

public class BuildGeometry extends AppCompatActivity {

	String nacaD;
	GlobalFunctions myGlobal;
	double chLn;
	double stagg;
	int nPan;
	double pr = 0.0;
	double den = 1;
	double spd = 30.0;
	double aOA;
	private FirebaseAnalytics mFirebaseAnalytics;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_build_geometry);
		mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);


		// Receiving the values typed by the user in the MainActivity

		nacaD = "NACA";
		nacaD = nacaD + getIntent().getStringExtra("NACA_D");
		chLn = getIntent().getDoubleExtra("Chord", 1.0);
		stagg = getIntent().getDoubleExtra("Stagger", 0.0);
		nPan = getIntent().getIntExtra("nPan", 20);
		pr = getIntent().getDoubleExtra("Pressure", 101325.0);
		den = getIntent().getDoubleExtra("Density", 1.225);
		spd = getIntent().getDoubleExtra("Speed", 30);
		aOA = getIntent().getDoubleExtra("AngleOA", 0.0);

		//From here on all the code it's quite similar to the one we implemented on MATLAB.

		double[] refPoint = {0,0};

		// Now a FreeStream type object, freeStream, is initialized, using the values typed by the user
		// in the MainActivity screen. The same is done with an Airfoil type object, named airfoil.
		FreeStream freeStream = new FreeStream(pr, den, spd, Math.toRadians(aOA));

		Airfoil airfoil = new Airfoil(1, nacaD, chLn, -Math.toRadians(stagg),
				0.25,refPoint, nPan);

		// Since many functions are meant to be used many times and shared with other activities,
		// I decided to put all of them inside a new class, named GlobalFunctions. In order to use
		// them a GlobalFunctions object, myGlobal, initialized as follows, will be called every time
		// it will be necessary.

		GlobalFunctions myGlobal = new GlobalFunctions();

		//I tried to use the same variable names as we already wrote on MATLAB.

		double[][] rr = myGlobal.getCoordinates(airfoil);
		int[][] ee = myGlobal.connectivityMatrix(rr[0]);
		int nelems = Array.getLength(ee[0]);
		Elems[] elems = new Elems[nelems];

			for(int i = 0; i<nelems;i++){
				elems[i] = new Elems();
				elems[i].airfoilId = airfoil.id;
				elems[i].id = i;
				elems[i].ver1 = new double[] {rr[0][ee[0][i]-1], rr[1][ee[0][i]-1]};
				elems[i].ver2 = new double[] {rr[0][ee[1][i]-1], rr[1][ee[1][i]-1]};
				elems[i].cen = myGlobal.midpoint(elems[i].ver1, elems[i].ver2);
				elems[i].len = myGlobal.len_norm(elems[i].ver2, elems[i].ver1);
				elems[i].tver = myGlobal.tver(elems[i].ver2, elems[i].ver1, elems[i].len);
				elems[i].nver = new double[]{-elems[i].tver[1], elems[i].tver[0]};
			}

			double[][] rrc = new double[2][nelems];
			double[] len = new double[nelems];
			double[][] tvers = new double[2][nelems];
			double[][] nvers = new double[2][nelems];

			for(int i = 0; i<nelems; i++){
				rrc[0][i] = elems[i].cen[0];
				rrc[1][i] = elems[i].cen[1];
				len[i] = elems[i].len;
				tvers[0][i] = elems[i].tver[0];
				tvers[1][i] = elems[i].tver[1];
				nvers[0][i] = elems[i].nver[0];
				nvers[1][i] = elems[i].nver[1];
			}

			double[][] A = myGlobal.coefficientMatrix(elems, nPan);
			double[] b = myGlobal.coefficientArray(freeStream, elems, nPan);
			double[][] Au = myGlobal.onBodyUMatrix(elems, nPan);
			double[][] Av = myGlobal.onBodyVMatrix(elems, nPan);

			// In order to solve the linear system I discovered the math3 library. In the following
			// lines of code I've "transformed" the matrices and arrays of type double[][] and double[]
			// into objects I can operate with using the library.
			// Eventually the solution of the linear system is transformed again in a double[] type
			// variable.

			RealMatrix Acoefficients = new Array2DRowRealMatrix(A,false);

			// LU decomposition of the matrix A.
			DecompositionSolver solver = new LUDecomposition(Acoefficients).getSolver();

			RealVector bconstants = new ArrayRealVector(b,false);
			RealVector x = solver.solve(bconstants);

			RealVector u,v;

			RealMatrix AuCoeff = new Array2DRowRealMatrix(Au,false);
			RealMatrix AvCoeff = new Array2DRowRealMatrix(Av, false);

			RealMatrix uvec = new Array2DRowRealMatrix(nelems, 2);
			RealMatrix uvecTi = new Array2DRowRealMatrix(2,nelems);
			RealMatrix uvecNi = new Array2DRowRealMatrix(2,nelems);

			u = AuCoeff.operate(x);
			v = AvCoeff.operate(x);
			u = u.mapAdd(freeStream.v*Math.cos(freeStream.alpha));
			v = v.mapAdd(freeStream.v*Math.sin(freeStream.alpha));

			uvec.setColumnVector(0,u);
			uvec.setColumnVector(1,v);


			uvecTi = uvec.transpose().copy();
			uvecNi = uvec.transpose().copy();

			double[][] uvec_vTi = new double[2][nelems];
			double[][] uvec_vNi = new double[2][nelems];

			double[] vTi = new double[nelems];
			double[] vNi = new double[nelems];
			double[] vNi_check = new double[nelems];


			for (int j = 0; j < nelems; j++) {
				for (int i = 0; i< 2; i++){
					uvecTi.multiplyEntry(i, j, tvers[i][j]);
					uvecNi.multiplyEntry(i, j, nvers[i][j]);
					uvec_vTi[i][j] = uvecTi.getEntry(i, j);
					uvec_vNi[i][j] = uvecNi.getEntry(i, j);
				}
				vTi[j] = uvec_vTi[0][j] + uvec_vTi[1][j];
				vNi[j] = uvec_vNi[0][j] + uvec_vNi[1][j];
			}

			for(int i = 0; i<nelems; i++){
				if(vNi[i]<0)
					vNi_check[i] = -vNi[i];
				else
					vNi_check[i] = vNi[i];

				if(vNi_check[i] > Math.pow(10,-6))
					Toast.makeText(this, "Ho sbagliato qualcosa", Toast.LENGTH_LONG).show();
			}

			double[] cP = new double[nelems];
			for(int i = 0; i<nelems;i++){
				cP[i] = 1.0 - vTi[i]*vTi[i]/(freeStream.v*freeStream.v);
			}

			double[][] cp_len = new double[2][cP.length];

			for(int i = 0;  i<cP.length; i++){
				cp_len[0][i] = cP[i]*len[i];
				cp_len[1][i] = cP[i]*len[i];
			}

			double[][] cp_len_nvers = new double[2][cP.length];

			for(int i = 0;  i<cP.length; i++){
				cp_len_nvers[0][i] = cp_len[0][i]*nvers[0][i];
				cp_len_nvers[1][i] = cp_len[1][i]*nvers[1][i];
			}

			double cfx = -DoubleStream.of(cp_len_nvers[0]).sum();
			double cfy = -DoubleStream.of(cp_len_nvers[1]).sum();

			double cl = -Math.sin(freeStream.alpha)*cfx + Math.cos(freeStream.alpha)*cfy;

			double sum = DoubleStream.of(len).sum();
			double Gamma = x.getEntry(nelems)*sum;

			double liftKJ = -freeStream.rho*freeStream.v*Gamma;

			//CL calculus using the Kutta-Joukowski theorem.
			double cLKJ = liftKJ / (0.5*freeStream.rho*Math.pow(freeStream.v,2)*airfoil.chord);

			double[] rrc_scaledx = new double[nelems];
			for(int i = 0; i <nelems; i++){
				rrc_scaledx[i] = rrc[0][i]/chLn;
			}

			double[] rrc_scaledy = new double[nelems];
			for(int i = 0; i <nelems; i++){
				rrc_scaledy[i] = rrc[1][i]/chLn;
			}

			double[] rrc_scaledx1 = new double[nelems/2];
			double[] rrc_scaledx2 = new double[nelems/2];
			double[] rrc_scaledy1 = new double[nelems/2];
			double[] rrc_scaledy2 = new double[nelems/2];
			double[] cP1 = new double[nelems/2];
			double[] cP2 = new double[nelems/2];
			System.arraycopy(rrc_scaledx, 0, rrc_scaledx1, 0, nelems / 2);
			for(int i = nelems/2; i < nelems; i++)
				rrc_scaledx2[i - nelems / 2] = rrc_scaledx[i];
			for(int i = 0; i <nelems/2;i++)
					rrc_scaledy1[i] = rrc_scaledy[i];
			for(int i = nelems/2; i < nelems; i++)
				rrc_scaledy2[i - nelems / 2] = rrc_scaledy[i];
			for(int i = 0; i <nelems/2;i++)
				cP1[i] = cP[i];
			for(int i = nelems/2; i < nelems; i++)
				cP2[i - nelems / 2] = cP[i];

			// Unfortunately, all the graph libraries I have found don't allow to plot anything if
			// the x array is not ordered from the smallest to the greatest value.

			// In the following lines of code I reordered the x array keeping unaltered all the
			// correspondences with the other arrays (y array for the airfoil drawing and cp array
			// for the CP curve graph).

			// This means that if two elements of the x array are switched in position, the same
			// happens for the correspondent elements of the other arrays.

			// I have to admit that this is a cumbersome, complicated and probably not efficient way
			// to rearrange the arrays, but this is the only way I've been able to overcome this problem.

			boolean swapped = true;
			int j = 0;
			double temp1, temp2, temp3;

			while (swapped){
				swapped = false;
				j++;
				for(int i = 0; i<nelems/2 - j; i++){
					if(rrc_scaledx1[i]>rrc_scaledx1[i+1]){
						temp1 = rrc_scaledx1[i];
						temp2 = rrc_scaledy1[i];
						temp3 = cP1[i];
						rrc_scaledx1[i] = rrc_scaledx1[i+1];
						rrc_scaledy1[i] = rrc_scaledy1[i+1];
						cP1[i] = cP1[i+1];
						rrc_scaledx1[i+1] = temp1;
						rrc_scaledy1[i+1] = temp2;
						cP1[i+1] = temp3;
						swapped = true;
					}
					if(rrc_scaledx2[i]>rrc_scaledx2[i+1]){
						temp1 = rrc_scaledx2[i];
						temp2 = rrc_scaledy2[i];
						temp3 = cP2[i];
						rrc_scaledx2[i] = rrc_scaledx2[i+1];
						rrc_scaledy2[i] = rrc_scaledy2[i+1];
						cP2[i] = cP2[i+1];
						rrc_scaledx2[i+1] = temp1;
						rrc_scaledy2[i+1] = temp2;
						cP2[i+1] = temp3;
						swapped = true;
					}
				}
			}

			// The following 56 lines of code refer to the airfoil and CP graph construction.

			DataPoint[] values = new DataPoint[nelems/2];
			DataPoint[] values1 = new DataPoint[nelems/2];

			DataPoint[] values2 = new DataPoint[nelems/2];
			DataPoint[] values3 = new DataPoint[nelems/2];

			for(int i =0; i<nelems/2; i++){

				DataPoint value = new DataPoint(rrc_scaledx1[i], rrc_scaledy1[i]);
				DataPoint value1 = new DataPoint(rrc_scaledx2[i], rrc_scaledy2[i]);

				values[i] = value;
				values1[i] = value1;

				DataPoint cp = new DataPoint(rrc_scaledx1[i], cP1[i]);
				DataPoint cp1 = new DataPoint(rrc_scaledx2[i], cP2[i]);

				values2[i] = cp;
				values3[i] = cp1;
			}

			final GraphView graph = findViewById(R.id.graph);

			LineGraphSeries<DataPoint> series = new LineGraphSeries<DataPoint>(values);
			LineGraphSeries<DataPoint> series1 = new LineGraphSeries<DataPoint>(values1);
			LineGraphSeries<DataPoint> series2 = new LineGraphSeries<DataPoint>(values2);
			LineGraphSeries<DataPoint> series3 = new LineGraphSeries<DataPoint>(values3);
			graph.addSeries(series);
			graph.addSeries(series1);
			double max_x = 1.0;
			double max_y = Arrays.stream(rr[1]).max().getAsDouble();
			double max_cp = Arrays.stream(cP).max().getAsDouble();
			double min_cp = Arrays.stream(cP).min().getAsDouble();
			graph.getViewport().setMaxX(max_x);
			graph.getViewport().setMaxY(3*max_y);
			graph.getViewport().setMinY(-3*max_y);
			series.setDrawDataPoints(true);
			series.setDataPointsRadius(8);
			series1.setDrawDataPoints(true);
			series1.setDataPointsRadius(8);
			graph.getViewport().setXAxisBoundsManual(true);
			graph.getViewport().setYAxisBoundsManual(true);
			graph.getViewport().setScalableY(false);
			graph.setTitle("Geometry and CP");
			graph.setTitleTextSize(40);
			graph.getSecondScale().addSeries(series2);
			graph.getSecondScale().addSeries(series3);
			graph.getSecondScale().setMinY(1*min_cp);
			graph.getSecondScale().setMaxY(1*max_cp);
			series.setColor(Color.rgb(0,145,234));
			series1.setColor(Color.rgb(0,145,234));
			series2.setColor(Color.rgb(230,81,0));
			series3.setColor(Color.rgb(126,87,194));
			graph.getGridLabelRenderer().setVerticalLabelsSecondScaleColor(Color.BLACK);

			// The code below refers to 2 tables, which show the CP behaviour as a function of the position
			// respectively on the lower and upper surface of the airfoil.

			TableLayout table = (TableLayout) findViewById(R.id.table);
			TableLayout table1 = (TableLayout) findViewById(R.id.table1);

			TableRow row_header = new TableRow(this);
			TableRow row_header1 = new TableRow(this);
			TextView x_c = new TextView(this);
			x_c.setText(R.string.x_c);
			TextView x_c1 = new TextView(this);
			x_c1.setText(R.string.x_c);
			TextView cP_text = new TextView(this);
			cP_text.setText(R.string.cp);
			TextView cP_text1 = new TextView(this);
			cP_text1.setText(R.string.cp);
			row_header.addView(x_c);
			row_header.addView(cP_text);
			row_header1.addView(x_c1);
			row_header1.addView(cP_text1);
			table.addView(row_header);
			table1.addView(row_header1);

			for(int i=0;i<nelems/2;i++){
				TableRow row = new TableRow(this);
				TableRow row1 = new TableRow(this);
				TextView lowerSurface = new TextView(this);
				TextView upperSurface = new TextView(this);
				lowerSurface.setText(String.format("%.3f", rrc_scaledx[nelems/2-1-i]));
				lowerSurface.append("          ");
				upperSurface.setText(String.format("%.3f", rrc_scaledx[nelems/2+i]));
				upperSurface.append("          ");
				TextView cpLow = new TextView(this);
				TextView cpUp = new TextView(this);
				cpLow.setText(String.format("%.3f", cP[nelems/2-1-i]));
				cpUp.setText(String.format("%.3f", cP[nelems/2+i]));
				row.addView(lowerSurface);
				row.addView(cpLow);
				row1.addView(upperSurface);
				row1.addView(cpUp);
				table.addView(row);
				table1.addView(row1);
			}

			//The following lines of code contain the implementation of the Thwaites method for the
			//viscous flow corrections.

			int n_stg = 0;
			double csi_stg = 0.0;
			int i_stg = 0;

			for(int i = 0; i<nelems-1; i++){
				if(vTi[i]*vTi[i+1]<=0) {
					n_stg++;
					i_stg = i;
					csi_stg = -vTi[i]/(vTi[i+1] - vTi[i]);
				}
			}

			double[] s = new double[nelems];

			s[i_stg+1] = csi_stg*(len[i_stg] + len[i_stg+1]) * 0.5;
			s[i_stg] = (1.0-csi_stg) * (len[i_stg] + len[i_stg+1]) * 0.5;

			for(int i = i_stg+2; i < nelems; i++){
				s[i] = s[i-1] + 0.5*(len[i] + len[i-1]);
			}
			for(int i = i_stg-1; i>-1; i--){
				s[i] = s[i+1] + 0.5*(len[i] + len[i+1]);
			}

			double th = Math.acos(nvers[0][i_stg]*nvers[0][i_stg+1]+nvers[1][i_stg]*nvers[1][i_stg+1]);
			double r0 = len[1]*Math.sin(th) + (len[0] + len[1]*Math.cos(th))/Math.tan(th);
			double k = freeStream.v/r0;
			double theta0 = Math.sqrt(0.075*freeStream.kin_visc/k);  //never used

			double[] Ue = new double[vTi.length];
			for(int i = 0; i<vTi.length; i++){
				if(vTi[i] > 0)
					Ue[i] = vTi[i];
				else
					Ue[i] = -vTi[i];
			}

			double[] theta2Ue6 			= 		new double[nelems];
			double[] theta 				= 		new double[nelems];
			double[] ReTheta			=		new double[nelems];
			double[] lambda 			=		new double[nelems];
			double[] ell 				= 		new double[nelems];
			double[] H					=		new double[nelems];
			double[] delta				= 		new double[nelems];
			double[] cf					= 		new double[nelems];
			double[] UeThetaH1 			= 		new double[nelems];
			double[] H1					= 		new double[nelems];

			double[] Res				=		new double[Ue.length];
			double[] P					=		new double[cP.length];
			for(int i = 0; i<Ue.length;i++){
				Res[i] = Ue[i]*s[i] / freeStream.kin_visc;
				P[i] = cP[i]*(0.5*freeStream.rho*freeStream.v*freeStream.v);
			}

			double[] dPdx			= 		new double[nelems];
			double[] dUedx			= 		new double[nelems];

			double Ptot = freeStream.P + 0.5*freeStream.rho * freeStream.v * freeStream.v;

			dPdx[i_stg+1]	 = 	(P[i_stg+1] - Ptot) / (0.5*csi_stg*(len[i_stg]+ len[i_stg+1]));
			dPdx[i_stg]		 = 	(P[i_stg] - Ptot) / (0.5*(1.0-csi_stg)*(len[i_stg]+ len[i_stg+1]));
			dUedx[i_stg+1]	 = 	(Ue[i_stg+1] - 0.0) / (0.5*csi_stg*(len[i_stg]+ len[i_stg+1]));
			dUedx[i_stg]	 = 	(Ue[i_stg] -   0.0) / (0.5*(1.0-csi_stg)*(len[i_stg]+ len[i_stg+1]));

			theta2Ue6[i_stg+1] = 0.45*freeStream.kin_visc*Math.pow(Ue[i_stg+1],5.0)*
								 (0.25*csi_stg*(len[i_stg]+len[i_stg+1]));
			theta2Ue6[i_stg]   = 0.45*freeStream.kin_visc*Math.pow(Ue[i_stg],5.0)*
								 (0.25*(1.0-csi_stg)*(len[i_stg]+len[i_stg+1]));
			theta[i_stg+1]     = Math.sqrt(theta2Ue6[i_stg+1]/Math.pow(Ue[i_stg+1],6.0));
			theta[i_stg]       = Math.sqrt(theta2Ue6[i_stg]/Math.pow(Ue[i_stg],6.0));
			ReTheta[i_stg+1]   = Ue[i_stg+1]*theta[i_stg+1]/freeStream.kin_visc;
			ReTheta[i_stg]     = Ue[i_stg]*theta[i_stg]/freeStream.kin_visc;
			lambda[i_stg+1]    = Math.pow(theta[i_stg+1],2.0)*dUedx[i_stg+1]/freeStream.kin_visc;
			lambda[i_stg]      = Math.pow(theta[i_stg],2.0)*dUedx[i_stg]/freeStream.kin_visc;
			ell[i_stg+1] 	   = myGlobal.thwaites_ell(lambda[i_stg+1]);
			ell[i_stg] 	       = myGlobal.thwaites_ell(lambda[i_stg]);
			H[i_stg+1]		   = myGlobal.thwaites_H(lambda[i_stg+1]);
			H[i_stg]		   = myGlobal.thwaites_H(lambda[i_stg]);
			delta[i_stg+1]     = theta[i_stg+1] * H[i_stg+1];
			delta[i_stg]       = theta[i_stg] * H[i_stg];
			cf[i_stg+1]		   = 2.0*ell[i_stg+1]/ReTheta[i_stg+1];
			cf[i_stg]		   = 2.0*ell[i_stg]/ReTheta[i_stg];


			String regime = "laminar";
			String transition = "no";

			for(int i = i_stg+2; i<nelems; i++){

				if(i<nelems-1){
					dPdx[i]  = (P[i+1] - P[i-1]) / (0.5 * (len[i+1] + len[i-1]) + len[i]);
					dUedx[i] = ( Ue[i+1] - Ue[i-1] ) / ( 0.5 * ( len[i+1] + len[i-1] ) + len[i]);
				}
				else{
					dPdx[i]  = (P[i] - P[i-1]) / (0.5 * (len[i] + len[i-1]));
					dUedx[i] = ( Ue[i] - Ue[i-1] ) / ( 0.5 * ( len[i] + len[i-1] ));
				}

				if(regime.equals("laminar")) {
					theta2Ue6[i] = theta2Ue6[i - 1] + 0.45 * freeStream.kin_visc * (Math.pow(Ue[i], 5.0) +
							       Math.pow(Ue[i - 1], 5.0)) * (0.25 * (len[i] + len[i - 1]));

					theta[i]     = Math.sqrt(theta2Ue6[i] / Math.pow(Ue[i], 6.0));
					ReTheta[i]   = Ue[i] * theta[i] / freeStream.kin_visc;
					lambda[i]    = Math.pow(theta[i], 2.0) * dUedx[i] / freeStream.kin_visc;
					ell[i]       = myGlobal.thwaites_ell(lambda[i]);
					H[i]         = myGlobal.thwaites_H(lambda[i]);
					delta[i]     = theta[i] * H[i];
					cf[i]        = 2.0 * ell[i] / ReTheta[i];

					if (ReTheta[i] > 1.174 * (1.0 + 22400.0 / Res[i]) * Math.pow(Res[i], 0.46)) {
						regime = "turbulent";
						transition = "transition";
						H[i] = 1.35;
						delta[i] = theta[i] * H[i];
						H1[i] = myGlobal.head_HtoH1(H[i]);
						UeThetaH1[i] = Ue[i] * theta[i] * H1[i];
						int i_trans = i;
					}

					} else {

						double dx    = (0.5 * (len[i] + len[i - 1]));
						UeThetaH1[i] = UeThetaH1[i - 1] + dx * Ue[i - 1] * 0.0306 / Math.pow(H1[i - 1]
								       - 3.0, 0.6169);
						theta[i]     = theta[i - 1] + dx * (0.5 * cf[i - 1] - dUedx[i - 1] / Ue[i - 1]
								       * (2.0 + H[i - 1]) * theta[i - 1]);
						H1[i]        = Math.max(UeThetaH1[i] / (Ue[i] * theta[i]), 3.0 + 0.001);
						H[i]         = myGlobal.head_H1toH(H1[i]);
						delta[i]     = H[i] * theta[i];

						ReTheta[i]   = Ue[i] * theta[i] / freeStream.kin_visc;
						cf[i]        = myGlobal.ludweig_tillman_cf(H[i], ReTheta[i]);

						transition = "occurred";
					}

			}

			regime = "laminar";
			transition = "no";

			for(int i = i_stg-1;i>-1;i--){
				if(i>0){
					dPdx[i]  = (P[i-1]-P[i+1])/(0.5*(len[i+1]+len[i-1])+len[i]);
					dUedx[i] = (Ue[i-1]-Ue[i+1])/(0.5*(len[i+1]+len[i-1])+len[i]);
				}
				else{
					dPdx[i]  = (P[i]-P[i+1])/(0.5*(len[i]+len[i+1]));
					dUedx[i] = (Ue[i]-Ue[i+1])/(0.5*(len[i]+len[i+1]));
				}

				if(regime.equals("laminar")){
					theta2Ue6[i] = theta2Ue6[i + 1] + 0.45 * freeStream.kin_visc * (Math.pow(Ue[i], 5.0) +
							       Math.pow(Ue[i + 1], 5.0)) * (0.25 * (len[i] + len[i + 1]));
					theta[i]     = Math.sqrt(theta2Ue6[i] / Math.pow(Ue[i], 6.0));
					ReTheta[i]   = Ue[i] * theta[i] / freeStream.kin_visc;
					lambda[i]    = Math.pow(theta[i], 2.0) * dUedx[i] / freeStream.kin_visc;
					ell[i]       = myGlobal.thwaites_ell(lambda[i]);
					H[i]         = myGlobal.thwaites_H(lambda[i]);
					delta[i]     = theta[i] * H[i];
					cf[i]        = 2.0 * ell[i] / ReTheta[i];

				if (ReTheta[i] > 1.174 * (1.0 + 22400.0 / Res[i]) * Math.pow(Res[i], 0.46)) {
					regime        = "turbulent";
					transition    = "transition";
					H[i]          = 1.35;
					delta[i]      = theta[i] * H[i];
					H1[i]         = myGlobal.head_HtoH1(H[i]);
					UeThetaH1[i]  = Ue[i] * theta[i] * H1[i];
					int i_trans   = i;
				}

				}  else {

					double dx     = (0.5 * (len[i] + len[i + 1]));
					UeThetaH1[i]  = UeThetaH1[i + 1] + dx * Ue[i + 1] * 0.0306 / Math.pow(H1[i + 1]
							        - 3.0, 0.6169);
					theta[i]      = theta[i + 1] + dx * (0.5 * cf[i + 1] - dUedx[i + 1] / Ue[i + 1]
							        * (2.0 + H[i + 1]) * theta[i + 1]);
					H1[i]         = Math.max(UeThetaH1[i] / (Ue[i] * theta[i]), 3.0 + 0.001);
					H[i]          = myGlobal.head_H1toH(H1[i]);
					delta[i]      = H[i] * theta[i];

					ReTheta[i]    = Ue[i] * theta[i] / freeStream.kin_visc;
					cf[i]         = myGlobal.ludweig_tillman_cf(H[i], ReTheta[i]);

					transition    = "occurred";
				}
			}

			double[] tauW = new double[cf.length];

			for(int i = 0; i<cf.length;i++){
				tauW[i] = 0.5*freeStream.rho*Math.pow(freeStream.v,2)*cf[i];
			}

			double[][] temp = new double[2][tauW.length];
			double[][] dF_visc = new double[2][tauW.length];

			for(int i = 0; i<tauW.length;i++){
				temp[0][i] = tauW[i]*len[i]*vTi[i]/Ue[i];
				temp[1][i] = tauW[i]*len[i]*vTi[i]/Ue[i];
			}

			for(int i = 0; i<tauW.length;i++){
				dF_visc[0][i] = temp[0][i]*tvers[0][i];
				dF_visc[1][i] = temp[1][i]*tvers[1][i];
			}
			double[] F_visc = new double[2];

			F_visc[0] = DoubleStream.of(dF_visc[0]).sum();
			F_visc[1] = DoubleStream.of(dF_visc[1]).sum();

			double L_visc = -F_visc[0]*Math.sin(freeStream.alpha)+F_visc[1]*Math.cos(freeStream.alpha);
			double D_visc = F_visc[0]*Math.cos(freeStream.alpha) + F_visc[1]*Math.sin(freeStream.alpha);

			double[][] temp0 = new double[2][tauW.length];
			double[][] dF_Pres = new double[2][tauW.length];

			for(int i = 0; i<tauW.length;i++){
				temp0[0][i] = -P[i]*len[i];
				temp0[1][i] = -P[i]*len[i];
			}
			for(int i = 0; i<tauW.length;i++){
				dF_Pres[0][i] = temp0[0][i]*nvers[0][i];
				dF_Pres[1][i] = temp0[1][i]*nvers[1][i];
			}
			double[] F_pres = new double[2];

			F_pres[0] = DoubleStream.of(dF_Pres[0]).sum();
			F_pres[1] = DoubleStream.of(dF_Pres[1]).sum();

			double L_pres = -F_pres[0]*Math.sin(freeStream.alpha)+F_pres[1]*Math.cos(freeStream.alpha);
			double D_pres = F_pres[0]*Math.cos(freeStream.alpha) + F_pres[1]*Math.sin(freeStream.alpha);

			double L = L_pres + L_visc;
			double D = D_pres + D_visc;

			double cL = L/(0.5*freeStream.rho*freeStream.v*freeStream.v);
			double cD = D/(0.5*freeStream.rho*freeStream.v*freeStream.v);

			//Stampa a schermo dei risultati: CL, CLKJ, CL_thwaites, CD_thwaites

			TextView cltext = findViewById(R.id.textCl);
			cltext.setText(R.string.cl);
			cltext.append("   " + String.format("%.5f", cl) + " [-]");
			TextView cltextzukowski = findViewById(R.id.textClzukowski);
			cltextzukowski.setText(R.string.cl_zukowski);
			cltextzukowski.append("   " + String.format("%.5f", cLKJ) + "[-]");
			TextView cltextthwaites = findViewById(R.id.textCltwhaites);
			cltextthwaites.setText(R.string.cl_thwaites);
			cltextthwaites.append("   " + String.format("%.5f", cL)+ " [-]");
			TextView cdtextthwaites = findViewById(R.id.textCdtwhaites);
			cdtextthwaites.setText(R.string.cd_thwaites);
			cdtextthwaites.append("   " + String.format("%.5f", cD)+ " [-]");

	}
}
