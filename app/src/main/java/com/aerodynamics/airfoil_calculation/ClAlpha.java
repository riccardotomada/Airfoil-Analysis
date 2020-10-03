package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//This class contains all the Java code for the Cl_alpha screen                  			      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import android.graphics.Color;
import android.os.Bundle;
import com.google.firebase.analytics.FirebaseAnalytics;
import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;
import androidx.appcompat.app.AppCompatActivity;
import android.widget.TableLayout;
import android.widget.TableRow;
import android.widget.TextView;
import android.widget.Toast;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import java.lang.reflect.Array;
import java.util.stream.DoubleStream;

public class ClAlpha extends AppCompatActivity {

	String nacaD;
	GlobalFunctions myGlobal;
	double chLn;
	int nPan;
	double pr;
	double den;
	double spd;
	double aOA;
	private FirebaseAnalytics mFirebaseAnalytics;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_cl_alpha);
		mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		// The code written in this class is mostly copied from the BuildGeometry class. The only
		// difference is that no graph is drawn and the operations are repeated inside a for cycle
		// in which the angle of attack changes on every iteration.

		nacaD = "NACA";
		nacaD = nacaD + getIntent().getStringExtra("NACA_D");
		chLn = getIntent().getDoubleExtra("Chord", 1.0);
		nPan = getIntent().getIntExtra("nPan", 20);
		pr = getIntent().getDoubleExtra("Pressure", 10325.15);
		den = getIntent().getDoubleExtra("Density", 1.225);
		spd = getIntent().getDoubleExtra("Speed", 10);
		aOA = getIntent().getDoubleExtra("AngleOA", 0.0);

		double[] refPoint = {0, 0};

		Airfoil airfoil = new Airfoil(1, nacaD, chLn, -Math.toRadians(0.0),
				0.25, refPoint, nPan);
		GlobalFunctions myGlobal = new GlobalFunctions();
		double[][] rr = myGlobal.getCoordinates(airfoil);

		int[][] ee = myGlobal.connectivityMatrix(rr[0]);

		int nelems = Array.getLength(ee[0]);

		Elems[] elems = new Elems[nelems];
		for (int i = 0; i < nelems; i++) {
			elems[i] = new Elems();
			elems[i].airfoilId = airfoil.id;
			elems[i].id = i;
			elems[i].ver1 = new double[]{rr[0][ee[0][i] - 1], rr[1][ee[0][i] - 1]};
			elems[i].ver2 = new double[]{rr[0][ee[1][i] - 1], rr[1][ee[1][i] - 1]};
			elems[i].cen = myGlobal.midpoint(elems[i].ver1, elems[i].ver2);
			elems[i].len = myGlobal.len_norm(elems[i].ver2, elems[i].ver1);
			elems[i].tver = myGlobal.tver(elems[i].ver2, elems[i].ver1, elems[i].len);
			elems[i].nver = new double[]{-elems[i].tver[1], elems[i].tver[0]};
		}
		double[] len = new double[nelems];
		double[][] tvers = new double[2][nelems];
		double[][] nvers = new double[2][nelems];

		for (int i = 0; i < nelems; i++) {
			len[i] = elems[i].len;
			tvers[0][i] = elems[i].tver[0];
			tvers[1][i] = elems[i].tver[1];
			nvers[0][i] = elems[i].nver[0];
			nvers[1][i] = elems[i].nver[1];
		}


		double[] cL = new double[21];
		double[] cL_zukowski = new double[21];
		int count = 0;

		for (int aOa = -10; aOa < 11; aOa++) {
			FreeStream freeStream = new FreeStream(pr, den, spd, Math.toRadians(aOa));

			double[][] A = myGlobal.coefficientMatrix(elems, nPan);
			double[] b = myGlobal.coefficientArray(freeStream, elems, nPan);
			double[][] Au = myGlobal.onBodyUMatrix(elems, nPan);
			double[][] Av = myGlobal.onBodyVMatrix(elems, nPan);

			RealMatrix Acoefficients = new Array2DRowRealMatrix(A, false);

			DecompositionSolver solver = new LUDecomposition(Acoefficients).getSolver();

			RealVector bconstants = new ArrayRealVector(b, false);
			RealVector x = solver.solve(bconstants);

			RealVector u, v;

			RealMatrix AuCoeff = new Array2DRowRealMatrix(Au, false);
			RealMatrix AvCoeff = new Array2DRowRealMatrix(Av, false);

			RealMatrix uvec = new Array2DRowRealMatrix(nelems, 2);
			RealMatrix uvecTi = new Array2DRowRealMatrix(2, nelems);
			RealMatrix uvecNi = new Array2DRowRealMatrix(2, nelems);

			u = AuCoeff.operate(x);
			v = AvCoeff.operate(x);
			u = u.mapAdd(freeStream.v * Math.cos(freeStream.alpha));
			v = v.mapAdd(freeStream.v * Math.sin(freeStream.alpha));

			uvec.setColumnVector(0, u);
			uvec.setColumnVector(1, v);


			uvecTi = uvec.transpose().copy();
			uvecNi = uvec.transpose().copy();

			double[][] uvec_vTi = new double[2][nelems];
			double[][] uvec_vNi = new double[2][nelems];

			double[] vTi = new double[nelems];
			double[] vNi = new double[nelems];
			double[] vNi_check = new double[nelems];


			for (int j = 0; j < nelems; j++) {
				for (int i = 0; i < 2; i++) {
					uvecTi.multiplyEntry(i, j, tvers[i][j]);
					uvecNi.multiplyEntry(i, j, nvers[i][j]);
					uvec_vTi[i][j] = uvecTi.getEntry(i, j);
					uvec_vNi[i][j] = uvecNi.getEntry(i, j);
				}
				vTi[j] = uvec_vTi[0][j] + uvec_vTi[1][j];
				vNi[j] = uvec_vNi[0][j] + uvec_vNi[1][j];
			}

			for (int i = 0; i < nelems; i++) {
				if (vNi[i] < 0)
					vNi_check[i] = -vNi[i];
				else
					vNi_check[i] = vNi[i];

				if (vNi_check[i] > Math.pow(10, -6))
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
			double Gamma = x.getEntry(nelems) * sum;

			double liftKJ = -freeStream.rho * freeStream.v * Gamma;
			double cLKJ = liftKJ / (0.5 * freeStream.rho * Math.pow(freeStream.v, 2) * airfoil.chord);

			cL_zukowski[count] = cl;
			cL[count] = cLKJ;

			count++;
		}

		int[] alpha = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};

		final GraphView graph = findViewById(R.id.graph);

		DataPoint[] values = new DataPoint[21];
		DataPoint[] values1 = new DataPoint[21];
		for(int i =0; i<21; i++){
			DataPoint value = new DataPoint(alpha[i], cL[i]);
			DataPoint value1 = new DataPoint(alpha[i], cL_zukowski[i]);
			values[i] = value;
			values1[i] = value1;
		}

		LineGraphSeries<DataPoint> series = new LineGraphSeries<DataPoint>(values);
		LineGraphSeries<DataPoint> series1 = new LineGraphSeries<DataPoint>(values1);
		graph.addSeries(series);
		graph.addSeries(series1);
		series.setColor(Color.rgb(0,151,167));
		series1.setColor(Color.rgb(255,109,0));

		TableLayout table = findViewById(R.id.table);
		TableLayout table1 = findViewById(R.id.table1);

		TableRow row_header = new TableRow(this);
		TableRow row_header1 = new TableRow(this);
		TextView angle = new TextView(this);
		TextView angle0 = new TextView(this);
		angle.setText(R.string.angle1);
		angle0.setText(R.string.angle1);
		TextView cL_text = new TextView(this);
		TextView cLKJ_text = new TextView(this);
		cL_text.setText(R.string.cl1);
		cLKJ_text.setText(R.string.cl2);
		row_header.addView(angle);
		row_header1.addView(angle0);
		row_header.addView(cL_text);
		row_header1.addView(cLKJ_text);
		table.addView(row_header);
		table1.addView(row_header1);

		for(int i=0;i<cL.length;i++)
		{
			TableRow row1 = new TableRow(this);
			TableRow row2 = new TableRow(this);
			TextView angle1 =new TextView(this);
			TextView angle2 =new TextView(this);
			angle1.setText("" + alpha[i]);
			angle2.setText("" + alpha[i]);
			angle1.append("                    ");
			angle2.append("                    ");
			TextView clalpha1 =new TextView(this);
			TextView clalpha2 =new TextView(this);
			clalpha1.setText(String.format("%.5f", cL[i]));
			clalpha2.setText(String.format("%.5f", cL_zukowski[i]));
			row1.addView(angle1);
			row2.addView(angle2);
			row1.addView(clalpha1);
			row2.addView(clalpha2);
			table.addView(row1);
			table1.addView(row2);
		}
	}
}
