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
import androidx.constraintlayout.widget.ConstraintLayout;

import android.view.View;
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

	String nacaDW;
	double chLnW, pr, den, spd, grndH;
	int nPan;
	boolean isGround;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_cl_alpha);
		FirebaseAnalytics mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		// The code written in this class is mostly copied from the BuildGeometry class. The only
		// difference is that no graph is drawn and the operations are repeated inside a for cycle
		// in which the angle of attack changes on every iteration.

		nacaDW   = "NACA";
		nacaDW   = nacaDW + getIntent().getStringExtra("NACA_DW");
		chLnW    = getIntent().getDoubleExtra("ChordW", 1.0);
		nPan     = getIntent().getIntExtra("nPan", 20);
		isGround = getIntent().getBooleanExtra("checkGround", false);
		grndH    = getIntent().getDoubleExtra("GroundH", 1.0);
		spd      = getIntent().getDoubleExtra("Speed", 10);

		pr  = 103125; //Pa
		den = 1.225; //kg/m3

		double[] refPoint = {0, 0};

		if(isGround)
			refPoint[1] = grndH;

		Airfoil airfoilW = new Airfoil(1, nacaDW, chLnW, -Math.toRadians(0.0),
				0.25, refPoint, nPan);

		GlobalFunctions myGlobal = new GlobalFunctions();

		double[][] rrW      = myGlobal.getCoordinates(airfoilW, false);
		int[][]    eeW      = myGlobal.connectivityMatrix(rrW[0]);
		int nelems          = Array.getLength(eeW[0]);
		Elems[] elemsW      = myGlobal.elemsMethod(airfoilW,rrW,eeW,nelems);
		Elems[] elemsW_mirr = elemsW;

		double[] lenW     = new double[nelems];
		double[][] tversW = new double[2][nelems];
		double[][] nversW = new double[2][nelems];

		for (int i = 0; i < nelems; i++) {
			lenW[i]      = elemsW[i].len;
			tversW[0][i] = elemsW[i].tver[0];
			tversW[1][i] = elemsW[i].tver[1];
			nversW[0][i] = elemsW[i].nver[0];
			nversW[1][i] = elemsW[i].nver[1];
		}

		if(isGround){
			double[][] rrW_mirr         = myGlobal.getCoordinates(airfoilW, true);
			int[][]    eeW_mirr         = myGlobal.connectivityMatrix(rrW_mirr[0]);
			elemsW_mirr                 = myGlobal.elemsMethod(airfoilW,rrW_mirr,eeW_mirr,nelems);
		}


		double[] cLW               = new double[21];
		double[] cL_zukowskiW      = new double[21];
		double[] cLW_mirr          = new double[21];
		double[] cL_zukowskiW_mirr = new double[21];
		int count = 0;

		for (int aOa = -10; aOa < 11; aOa++) {
			FreeStream freeStream = new FreeStream(pr, den, spd, Math.toRadians(aOa));

			double[][]  A, Au, Av;
			double[][] A_mirr, Au_mirr, Av_mirr;
			double[]    b;
			double[]    b_mirr;

			A  = myGlobal.coefficientMatrix(elemsW, nPan);
			b  = myGlobal.coefficientArray(freeStream, elemsW, nPan);
			Au = myGlobal.onBodyUMatrix(elemsW, nPan);
			Av = myGlobal.onBodyVMatrix(elemsW, nPan);

			A_mirr = A;
			b_mirr = b;
			Au_mirr = Au;
			Av_mirr = Av;

			if(isGround) {
				A_mirr = myGlobal.coefficientMatrixG(elemsW, elemsW_mirr, nPan);
				b_mirr = myGlobal.coefficientArray(freeStream, elemsW, nPan);
				Au_mirr = myGlobal.onBodyUMatrixG(elemsW, elemsW_mirr, nPan);
				Av_mirr = myGlobal.onBodyVMatrixG(elemsW, elemsW_mirr, nPan);
			}

			RealMatrix Acoefficients = new Array2DRowRealMatrix(A, false);
			DecompositionSolver solver = new LUDecomposition(Acoefficients).getSolver();
			RealVector bconstants = new ArrayRealVector(b, false);
			RealVector x = solver.solve(bconstants);
			RealVector u, v;
			RealMatrix AuCoeff = new Array2DRowRealMatrix(Au, false);
			RealMatrix AvCoeff = new Array2DRowRealMatrix(Av, false);

			RealMatrix uvec, uvecTi, uvecNi;
			uvec = new Array2DRowRealMatrix(nelems, 2);

			u = AuCoeff.operate(x);
			v = AvCoeff.operate(x);
			u = u.mapAdd(freeStream.v * Math.cos(freeStream.alpha));
			v = v.mapAdd(freeStream.v * Math.sin(freeStream.alpha));

			uvec.setColumnVector(0, u);
			uvec.setColumnVector(1, v);


			uvecTi = uvec.transpose().copy();
			uvecNi = uvec.transpose().copy();

			double[][] uvec_vTiW = new double[2][nelems];
			double[][] uvec_vNiW = new double[2][nelems];

			double[] vTiW = new double[nelems];
			double[] vNiW = new double[nelems];
			double[] vNi_checkW = new double[nelems];


			for (int j = 0; j < nelems; j++) {
				for (int i = 0; i < 2; i++) {
					uvecTi.multiplyEntry(i, j, tversW[i][j]);
					uvecNi.multiplyEntry(i, j, nversW[i][j]);
					uvec_vTiW[i][j] = uvecTi.getEntry(i, j);
					uvec_vNiW[i][j] = uvecNi.getEntry(i, j);
				}
				vTiW[j] = uvec_vTiW[0][j] + uvec_vTiW[1][j];
				vNiW[j] = uvec_vNiW[0][j] + uvec_vNiW[1][j];
			}

			for (int i = 0; i < nelems; i++) {
				if (vNiW[i] < 0)
					vNi_checkW[i] = -vNiW[i];
				else
					vNi_checkW[i] = vNiW[i];

				if (vNi_checkW[i] > Math.pow(10, -6))
					Toast.makeText(this, "Something is wrong", Toast.LENGTH_LONG).show();
			}

			double[] cPW = new double[nelems];

			for(int i = 0; i<nelems;i++){
				cPW[i] = 1.0 - vTiW[i]*vTiW[i]/(freeStream.v*freeStream.v);
			}

			double[][] cp_len = new double[2][cPW.length];

			for(int i = 0;  i<cPW.length; i++){
				cp_len[0][i] = cPW[i]*lenW[i];
				cp_len[1][i] = cPW[i]*lenW[i];
			}

			double[][] cp_len_nversW = new double[2][cPW.length];

			for(int i = 0;  i<cPW.length; i++){
				cp_len_nversW[0][i] = cp_len[0][i]*nversW[0][i];
				cp_len_nversW[1][i] = cp_len[1][i]*nversW[1][i];
			}

			double cfx = -DoubleStream.of(cp_len_nversW[0]).sum()/airfoilW.chord;
			double cfy = -DoubleStream.of(cp_len_nversW[1]).sum()/airfoilW.chord;

			double clW = -Math.sin(freeStream.alpha)*cfx + Math.cos(freeStream.alpha)*cfy;

			double sumW = DoubleStream.of(lenW).sum();
			double GammaW;

			GammaW =x.getEntry(nelems) * sumW;

			double liftKJW = -freeStream.rho * freeStream.v * GammaW;

			double cLKJW = liftKJW / (0.5 * freeStream.rho * Math.pow(freeStream.v, 2) * airfoilW.chord);

			cL_zukowskiW[count] = cLKJW;
			cLW[count] = clW;

			if(isGround){
				RealMatrix Acoefficients_mirr = new Array2DRowRealMatrix(A_mirr, false);
				DecompositionSolver solver_mirr = new LUDecomposition(Acoefficients_mirr).getSolver();
				RealVector bconstants_mirr = new ArrayRealVector(b_mirr, false);
				RealVector x_mirr = solver_mirr.solve(bconstants_mirr);
				RealVector u_mirr, v_mirr;
				RealMatrix AuCoeff_mirr = new Array2DRowRealMatrix(Au_mirr, false);
				RealMatrix AvCoeff_mirr = new Array2DRowRealMatrix(Av_mirr, false);

				RealMatrix uvec_mirr, uvecTi_mirr, uvecNi_mirr;
				uvec_mirr = new Array2DRowRealMatrix(nelems, 2);

				u_mirr = AuCoeff_mirr.operate(x_mirr);
				v_mirr = AvCoeff_mirr.operate(x_mirr);
				u_mirr = u_mirr.mapAdd(freeStream.v * Math.cos(freeStream.alpha));
				v_mirr = v_mirr.mapAdd(freeStream.v * Math.sin(freeStream.alpha));

				uvec_mirr.setColumnVector(0, u_mirr);
				uvec_mirr.setColumnVector(1, v_mirr);


				uvecTi_mirr = uvec_mirr.transpose().copy();
				uvecNi_mirr = uvec_mirr.transpose().copy();

				double[][] uvec_vTiW_mirr = new double[2][nelems];
				double[][] uvec_vNiW_mirr = new double[2][nelems];

				double[] vTiW_mirr       = new double[nelems];
				double[] vNiW_mirr       = new double[nelems];
				double[] vNi_checkW_mirr = new double[nelems];


				for (int j = 0; j < nelems; j++) {
					for (int i = 0; i < 2; i++) {
						uvecTi_mirr.multiplyEntry(i, j, tversW[i][j]);
						uvecNi_mirr.multiplyEntry(i, j, nversW[i][j]);
						uvec_vTiW_mirr[i][j] = uvecTi_mirr.getEntry(i, j);
						uvec_vNiW_mirr[i][j] = uvecNi_mirr.getEntry(i, j);
					}
					vTiW_mirr[j] = uvec_vTiW_mirr[0][j] + uvec_vTiW_mirr[1][j];
					vNiW_mirr[j] = uvec_vNiW_mirr[0][j] + uvec_vNiW_mirr[1][j];
				}

				for (int i = 0; i < nelems; i++) {
					if (vNiW_mirr[i] < 0)
						vNi_checkW_mirr[i] = -vNiW_mirr[i];
					else
						vNi_checkW_mirr[i] = vNiW_mirr[i];

					if (vNi_checkW_mirr[i] > Math.pow(10, -6)) {
						Toast.makeText(this, "Something is wrong!!!!!", Toast.LENGTH_SHORT).show();
					}
				}

				double[] cPW_mirr = new double[nelems];

				for(int i = 0; i<nelems;i++){
					cPW_mirr[i] = 1.0 - vTiW_mirr[i]*vTiW_mirr[i]/(freeStream.v*freeStream.v);
				}

				double[][] cp_len_mirr = new double[2][cPW.length];

				for(int i = 0;  i<cPW.length; i++){
					cp_len_mirr[0][i] = cPW_mirr[i]*lenW[i];
					cp_len_mirr[1][i] = cPW_mirr[i]*lenW[i];
				}

				double[][] cp_len_nversW_mirr = new double[2][cPW_mirr.length];

				for(int i = 0;  i<cPW_mirr.length; i++){
					cp_len_nversW_mirr[0][i] = cp_len_mirr[0][i]*nversW[0][i];
					cp_len_nversW_mirr[1][i] = cp_len_mirr[1][i]*nversW[1][i];
				}

				double cfx_mirr = -DoubleStream.of(cp_len_nversW_mirr[0]).sum()/airfoilW.chord;
				double cfy_mirr = -DoubleStream.of(cp_len_nversW_mirr[1]).sum()/airfoilW.chord;

				double clW_mirr = -Math.sin(freeStream.alpha)*cfx_mirr + Math.cos(freeStream.alpha)*cfy_mirr;

				double sumW_mirr = DoubleStream.of(lenW).sum();
				double GammaW_mirr;

				GammaW_mirr =x_mirr.getEntry(nelems) * sumW_mirr;

				double liftKJW_mirr = -freeStream.rho * freeStream.v * GammaW_mirr;

				double cLKJW_mirr = liftKJW_mirr / (0.5 * freeStream.rho * Math.pow(freeStream.v, 2) * airfoilW.chord);

				cL_zukowskiW_mirr[count] = cLKJW_mirr;
				cLW_mirr[count] = clW_mirr;
				}
			count++;
		}

		int[] alpha = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};

		final GraphView graph = findViewById(R.id.graph);

		DataPoint[] values   = new DataPoint[21];
		DataPoint[] values1  = new DataPoint[21];
		DataPoint[] values2  = new DataPoint[21];
		DataPoint[] values3  = new DataPoint[21];
		for(int i =0; i<21; i++){
			DataPoint value = new DataPoint(alpha[i], cLW[i]);
			DataPoint value1 = new DataPoint(alpha[i], cL_zukowskiW[i]);
			values[i] = value;
			values1[i] = value1;
			if(isGround){
				DataPoint value2 = new DataPoint(alpha[i], cLW_mirr[i]);
				DataPoint value3 = new DataPoint(alpha[i], cL_zukowskiW_mirr[i]);
				values2[i] = value2;
				values3[i] = value3;
			}
		}

		LineGraphSeries<DataPoint> series = new LineGraphSeries<>(values);
		LineGraphSeries<DataPoint> series1 = new LineGraphSeries<>(values1);
		graph.addSeries(series);
		graph.addSeries(series1);
		series.setColor(Color.rgb(0,96,100));
		series1.setColor(Color.rgb(230,81,0));
		series.setThickness(3);
		series1.setThickness(3);

		if(isGround){
			LineGraphSeries<DataPoint> series2 = new LineGraphSeries<>(values2);
			LineGraphSeries<DataPoint> series3 = new LineGraphSeries<>(values3);
			graph.addSeries(series2);
			graph.addSeries(series3);
			series2.setColor(Color.rgb(124,179,66));
			series3.setColor(Color.rgb(255,171,0));
			series2.setThickness(3);
			series3.setThickness(3);
			}

		TableLayout table = findViewById(R.id.table);
		TableLayout table1 = findViewById(R.id.table1);

		TableRow row_header = new TableRow(this);
		TableRow row_header1 = new TableRow(this);
		TextView angle = new TextView(this);
		TextView angle0 = new TextView(this);
		angle.setText(R.string.angle1);
		angle.append("   ");
		angle.setPadding(0,0,0,10);
		angle0.setText(R.string.angle1);
		angle0.append("   ");
		angle0.setPadding(0,0,0,10);
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

		for(int i=0;i<cLW.length;i++)
		{
			TableRow row1 = new TableRow(this);
			TableRow row2 = new TableRow(this);
			TextView angle1 =new TextView(this);
			TextView angle2 =new TextView(this);
			angle1.setText("" + alpha[i]);
			angle2.setText("" + alpha[i]);
			TextView clalpha1 =new TextView(this);
			TextView clalpha2 =new TextView(this);
			clalpha1.setText(String.format("%.5f", cLW[i]));
			clalpha2.setText(String.format("%.5f", cL_zukowskiW[i]));
			row1.addView(angle1);
			row2.addView(angle2);
			row1.addView(clalpha1);
			row2.addView(clalpha2);
			table.addView(row1);
			table1.addView(row2);
		}

		if(isGround){
			TextView textView = findViewById(R.id.titolo);
			textView.setVisibility(View.VISIBLE);
			textView = findViewById(R.id.titolo1);
			textView.setVisibility(View.VISIBLE);
			ConstraintLayout constraintLayout = findViewById(R.id.all01);
			constraintLayout.setVisibility(View.VISIBLE);
			constraintLayout = findViewById(R.id.all11);
			constraintLayout.setVisibility(View.VISIBLE);
			constraintLayout = findViewById(R.id.s2);
			constraintLayout.setVisibility(View.VISIBLE);
			constraintLayout = findViewById(R.id.d2);
			constraintLayout.setVisibility(View.VISIBLE);

			table = findViewById(R.id.table_mirr);
			table1 = findViewById(R.id.table1_mirr);

			row_header = new TableRow(this);
			row_header1 = new TableRow(this);
			angle = new TextView(this);
			angle0 = new TextView(this);
			angle.setText(R.string.angle1);
			angle.append("   ");
			angle.setPadding(0,0,0,10);
			angle0.setText(R.string.angle1);
			angle0.append("   ");
			angle0.setPadding(0,0,0,10);
			cL_text = new TextView(this);
			cLKJ_text = new TextView(this);
			cL_text.setText(R.string.cl1);
			cLKJ_text.setText(R.string.cl2);
			row_header.addView(angle);
			row_header1.addView(angle0);
			row_header.addView(cL_text);
			row_header1.addView(cLKJ_text);
			table.addView(row_header);
			table1.addView(row_header1);

			for(int i=0;i<cLW.length;i++)
			{
				TableRow row1 = new TableRow(this);
				TableRow row2 = new TableRow(this);
				TextView angle1 =new TextView(this);
				TextView angle2 =new TextView(this);
				angle1.setText("" + alpha[i]);
				angle2.setText("" + alpha[i]);
				TextView clalpha1 =new TextView(this);
				TextView clalpha2 =new TextView(this);
				clalpha1.setText(String.format("%.5f", cLW_mirr[i]));
				clalpha2.setText(String.format("%.5f", cL_zukowskiW_mirr[i]));
				row1.addView(angle1);
				row2.addView(angle2);
				row1.addView(clalpha1);
				row2.addView(clalpha2);
				table.addView(row1);
				table1.addView(row2);
			}
		}
	}
}
