package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//This class contains all the Java code about the Geometry screen   					          //
////////////////////////////////////////////////////////////////////////////////////////////////////

import static java.lang.Double.max;
import static java.lang.Double.min;
import static java.lang.Math.abs;

import android.view.View;

////////////////////////////////////////////////////////////////////////////////////////////////////
//This class contains all the Java code about the "CP, CL, CD" screen   					      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import androidx.appcompat.app.AppCompatActivity;
import androidx.constraintlayout.widget.ConstraintLayout;

import android.graphics.Color;
import android.os.Bundle;
import android.widget.TableLayout;
import android.widget.TableRow;
import android.widget.TextView;

import com.google.firebase.analytics.FirebaseAnalytics;
import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;

import java.lang.reflect.Array;
import java.util.Arrays;

public class BuildGeometry extends AppCompatActivity {

	String          nacaDW, nacaDT;
	double          chLnW, chLnT, staggW, staggT, grndH, wTX, wTY;
	int             nPan;
	boolean         isWT, isGround;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_build_geometry);
		FirebaseAnalytics mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		// Receiving the values typed by the user in the MainActivity

		nacaDW  = "NACA";
		nacaDW  = nacaDW + getIntent().getStringExtra("NACA_DW");
		chLnW   = getIntent().getDoubleExtra("ChordW", 1.0);
		staggW  = getIntent().getDoubleExtra("StaggerW", 0.0);
		nPan    = getIntent().getIntExtra("nPan", 20);
		nacaDT  = "NACA";
		nacaDT  = nacaDT + getIntent().getStringExtra("NACA_DT");
		chLnT   = getIntent().getDoubleExtra("ChordT", 1.0);
		staggT  = getIntent().getDoubleExtra("StaggerT", 0.0);
		isGround = getIntent().getBooleanExtra("checkGround", false);
		grndH    = getIntent().getDoubleExtra("GroundH", 0.0);
		isWT    = getIntent().getBooleanExtra("checkTail", false);
		wTX     = getIntent().getDoubleExtra("wTX",1.5);
		wTY     = getIntent().getDoubleExtra("wTY",0.2);

		//From here on all the code it's quite similar to the one implemented on MATLAB during
		//fluid dynamics and aerodynamics courses.

		double[] refPoint = {0, 0};

		Airfoil airfoilW = new Airfoil(1, nacaDW, chLnW, -Math.toRadians(staggW),
				0.25, refPoint, nPan);

		// Since many functions are meant to be used many times and shared with other activities,
		// I decided to put all of them inside a new class, named GlobalFunctions. In order to use
		// them a GlobalFunctions object, myGlobal, initialized as follows, will be called every time
		// it will be necessary.

		GlobalFunctions myGlobal = new GlobalFunctions();

		//I tried to use the same variable names as we already wrote on MATLAB.

		double[][] rr_wing     = myGlobal.getCoordinates(airfoilW, false);
		int[][]    ee_wing     = myGlobal.connectivityMatrix(rr_wing[0]);
		int        nelems_wing = Array.getLength(ee_wing[0]);
		Elems[]    elems_wing  = new Elems[nelems_wing];

		for (int i = 0; i < nelems_wing; i++) {
			elems_wing[i]           = new Elems();
			elems_wing[i].airfoilId = airfoilW.id;
			elems_wing[i].id        = i;
			elems_wing[i].ver1      = new double[]{rr_wing[0][ee_wing[0][i] - 1], rr_wing[1][ee_wing[0][i] - 1]};
			elems_wing[i].ver2      = new double[]{rr_wing[0][ee_wing[1][i] - 1], rr_wing[1][ee_wing[1][i] - 1]};
			elems_wing[i].cen       = myGlobal.midpoint(elems_wing[i].ver1, elems_wing[i].ver2);
			elems_wing[i].len       = myGlobal.len_norm(elems_wing[i].ver2, elems_wing[i].ver1);
			elems_wing[i].tver      = myGlobal.tver(elems_wing[i].ver2, elems_wing[i].ver1, elems_wing[i].len);
			elems_wing[i].nver      = new double[]{-elems_wing[i].tver[1], elems_wing[i].tver[0]};
		}

		double[][] rrc   = new double[2][nelems_wing];
		double[] len     = new double[nelems_wing];
		double[][] tvers = new double[2][nelems_wing];
		double[][] nvers = new double[2][nelems_wing];

		for (int i = 0; i < nelems_wing; i++) {
			rrc[0][i]   = elems_wing[i].cen[0];
			rrc[1][i]   = elems_wing[i].cen[1];
			len[i]      = elems_wing[i].len;
			tvers[0][i] = elems_wing[i].tver[0];
			tvers[1][i] = elems_wing[i].tver[1];
			nvers[0][i] = elems_wing[i].nver[0];
			nvers[1][i] = elems_wing[i].nver[1];
		}

		double[] rrc_scaledx = new double[nelems_wing];
		for (int i = 0; i < nelems_wing; i++) {
			rrc_scaledx[i] = rrc[0][i] / chLnW;
		}

		double[] rrc_scaledy = new double[nelems_wing];
		for (int i = 0; i < nelems_wing; i++) {
			rrc_scaledy[i] = rrc[1][i] / chLnW;
		}

		double[] rrc_scaledx1 = new double[nelems_wing / 2];
		double[] rrc_scaledx2 = new double[nelems_wing / 2];
		double[] rrc_scaledy1 = new double[nelems_wing / 2];
		double[] rrc_scaledy2 = new double[nelems_wing / 2];
		System.arraycopy(rrc_scaledx, 0, rrc_scaledx1, 0, nelems_wing / 2);
		if (nelems_wing - nelems_wing / 2 >= 0)
			System.arraycopy(rrc_scaledx, nelems_wing / 2, rrc_scaledx2, nelems_wing / 2 - nelems_wing / 2, nelems_wing - nelems_wing / 2);
		System.arraycopy(rrc_scaledy, 0, rrc_scaledy1, 0, nelems_wing / 2);
		if (nelems_wing - nelems_wing / 2 >= 0)
			System.arraycopy(rrc_scaledy, nelems_wing / 2, rrc_scaledy2, nelems_wing / 2 - nelems_wing / 2, nelems_wing - nelems_wing / 2);

		// Unfortunately, all the graph libraries I have found don't allow to plot anything if
		// the x array is not ordered from the smallest to the greatest value.

		// In the following lines of code I reordered the x array keeping unaltered all the
		// correspondences with the other arrays (y array for the airfoil drawing).

		// This means that if two elements of the x array are switched in position, the same
		// happens for the correspondent elements of the other arrays.

		// I have to admit that this is a cumbersome, complicated and probably not efficient way
		// to rearrange the arrays, but this is the only way I've been able to overcome this problem.

		boolean swapped = true;
		int j = 0;
		double temp1, temp2;

		while (swapped) {
			swapped = false;
			j++;
			for (int i = 0; i < nelems_wing / 2 - j; i++) {
				if (rrc_scaledx1[i] > rrc_scaledx1[i + 1]) {
					temp1               = rrc_scaledx1[i];
					temp2               = rrc_scaledy1[i];
					rrc_scaledx1[i]     = rrc_scaledx1[i + 1];
					rrc_scaledy1[i]     = rrc_scaledy1[i + 1];
					rrc_scaledx1[i + 1] = temp1;
					rrc_scaledy1[i + 1] = temp2;
					swapped             = true;
				}
				if (rrc_scaledx2[i] > rrc_scaledx2[i + 1]) {
					temp1               = rrc_scaledx2[i];
					temp2               = rrc_scaledy2[i];
					rrc_scaledx2[i]     = rrc_scaledx2[i + 1];
					rrc_scaledy2[i]     = rrc_scaledy2[i + 1];
					rrc_scaledx2[i + 1] = temp1;
					rrc_scaledy2[i + 1] = temp2;
					swapped             = true;
				}
			}
		}

		// The following lines of code refer to the airfoil graph construction.

		DataPoint[] values  = new DataPoint[nelems_wing / 2];
		DataPoint[] values1 = new DataPoint[nelems_wing / 2];

		for (int i = 0; i < nelems_wing / 2; i++) {

			DataPoint value, value1;

			value = new DataPoint(rrc_scaledx1[i], rrc_scaledy1[i]);
			value1 = new DataPoint(rrc_scaledx2[i], rrc_scaledy2[i]);

			values[i]  = value;
			values1[i] = value1;
		}

		GraphView graph = findViewById(R.id.graph);

		LineGraphSeries<DataPoint> series  = new LineGraphSeries<>(values);
		LineGraphSeries<DataPoint> series1 = new LineGraphSeries<>(values1);
		graph.addSeries(series);
		graph.addSeries(series1);
		double max_x = 1.0;
		double max_y = Arrays.stream(rr_wing[1]).max().getAsDouble()/chLnW;
		double min_y = Arrays.stream(rr_wing[1]).min().getAsDouble()/chLnW;
		graph.getViewport().setMinY(-3 * abs(min_y));
		graph.getViewport().setMaxY(3 * abs(max_y));
		graph.getViewport().setMaxX(max_x);
		series.setDrawDataPoints(true);
		series.setDataPointsRadius(6);
		series1.setDrawDataPoints(true);
		series1.setDataPointsRadius(6);
		graph.getViewport().setXAxisBoundsManual(true);
		graph.getViewport().setYAxisBoundsManual(true);
		graph.getViewport().setScalableY(false);
		if(isWT)
			graph.setTitle(getString(R.string.wing) + ": " + nacaDW);
		else
			graph.setTitle(getString(R.string.airfoilchosen) + ": " + nacaDW);
		graph.setTitleTextSize(40);
		series.setColor(Color.rgb(0, 145, 234));
		series1.setColor(Color.rgb(0, 145, 234));

		// The code below refers to 2 tables, which show the y value of the airfoil as a function
		// of the position respectively on the lower and upper surface of the airfoil.

		TableLayout table  = findViewById(R.id.table);
		TableLayout table1 = findViewById(R.id.table1);

		TableRow row_header  = new TableRow(this);
		TableRow row_header1 = new TableRow(this);
		TextView x_c         = new TextView(this);
		x_c.setText(R.string.x_c);
		TextView x_c1        = new TextView(this);
		x_c1.setText(R.string.x_c);
		TextView yc_text     = new TextView(this);
		yc_text.setText(R.string.y_c);
		TextView yc_text1    = new TextView(this);
		yc_text1.setText(R.string.y_c);
		row_header.addView(x_c);
		row_header.addView(yc_text);
		row_header1.addView(x_c1);
		row_header1.addView(yc_text1);
		table.addView(row_header);
		table1.addView(row_header1);

		for (int i = 0; i < nelems_wing / 2; i++) {
			TableRow row          = new TableRow(this);
			TableRow row1         = new TableRow(this);
			TextView lowerSurface = new TextView(this);
			TextView upperSurface = new TextView(this);
			lowerSurface.setText(String.format("%.3f", rrc_scaledx[nelems_wing / 2 - 1 - i]));
			lowerSurface.append("          ");
			upperSurface.setText(String.format("%.3f", rrc_scaledx[nelems_wing / 2 + i]));
			upperSurface.append("          ");
			TextView ycLow        = new TextView(this);
			TextView ycUp         = new TextView(this);
			ycLow.setText(String.format("%.3f", rrc_scaledy[nelems_wing / 2 - 1 - i]));
			ycUp.setText(String.format("%.3f", rrc_scaledy[nelems_wing / 2 + i]));
			row.addView(lowerSurface);
			row.addView(ycLow);
			row1.addView(upperSurface);
			row1.addView(ycUp);
			table.addView(row);
			table1.addView(row1);
		}
		if(isWT) {
			Airfoil airfoilT = new Airfoil(2, nacaDT, chLnT, -Math.toRadians(staggT),
					0.25, refPoint, nPan);

			graph = findViewById(R.id.graph_tail);
			graph.setVisibility(View.VISIBLE);
			GraphView graph_assembly= findViewById(R.id.graph_assembly);
			graph_assembly.setVisibility(View.VISIBLE);
			ConstraintLayout constraintLayout = findViewById(R.id.s2);
			constraintLayout.setVisibility(View.VISIBLE);
			ConstraintLayout constraintLayout1 = findViewById(R.id.d2);
			constraintLayout1.setVisibility(View.VISIBLE);
			TextView textView = findViewById(R.id.titolo);
			textView.setVisibility(View.VISIBLE);
			textView = findViewById(R.id.titolo1);
			textView.setVisibility(View.VISIBLE);

			double[][] rr_tail = myGlobal.getCoordinates(airfoilT, false);
			int[][] ee_tail    = myGlobal.connectivityMatrix(rr_tail[0]);
			int nelems_tail    = Array.getLength(ee_tail[0]);
			Elems[] elems_tail = new Elems[nelems_tail];

			for (int i = 0; i < nelems_tail; i++) {
				elems_tail[i] = new Elems();
				elems_tail[i].airfoilId = airfoilT.id;
				elems_tail[i].id = i;
				elems_tail[i].ver1 = new double[]{rr_tail[0][ee_tail[0][i] - 1], rr_tail[1][ee_tail[0][i] - 1]};
				elems_tail[i].ver2 = new double[]{rr_tail[0][ee_tail[1][i] - 1], rr_tail[1][ee_tail[1][i] - 1]};
				elems_tail[i].cen = myGlobal.midpoint(elems_tail[i].ver1, elems_tail[i].ver2);
				elems_tail[i].len = myGlobal.len_norm(elems_tail[i].ver2, elems_tail[i].ver1);
				elems_tail[i].tver = myGlobal.tver(elems_tail[i].ver2, elems_tail[i].ver1, elems_tail[i].len);
				elems_tail[i].nver = new double[]{-elems_tail[i].tver[1], elems_tail[i].tver[0]};
			}

			double[][] rrc_tail   = new double[2][nelems_tail];
			double[]   len_tail   = new double[nelems_tail];
			double[][] tvers_tail = new double[2][nelems_tail];
			double[][] nvers_tail = new double[2][nelems_tail];

			for (int i = 0; i < nelems_tail; i++) {
				rrc_tail[0][i] = elems_tail[i].cen[0];
				rrc_tail[1][i] = elems_tail[i].cen[1];
				len_tail[i] = elems_tail[i].len;
				tvers_tail[0][i] = elems_tail[i].tver[0];
				tvers_tail[1][i] = elems_tail[i].tver[1];
				nvers_tail[0][i] = elems_tail[i].nver[0];
				nvers_tail[1][i] = elems_tail[i].nver[1];
			}

			double[] rrc_scaledx_tail = new double[nelems_tail];
			for (int i = 0; i < nelems_tail; i++) {
				rrc_scaledx_tail[i] = rrc_tail[0][i] / chLnT;
			}

			double[] rrc_scaledy_tail = new double[nelems_tail];
			for (int i = 0; i < nelems_tail; i++) {
				rrc_scaledy_tail[i] = rrc_tail[1][i] / chLnT;
			}

			double[] rrc_scaledx1_tail = new double[nelems_tail / 2];
			double[] rrc_scaledx2_tail = new double[nelems_tail / 2];
			double[] rrc_scaledy1_tail = new double[nelems_tail / 2];
			double[] rrc_scaledy2_tail = new double[nelems_tail / 2];
			System.arraycopy(rrc_scaledx_tail, 0, rrc_scaledx1_tail, 0, nelems_tail / 2);
			if (nelems_tail - nelems_tail / 2 >= 0)
				System.arraycopy(rrc_scaledx_tail, nelems_tail / 2, rrc_scaledx2_tail, nelems_tail / 2 - nelems_tail / 2, nelems_tail - nelems_tail / 2);
			System.arraycopy(rrc_scaledy_tail, 0, rrc_scaledy1_tail, 0, nelems_tail / 2);
			if (nelems_tail - nelems_tail / 2 >= 0)
				System.arraycopy(rrc_scaledy_tail, nelems_tail / 2, rrc_scaledy2_tail, nelems_tail / 2 - nelems_tail / 2, nelems_tail - nelems_tail / 2);

			swapped = true;
			j = 0;

			while (swapped) {
				swapped = false;
				j++;
				for (int i = 0; i < nelems_tail / 2 - j; i++) {
					if (rrc_scaledx1_tail[i] > rrc_scaledx1_tail[i + 1]) {
						temp1                    = rrc_scaledx1_tail[i];
						temp2                    = rrc_scaledy1_tail[i];
						rrc_scaledx1_tail[i]     = rrc_scaledx1_tail[i + 1];
						rrc_scaledy1_tail[i]     = rrc_scaledy1_tail[i + 1];
						rrc_scaledx1_tail[i + 1] = temp1;
						rrc_scaledy1_tail[i + 1] = temp2;
						swapped                  = true;
					}
					if (rrc_scaledx2_tail[i] > rrc_scaledx2_tail[i + 1]) {
						temp1                    = rrc_scaledx2_tail[i];
						temp2                    = rrc_scaledy2_tail[i];
						rrc_scaledx2_tail[i]     = rrc_scaledx2_tail[i + 1];
						rrc_scaledy2_tail[i]     = rrc_scaledy2_tail[i + 1];
						rrc_scaledx2_tail[i + 1] = temp1;
						rrc_scaledy2_tail[i + 1] = temp2;
						swapped                  = true;
					}
				}
			}

			// The following lines of code refer to the airfoil graph construction.

			values  = new DataPoint[nelems_tail / 2];
			values1 = new DataPoint[nelems_tail / 2];

			DataPoint[] valuesAssemblyWing  = new DataPoint[nelems_wing / 2];
			DataPoint[] valuesAssemblyWing1 = new DataPoint[nelems_wing / 2];
			DataPoint[] valuesAssemblyTail  = new DataPoint[nelems_tail / 2];
			DataPoint[] valuesAssemblyTail1 = new DataPoint[nelems_tail / 2];

			for (int i = 0; i < nelems_tail / 2; i++) {

				DataPoint value  = new DataPoint(rrc_scaledx1_tail[i], rrc_scaledy1_tail[i]);
				DataPoint value1 = new DataPoint(rrc_scaledx2_tail[i], rrc_scaledy2_tail[i]);

				DataPoint valueAssemblyWing   = new DataPoint(rrc_scaledx1[i]*chLnW, rrc_scaledy1[i]*chLnW);
				DataPoint valueAssemblyWing1  = new DataPoint(rrc_scaledx2[i]*chLnW, rrc_scaledy2[i]*chLnW);
				DataPoint valueAssemblyTail   = new DataPoint(rrc_scaledx1_tail[i]*chLnT+wTX, rrc_scaledy1_tail[i]*chLnT+wTY);
				DataPoint valueAssemblyTail1  = new DataPoint(rrc_scaledx2_tail[i]*chLnT+wTX, rrc_scaledy2_tail[i]*chLnT+wTY);

				if(isGround){
					valueAssemblyWing   = new DataPoint(rrc_scaledx1[i]*chLnW, rrc_scaledy1[i]*chLnW+grndH);
					valueAssemblyWing1  = new DataPoint(rrc_scaledx2[i]*chLnW, rrc_scaledy2[i]*chLnW+grndH);
					valueAssemblyTail   = new DataPoint(rrc_scaledx1_tail[i]*chLnT+wTX, rrc_scaledy1_tail[i]*chLnT+wTY+grndH);
					valueAssemblyTail1  = new DataPoint(rrc_scaledx2_tail[i]*chLnT+wTX, rrc_scaledy2_tail[i]*chLnT+wTY+grndH);
				}

				values[i]  = value;
				values1[i] = value1;

				valuesAssemblyWing[i]  = valueAssemblyWing;
				valuesAssemblyWing1[i] = valueAssemblyWing1;
				valuesAssemblyTail[i]  = valueAssemblyTail;
				valuesAssemblyTail1[i] = valueAssemblyTail1;
			}

			series  = new LineGraphSeries<>(values);
			series1 = new LineGraphSeries<>(values1);
			graph.addSeries(series);
			graph.addSeries(series1);
			max_x = 1.0;
			max_y = max(Arrays.stream(rr_tail[1]).max().getAsDouble()/chLnT, Arrays.stream(rr_wing[1]).max().getAsDouble()/chLnW);
			min_y = min(Arrays.stream(rr_tail[1]).min().getAsDouble()/chLnT, Arrays.stream(rr_wing[1]).min().getAsDouble()/chLnW);
			graph.getViewport().setMaxX(max_x);
			graph.getViewport().setMaxY(3 * abs(max_y));
			graph.getViewport().setMinY(-3 * abs(min_y));
			series.setDrawDataPoints(true);
			series.setDataPointsRadius(6);
			series1.setDrawDataPoints(true);
			series1.setDataPointsRadius(6);
			graph.getViewport().setXAxisBoundsManual(true);
			graph.getViewport().setYAxisBoundsManual(true);
			graph.getViewport().setScalableY(false);
			graph.setTitle(getString(R.string.tail) + ": " + nacaDT);
			graph.setTitleTextSize(40);
			series.setColor(Color.rgb(0, 145, 234));
			series1.setColor(Color.rgb(0, 145, 234));

			LineGraphSeries<DataPoint> seriesWing  = new LineGraphSeries<>(valuesAssemblyWing);
			LineGraphSeries<DataPoint> seriesWing1 = new LineGraphSeries<>(valuesAssemblyWing1);
			LineGraphSeries<DataPoint> seriesTail  = new LineGraphSeries<>(valuesAssemblyTail);
			LineGraphSeries<DataPoint> seriesTail1 = new LineGraphSeries<>(valuesAssemblyTail1);
			graph_assembly.addSeries(seriesWing);
			graph_assembly.addSeries(seriesWing1);
			graph_assembly.addSeries(seriesTail);
			graph_assembly.addSeries(seriesTail1);
			max_x = wTX+chLnT;
			if(!isGround){
				max_y = max(Arrays.stream(rr_tail[1]).max().getAsDouble()+wTY, Arrays.stream(rr_wing[1]).max().getAsDouble());
				min_y = min(Arrays.stream(rr_tail[1]).min().getAsDouble()+wTY, Arrays.stream(rr_wing[1]).min().getAsDouble());
			}else{
				max_y = max(Arrays.stream(rr_tail[1]).max().getAsDouble(), Arrays.stream(rr_wing[1]).max().getAsDouble()+wTY)+grndH;
			}
			graph_assembly.getViewport().setMaxX(max_x);
			graph_assembly.getViewport().setMaxY(1.5 * max_y);
			if(!isGround) {
				graph_assembly.getViewport().setMinY(-1.5 * abs(min_y));
			}else{
				graph_assembly.getViewport().setMinY(0);
			}
			seriesWing.setDrawDataPoints(true);
			seriesWing.setDataPointsRadius(6);
			seriesWing1.setDrawDataPoints(true);
			seriesWing1.setDataPointsRadius(6);
			seriesTail.setDrawDataPoints(true);
			seriesTail.setDataPointsRadius(6);
			seriesTail1.setDrawDataPoints(true);
			seriesTail1.setDataPointsRadius(6);
			graph_assembly.getViewport().setXAxisBoundsManual(true);
			graph_assembly.getViewport().setYAxisBoundsManual(true);
			graph_assembly.getViewport().setScalableY(false);
			graph_assembly.setTitle(getString(R.string.assembly));
			graph_assembly.setTitleTextSize(40);
			graph_assembly.getViewport().setScalableY(false);
			seriesWing.setColor(Color.rgb(0, 145, 234));
			seriesWing1.setColor(Color.rgb(0, 145, 234));
			seriesTail.setColor(Color.rgb(0, 145, 234));
			seriesTail1.setColor(Color.rgb(0, 145, 234));

			// The code below refers to 2 tables, which show the y value of the airfoil as a function
			// of the position respectively on the lower and upper surface of the airfoil.

			table  = findViewById(R.id.table_tail);
			table1 = findViewById(R.id.table1_tail);

			row_header  = new TableRow(this);
			row_header1 = new TableRow(this);
			TextView x_cT         = new TextView(this);
			x_cT.setText(R.string.x_c);
			TextView x_c1T        = new TextView(this);
			x_c1T.setText(R.string.x_c);
			TextView yc_textT     = new TextView(this);
			yc_textT.setText(R.string.y_c);
			TextView yc_text1T    = new TextView(this);
			yc_text1T.setText(R.string.y_c);
			row_header.addView(x_cT);
			row_header.addView(yc_textT);
			row_header1.addView(x_c1T);
			row_header1.addView(yc_text1T);
			table.addView(row_header);
			table1.addView(row_header1);

			for (int i = 0; i < nelems_tail / 2; i++) {
				TableRow row  = new TableRow(this);
				TableRow row1 = new TableRow(this);
				TextView lowerSurface = new TextView(this);
				TextView upperSurface = new TextView(this);
				lowerSurface.setText(String.format("%.3f", rrc_scaledx_tail[nelems_tail / 2 - 1 - i]));
				lowerSurface.append("          ");
				upperSurface.setText(String.format("%.3f", rrc_scaledx_tail[nelems_tail / 2 + i]));
				upperSurface.append("          ");
				TextView ycLow = new TextView(this);
				TextView ycUp = new TextView(this);
				ycLow.setText(String.format("%.3f", rrc_scaledy_tail[nelems_tail / 2 - 1 - i]));
				ycUp.setText(String.format("%.3f", rrc_scaledy_tail[nelems_tail / 2 + i]));
				row.addView(lowerSurface);
				row.addView(ycLow);
				row1.addView(upperSurface);
				row1.addView(ycUp);
				table.addView(row);
				table1.addView(row1);
			}
		}
	}
}