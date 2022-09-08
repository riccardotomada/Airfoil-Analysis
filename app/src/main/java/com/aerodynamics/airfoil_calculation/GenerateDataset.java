package com.aerodynamics.airfoil_calculation;

import androidx.appcompat.app.AppCompatActivity;
import androidx.constraintlayout.widget.ConstraintLayout;

import android.app.ProgressDialog;
import android.content.Context;
import android.content.Intent;
import android.graphics.Color;
import android.graphics.PorterDuff;
import android.net.Uri;
import android.os.AsyncTask;
import android.os.Build;
import android.os.Bundle;
import android.os.Environment;
import android.provider.Settings;
import android.view.View;
import android.view.inputmethod.InputMethodManager;
import android.widget.Button;
import android.widget.EditText;
import android.widget.ProgressBar;
import android.widget.Toast;

import com.google.android.material.switchmaterial.SwitchMaterial;
import com.google.firebase.analytics.FirebaseAnalytics;
import com.opencsv.CSVWriter;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.DoubleStream;

public class GenerateDataset extends AppCompatActivity {

	SwitchMaterial naca4digits, naca5digits;
	ConstraintLayout all1, all2, s1, d1, s2, d2, s3, d3, s4, d4, s4_, d4_, s5, d5, s6, d6, s7, d7, s8, d8, s9, d9;
	EditText minMaxC, maxMaxC, minClDes, maxClDes, minMaxcPos4, maxMaxcPos4,minMaxcPos5,maxMaxcPos5,
			 minThick4, maxThick4, minThick5, maxThick5, minAoa4, maxAoa4, minAoa5, maxAoa5, speeD4, speeD5;
	String minfirstdig4, maxfirstdig4, minsecdig4, maxsecdig4, mintrqtdig4, maxtrqtdig4,
		   minfirstdig5, maxfirstdig5, minsectrdig5, maxsectrdig5, minqtcin5, maxqtcin5;
	String nacaid;
	Double minaoa4, maxaoa4, minaoa5, maxaoa5, speed4, speed5;
	Button button;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_generate_dataset);

		FirebaseAnalytics mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.R && false == Environment.isExternalStorageManager()) {
			Uri uri = Uri.parse("package:" + BuildConfig.APPLICATION_ID);
			startActivity(new Intent(Settings.ACTION_MANAGE_APP_ALL_FILES_ACCESS_PERMISSION, uri));
		}

		naca4digits = findViewById(R.id.switch4);
		naca5digits = findViewById(R.id.switch5);

		button = findViewById(R.id.button);
		int defCol = button.getTextColors().getDefaultColor();

		all1 = findViewById(R.id.all1);
		all2 = findViewById(R.id.all2);
		s1   = findViewById(R.id.s1);
		d1   = findViewById(R.id.d1);
		s2   = findViewById(R.id.s2);
		d2   = findViewById(R.id.d2);
		s3   = findViewById(R.id.s3);
		d3   = findViewById(R.id.d3);
		s4   = findViewById(R.id.s4);
		d4   = findViewById(R.id.d4);
		s4_  = findViewById(R.id.s4_);
		d4_  = findViewById(R.id.d4_);
		s5   = findViewById(R.id.s5);
		d5   = findViewById(R.id.d5);
		s6   = findViewById(R.id.s6);
		d6   = findViewById(R.id.d6);
		s7   = findViewById(R.id.s7);
		d7   = findViewById(R.id.d7);
		s8   = findViewById(R.id.s8);
		d8   = findViewById(R.id.d8);
		s9   = findViewById(R.id.s9);
		d9   = findViewById(R.id.d9);

		minMaxC = findViewById(R.id.editMinMaxCamber);
		minMaxC.setHint("0");
		maxMaxC = findViewById(R.id.editMaxMaxCamber);
		maxMaxC.setHint("9");
		minMaxcPos4 = findViewById(R.id.editMinMaxCamberpos);
		minMaxcPos4.setHint("0");
		maxMaxcPos4 = findViewById(R.id.editMaxMaxCamberpos);
		maxMaxcPos4.setHint("9");
		minThick4 = findViewById(R.id.editMinThickness);
		minThick4.setHint("01");
		maxThick4 = findViewById(R.id.editMaxThickness);
		maxThick4.setHint("40");
		minAoa4 = findViewById(R.id.editMinAOA);
		minAoa4.setHint("-8.0");
		maxAoa4 = findViewById(R.id.editMaxAOA);
		maxAoa4.setHint("8.0");
		speeD4 = findViewById(R.id.editSpeed4);
		speeD4.setHint("50.0");

		minClDes = findViewById(R.id.editMinCL);
		minClDes.setHint("1");
		maxClDes = findViewById(R.id.editMaxCL);
		maxClDes.setHint("6");
		minMaxcPos5 = findViewById(R.id.editMinCamberpos);
		minMaxcPos5.setHint("10");
		maxMaxcPos5 = findViewById(R.id.editMaxCamberpos);
		maxMaxcPos5.setHint("50");
		minThick5 = findViewById(R.id.editMinThick);
		minThick5.setHint("01");
		maxThick5 = findViewById(R.id.editMaxThick);
		maxThick5.setHint("30");
		minAoa5 = findViewById(R.id.editMinAOA5);
		minAoa5.setHint("-8.0");
		maxAoa5 = findViewById(R.id.editMaxAOA5);
		maxAoa5.setHint("8.0");
		speeD5 = findViewById(R.id.editSpeed5);
		speeD5.setHint("50.0");

		naca4digits.setOnCheckedChangeListener((buttonView, isChecked) -> {
			if (isChecked) {
				naca5digits.setChecked(false);

				button.setEnabled(true);
				button.getBackground().setColorFilter(null);
				button.setTextColor(defCol);

				all1.setVisibility(View.VISIBLE);
				s1.setVisibility(View.VISIBLE);
				d1.setVisibility(View.VISIBLE);
				s2.setVisibility(View.VISIBLE);
				d2.setVisibility(View.VISIBLE);
				s3.setVisibility(View.VISIBLE);
				d3.setVisibility(View.VISIBLE);
				s4.setVisibility(View.VISIBLE);
				d4.setVisibility(View.VISIBLE);
				s4_.setVisibility(View.VISIBLE);
				d4_.setVisibility(View.VISIBLE);

				all2.setVisibility(View.GONE);
				s5.setVisibility(View.GONE);
				d5.setVisibility(View.GONE);
				s6.setVisibility(View.GONE);
				d6.setVisibility(View.GONE);
				s7.setVisibility(View.GONE);
				d7.setVisibility(View.GONE);
				s8.setVisibility(View.GONE);
				d8.setVisibility(View.GONE);
				s9.setVisibility(View.GONE);
				d9.setVisibility(View.GONE);
			}
		});
		naca5digits.setOnCheckedChangeListener((buttonView, isChecked) -> {
			if (isChecked) {
				naca4digits.setChecked(false);

				button.setEnabled(true);
				button.getBackground().setColorFilter(null);
				button.setTextColor(defCol);

				all1.setVisibility(View.GONE);
				s1.setVisibility(View.GONE);
				d1.setVisibility(View.GONE);
				s2.setVisibility(View.GONE);
				d2.setVisibility(View.GONE);
				s3.setVisibility(View.GONE);
				d3.setVisibility(View.GONE);
				s4.setVisibility(View.GONE);
				d4.setVisibility(View.GONE);
				s4_.setVisibility(View.GONE);
				d4_.setVisibility(View.GONE);

				all2.setVisibility(View.VISIBLE);
				s5.setVisibility(View.VISIBLE);
				d5.setVisibility(View.VISIBLE);
				s6.setVisibility(View.VISIBLE);
				d6.setVisibility(View.VISIBLE);
				s7.setVisibility(View.VISIBLE);
				d7.setVisibility(View.VISIBLE);
				s8.setVisibility(View.VISIBLE);
				d8.setVisibility(View.VISIBLE);
				s9.setVisibility(View.VISIBLE);
				d9.setVisibility(View.VISIBLE);
			}
		});
		if(!naca4digits.isChecked() && !naca5digits.isChecked()){
			button.setEnabled(false);
			button.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
			button.setTextColor(Color.LTGRAY);
		}
	}

	public void onClick(View view) {

		view = this.getCurrentFocus();
		if (view != null) {
			InputMethodManager imm = (InputMethodManager)getSystemService(Context.INPUT_METHOD_SERVICE);
			imm.hideSoftInputFromWindow(view.getWindowToken(), 0);
		}

		if(naca4digits.isChecked()) {

			if (minMaxC.getText().toString().isEmpty()) {
				minfirstdig4 = "0";
			} else
				minfirstdig4 = minMaxC.getText().toString();

			if (maxMaxC.getText().toString().isEmpty()) {
				maxfirstdig4 = "9";
			} else
				maxfirstdig4 = maxMaxC.getText().toString();

			if (minMaxcPos4.getText().toString().isEmpty()) {
				minsecdig4 = "0";
			} else
				minsecdig4 = minMaxcPos4.getText().toString();

			if (maxMaxcPos4.getText().toString().isEmpty()) {
				maxsecdig4 = "9";
			} else
				maxsecdig4 = maxMaxcPos4.getText().toString();

			if (minThick4.getText().toString().isEmpty()) {
				mintrqtdig4 = "01";
			} else
				mintrqtdig4 = minThick4.getText().toString();

			if (maxThick4.getText().toString().isEmpty()) {
				maxtrqtdig4 = "30";
			} else
				maxtrqtdig4 = maxThick4.getText().toString();

			if (minAoa4.getText().toString().isEmpty()) {
				minaoa4 = -8.0;
			} else
				minaoa4 = Double.parseDouble(minAoa4.getText().toString());

			if (maxAoa4.getText().toString().isEmpty()) {
				maxaoa4 = 8.0;
			} else
				maxaoa4 = Double.parseDouble(maxAoa4.getText().toString());
			if (speeD4.getText().toString().isEmpty()){
				speed4 = 50.0;
			} else
				speed4 = Double.parseDouble(speeD4.getText().toString());
		}
		else if(naca5digits.isChecked()){

			if (minClDes.getText().toString().isEmpty()) {
				minfirstdig5 = "1";
			} else
				minfirstdig5 = minClDes.getText().toString();

			if (maxClDes.getText().toString().isEmpty()) {
				maxfirstdig5 = "6";
			} else
				maxfirstdig5 = maxClDes.getText().toString();

			if (minMaxcPos5.getText().toString().isEmpty()) {
				minsectrdig5 = "10";
			} else
				minsectrdig5 = minMaxcPos5.getText().toString();

			if (maxMaxcPos5.getText().toString().isEmpty()) {
				maxsectrdig5 = "50";
			} else
				maxsectrdig5 = maxMaxcPos5.getText().toString();

			if (minThick5.getText().toString().isEmpty()) {
				minqtcin5 = "01";
			} else
				minqtcin5 = minThick5.getText().toString();

			if (maxThick5.getText().toString().isEmpty()) {
				maxqtcin5 = "30";
			} else
				maxqtcin5 = maxThick5.getText().toString();

			if (minAoa5.getText().toString().isEmpty()) {
				minaoa5 = -8.0;
			} else
				minaoa5 = Double.parseDouble(minAoa5.getText().toString());

			if (maxAoa5.getText().toString().isEmpty()) {
				maxaoa5 = 8.0;
			} else
				maxaoa5 = Double.parseDouble(maxAoa5.getText().toString());

			if (speeD5.getText().toString().isEmpty()){
				speed5 = 50.0;
			} else
				speed5 = Double.parseDouble(speeD4.getText().toString());
		}

		if (naca4digits.isChecked()) {

			if (Integer.parseInt(maxfirstdig4) < Integer.parseInt(minfirstdig4)) {
				Toast.makeText(this, R.string.err1dg, Toast.LENGTH_LONG).show();
			} else if (Integer.parseInt(maxsecdig4) < Integer.parseInt(minsecdig4)) {
				Toast.makeText(this, R.string.err2dg, Toast.LENGTH_LONG).show();
			} else if (mintrqtdig4.length() != 2 || maxtrqtdig4.length() != 2) {
				Toast.makeText(this, R.string.err34dg, Toast.LENGTH_LONG).show();
			} else if (minaoa4 < -8.0) {
				Toast.makeText(this, R.string.aoamin, Toast.LENGTH_LONG).show();
			} else if (maxaoa4 > 8.0) {
				Toast.makeText(this, R.string.aoamax, Toast.LENGTH_LONG).show();
			} else if (minaoa4 > maxaoa4) {
				Toast.makeText(this, R.string.aoa, Toast.LENGTH_LONG).show();
			} else if (speed4 <= 0) {
				Toast.makeText(this, R.string.check_speed, Toast.LENGTH_LONG).show();
			} else if (speed4 >= 100) {
				Toast.makeText(this, R.string.check_comp, Toast.LENGTH_LONG).show();
			} else {

				String csv = (Environment.getExternalStorageDirectory().getAbsolutePath() + "/CL_Dataset4digits.csv");

				CSVWriter writer = null;
				try {
					writer = new CSVWriter(new FileWriter(csv));

					List<String[]> data = new ArrayList<String[]>();
					data.add(new String[]{"ID", "AOA", "Cl", "Cm", "Cl Thwaites", "Cd Thwaites"});
					writer.writeAll(data); // data is adding to csv

					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}

				new MyTask().execute(csv);
			}
		}

		if (naca5digits.isChecked()){
			if (Integer.parseInt(maxfirstdig5) < Integer.parseInt(minfirstdig5)) {
				Toast.makeText(this, R.string.err1dg, Toast.LENGTH_LONG).show();
			} else if (minsectrdig5.length() != 2 || maxsectrdig5.length() != 2) {
				Toast.makeText(this, R.string.err23dglen, Toast.LENGTH_LONG).show();
			} else if ((((Integer.parseInt(minsectrdig5) != 10)) && ((Integer.parseInt(minsectrdig5) != 20))
					&& ((Integer.parseInt(minsectrdig5) != 30)) && ((Integer.parseInt(minsectrdig5) != 40))
					&& ((Integer.parseInt(minsectrdig5) != 50))) || (((Integer.parseInt(maxsectrdig5) != 10))
					&& ((Integer.parseInt(maxsectrdig5) != 20)) && ((Integer.parseInt(maxsectrdig5) != 30))
					&& ((Integer.parseInt(maxsectrdig5) != 40)) && ((Integer.parseInt(maxsectrdig5) != 50)))){
				Toast.makeText(this, R.string.err23dgval, Toast.LENGTH_LONG).show();
			} else if (minqtcin5.length() != 2 || maxqtcin5.length() != 2) {
				Toast.makeText(this, R.string.err45dglen, Toast.LENGTH_LONG).show();
			} else if (minaoa5 < -8.0) {
				Toast.makeText(this, R.string.aoamin, Toast.LENGTH_LONG).show();
			} else if (maxaoa5 > 8.0) {
				Toast.makeText(this, R.string.aoamax, Toast.LENGTH_LONG).show();
			} else if (minaoa5 > maxaoa5) {
				Toast.makeText(this, R.string.aoa, Toast.LENGTH_LONG).show();
			} else if (speed5 <= 0) {
				Toast.makeText(this, R.string.check_speed, Toast.LENGTH_LONG).show();
			} else if (speed5 >= 100) {
				Toast.makeText(this, R.string.check_comp, Toast.LENGTH_LONG).show();
			} else {

				String csv = (Environment.getExternalStorageDirectory().getAbsolutePath() + "/CL_Dataset5digits.csv");

				CSVWriter writer = null;
				try {
					writer = new CSVWriter(new FileWriter(csv));

					List<String[]> data = new ArrayList<String[]>();
					data.add(new String[]{"ID", "AOA", "Cl", "Cm", "Cl Thwaites", "Cd Thwaites"});
					writer.writeAll(data); // data is adding to csv

					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}

				new MyTask().execute(csv);
			}
		}
}

public void generateDataset(String nacaid, Double aoa, Double speed, String csv) {

	double spd, pr, den;

	spd = speed;
	pr = 103125; //Pa
	den = 1.225; //kg/m3

	double[] refPoint = {0, 0};

	Airfoil airfoilW = new Airfoil(1, nacaid, 1.0, -Math.toRadians(0.0),
			0.25, refPoint, 50);

	GlobalFunctions myGlobal = new GlobalFunctions();

	double[][] rrW = myGlobal.getCoordinates(airfoilW, false);
	int[][] eeW = myGlobal.connectivityMatrix(rrW[0]);
	int nelems = Array.getLength(eeW[0]);
	Elems[] elemsW = myGlobal.elemsMethod(airfoilW, rrW, eeW, nelems);

	double[] lenW = new double[nelems];
	double[][] tversW = new double[2][nelems];
	double[][] nversW = new double[2][nelems];

	for (int i = 0; i < nelems; i++) {
		lenW[i] = elemsW[i].len;
		tversW[0][i] = elemsW[i].tver[0];
		tversW[1][i] = elemsW[i].tver[1];
		nversW[0][i] = elemsW[i].nver[0];
		nversW[1][i] = elemsW[i].nver[1];
	}

	FreeStream freeStream = new FreeStream(pr, den, spd, Math.toRadians(aoa));

	double[][] A, Au, Av;
	double[][] A_mirr, Au_mirr, Av_mirr;
	double[] b;
	double[] b_mirr;

	A = myGlobal.coefficientMatrix(elemsW, 50);
	b = myGlobal.coefficientArray(freeStream, elemsW, 50);
	Au = myGlobal.onBodyUMatrix(elemsW, 50);
	Av = myGlobal.onBodyVMatrix(elemsW, 50);

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

	for (int i = 0; i < nelems; i++) {
		cPW[i] = 1.0 - vTiW[i] * vTiW[i] / (freeStream.v * freeStream.v);
	}

	double[][] cp_len = new double[2][cPW.length];

	for (int i = 0; i < cPW.length; i++) {
		cp_len[0][i] = cPW[i] * lenW[i];
		cp_len[1][i] = cPW[i] * lenW[i];
	}

	double[][] cp_len_nversW = new double[2][cPW.length];

	for (int i = 0; i < cPW.length; i++) {
		cp_len_nversW[0][i] = cp_len[0][i] * nversW[0][i];
		cp_len_nversW[1][i] = cp_len[1][i] * nversW[1][i];
	}

	double cfx = -DoubleStream.of(cp_len_nversW[0]).sum() / airfoilW.chord;
	double cfy = -DoubleStream.of(cp_len_nversW[1]).sum() / airfoilW.chord;

	double clW = -Math.sin(freeStream.alpha) * cfx + Math.cos(freeStream.alpha) * cfy;

	double CmW = 0;
	double dxW,dyW;
	double x_acW = 0.25*airfoilW.chord;

	for(int i = 0; i<cPW.length;i++){

		dxW = elemsW[i].ver2[0]-elemsW[i].ver1[0];
		dyW = elemsW[i].ver2[1]-elemsW[i].ver1[1];
		CmW = CmW + cPW[i]*(dxW*(elemsW[i].cen[0]-0.25)+dyW*(elemsW[i].cen[1]-refPoint[1]));
	}



	//The following lines of code contain the implementation of the Thwaites method for the
	//viscous flow corrections.

	//Find the stagnation point

	int n_stg = 0;
	double csi_stg = 0.0;
	int i_stg = 0;

	for(int i = 0; i<nelems-1; i++){
		if(vTiW[i]*vTiW[i+1]<=0) {
			n_stg++;
			i_stg = i;
			csi_stg = -vTiW[i]/(vTiW[i+1] - vTiW[i]); // coordinate along the airfoil of the stagnation point
		}
	}

	// arc-length from the stagnation point

	double[] s = new double[nelems];

	s[i_stg+1] = csi_stg * (lenW[i_stg] + lenW[i_stg+1]) * 0.5;
	s[i_stg] = (1.0 - csi_stg) * (lenW[i_stg] + lenW[i_stg+1]) * 0.5;

	for(int i = i_stg+2; i < nelems; i++){
		s[i] = s[i-1] + 0.5 * (lenW[i] + lenW[i-1]);
	}
	for(int i = i_stg-1; i>-1; i--){
		s[i] = s[i+1] + 0.5 * (lenW[i] + lenW[i+1]);
	}

	// Initial condition of the momentum thickness
	double th = Math.acos(nversW[0][i_stg] * nversW[0][i_stg + 1] + nversW[1][i_stg] * nversW[1][i_stg + 1]);
	double r0 = lenW[1] * Math.sin(th) + (lenW[0] + lenW[1] * Math.cos(th)) / Math.tan(th);
	double k = freeStream.v / r0;
	double theta0 = Math.sqrt(0.075 * freeStream.kin_visc / k);  //never used

	// Integrate Thwaites' equation in both directions

	double[] Ue = new double[vTiW.length];
	for(int i = 0; i<vTiW.length; i++){
		if(vTiW[i] > 0)
			Ue[i] = vTiW[i];
		else
			Ue[i] = -vTiW[i];
	}

	double[] theta2Ue6 			= 		new double[nelems];
	double[] theta 				= 		new double[nelems];
	double[] ReTheta			=		new double[nelems];
	double[] lambda 			=		new double[nelems];
	double[] ell 				= 		new double[nelems];
	double[] H					=		new double[nelems];
	double[] delta				= 		new double[nelems];
	double[] cf					= 		new double[nelems];
	double[] UeThetaH1 			= 		new double[nelems]; // for turbulent b.l.
	double[] H1					= 		new double[nelems]; // for turbulent b.l.

	double[] Res				=		new double[Ue.length];

	// Pressure and its derivatives along the surface
	double[] P					=		new double[cPW.length];
	for(int i = 0; i < Ue.length; i++){
		Res[i] = Ue[i] * s[i] / freeStream.kin_visc;
		P[i] = cPW[i] * (0.5 * freeStream.rho * Math.pow(freeStream.v, 2));
	}

	double[] dPdx			= 		new double[nelems];
	double[] dUedx			= 		new double[nelems];

	double Ptot = freeStream.P + 0.5 * freeStream.rho * Math.pow(freeStream.v, 2);

	dPdx[i_stg+1]	 = 	(P[i_stg+1] - Ptot) / (0.5 * csi_stg * (lenW[i_stg] + lenW[i_stg+1]));
	dPdx[i_stg]		 = 	(P[i_stg] - Ptot) / (0.5 * (1.0 - csi_stg) * (lenW[i_stg] + lenW[i_stg+1]));
	dUedx[i_stg+1]	 = 	(Ue[i_stg+1] - 0.0) / (0.5 * csi_stg * (lenW[i_stg] + lenW[i_stg+1]));
	dUedx[i_stg]	 = 	(Ue[i_stg] -   0.0) / (0.5 * (1.0 - csi_stg) * (lenW[i_stg] + lenW[i_stg+1]));

	// Initial conditions

	theta2Ue6[i_stg + 1] = 0.45 * freeStream.kin_visc * Math.pow(Ue[i_stg+1], 5.0) *
			(0.25 * csi_stg * (lenW[i_stg] + lenW[i_stg + 1]));
	theta2Ue6[i_stg]   = 0.45 * freeStream.kin_visc * Math.pow(Ue[i_stg], 5.0)*
			(0.25 * (1.0 - csi_stg) * (lenW[i_stg] + lenW[i_stg + 1]));
	theta[i_stg + 1]   = Math.sqrt(theta2Ue6[i_stg + 1] / Math.pow(Ue[i_stg + 1], 6.0));
	theta[i_stg]       = Math.sqrt(theta2Ue6[i_stg] / Math.pow(Ue[i_stg], 6.0));
	ReTheta[i_stg+1]   = Ue[i_stg+1] * theta[i_stg+1] / freeStream.kin_visc;
	ReTheta[i_stg]     = Ue[i_stg] * theta[i_stg] / freeStream.kin_visc;
	lambda[i_stg+1]    = Math.pow(theta[i_stg + 1], 2.0) * dUedx[i_stg + 1] / freeStream.kin_visc;
	lambda[i_stg]      = Math.pow(theta[i_stg], 2.0) * dUedx[i_stg] / freeStream.kin_visc;
	ell[i_stg+1] 	   = myGlobal.thwaites_ell(lambda[i_stg + 1]);
	ell[i_stg] 	       = myGlobal.thwaites_ell(lambda[i_stg]);
	H[i_stg+1]		   = myGlobal.thwaites_H(lambda[i_stg + 1]);
	H[i_stg]		   = myGlobal.thwaites_H(lambda[i_stg]);
	delta[i_stg + 1]   = theta[i_stg + 1] * H[i_stg + 1];
	delta[i_stg]       = theta[i_stg] * H[i_stg];
	cf[i_stg + 1]      = 2.0 * ell[i_stg + 1] / ReTheta[i_stg + 1];
	cf[i_stg]		   = 2.0 * ell[i_stg] / ReTheta[i_stg];

	// Integrate equation

	String regime = "laminar";
	String transition = "no";

	for(int i = i_stg + 2; i<nelems; i++){

		if(i < nelems - 1){
			dPdx[i]  = (P[i + 1] - P[i - 1]) / (0.5 * (lenW[i + 1] + lenW[i - 1]) + lenW[i]);
			dUedx[i] = (Ue[i + 1] - Ue[i - 1]) / (0.5 * ( lenW[i + 1] + lenW[i - 1] ) + lenW[i]);
		}
		else{
			dPdx[i]  = (P[i] - P[i - 1]) / (0.5 * (lenW[i] + lenW[i - 1]));
			dUedx[i] = (Ue[i] - Ue[i - 1]) / ( 0.5 * ( lenW[i] + lenW[i - 1] ));
		}

		if(regime.equals("laminar")) {
			theta2Ue6[i] = theta2Ue6[i - 1] + 0.45 * freeStream.kin_visc * (Math.pow(Ue[i], 5.0) +
					Math.pow(Ue[i - 1], 5.0)) * (0.25 * (lenW[i] + lenW[i - 1]));

			theta[i]     = Math.sqrt(theta2Ue6[i] / Math.pow(Ue[i], 6.0));
			ReTheta[i]   = Ue[i] * theta[i] / freeStream.kin_visc;
			lambda[i]    = Math.pow(theta[i], 2.0) * dUedx[i] / freeStream.kin_visc;
			ell[i]       = myGlobal.thwaites_ell(lambda[i]);
			H[i]         = myGlobal.thwaites_H(lambda[i]);
			delta[i]     = theta[i] * H[i];
			cf[i]        = 2.0 * ell[i] / ReTheta[i];

			// Michel's criterion for transition
			if (ReTheta[i] > 1.174 * (1.0 + 22400.0 / Res[i]) * Math.pow(Res[i], 0.46)) {
				regime = "turbulent";
				transition = "transition";
				H[i] = 1.35;
				delta[i] = theta[i] * H[i];
				H1[i] = myGlobal.head_HtoH1(H[i]);
				UeThetaH1[i] = Ue[i] * theta[i] * H1[i];
			}

		} else {

			double dx1 = (0.5 * (lenW[i] + lenW[i - 1]));
			UeThetaH1[i] = UeThetaH1[i - 1] + dx1 * Ue[i - 1] * 0.0306 / Math.pow(H1[i - 1]
					- 3.0, 0.6169);
			theta[i]     = theta[i - 1] + dx1 * (0.5 * cf[i - 1] - dUedx[i - 1] / Ue[i - 1]
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

	for(int i = i_stg - 1; i > -1; i--){
		if(i > 0){
			dPdx[i]  = (P[i - 1] - P[i + 1]) / (0.5 * (lenW[i + 1] + lenW[i - 1]) + lenW[i]);
			dUedx[i] = (Ue[i - 1] - Ue[i + 1]) / (0.5 * (lenW[i+1] + lenW[i - 1]) + lenW[i]);
		}
		else{
			dPdx[i]  = (P[i] - P[i+1]) / (0.5 * (lenW[i] + lenW[i+1]));
			dUedx[i] = (Ue[i] - Ue[i+1]) / (0.5*(lenW[i] + lenW[i+1]));
		}

		if(regime.equals("laminar")){
			theta2Ue6[i] = theta2Ue6[i + 1] + 0.45 * freeStream.kin_visc * (Math.pow(Ue[i], 5.0) +
					Math.pow(Ue[i + 1], 5.0)) * (0.25 * (lenW[i] + lenW[i + 1]));
			theta[i]     = Math.sqrt(theta2Ue6[i] / Math.pow(Ue[i], 6.0));
			ReTheta[i]   = Ue[i] * theta[i] / freeStream.kin_visc;
			lambda[i]    = Math.pow(theta[i], 2.0) * dUedx[i] / freeStream.kin_visc;
			ell[i]       = myGlobal.thwaites_ell(lambda[i]);
			H[i]         = myGlobal.thwaites_H(lambda[i]);
			delta[i]     = theta[i] * H[i];
			cf[i]        = 2.0 * ell[i] / ReTheta[i];

			// Michel's criterion for transition
			if (ReTheta[i] > 1.174 * (1.0 + 22400.0 / Res[i]) * Math.pow(Res[i], 0.46)) {
				regime        = "turbulent";
				transition    = "transition";
				H[i]          = 1.35;
				delta[i]      = theta[i] * H[i];
				H1[i]         = myGlobal.head_HtoH1(H[i]);
				UeThetaH1[i]  = Ue[i] * theta[i] * H1[i];
			}

		}  else {

			double dx1 = (0.5 * (lenW[i] + lenW[i + 1]));
			UeThetaH1[i]  = UeThetaH1[i + 1] + dx1 * Ue[i + 1] * 0.0306 / Math.pow(H1[i + 1]
					- 3.0, 0.6169);
			theta[i]      = theta[i + 1] + dx1 * (0.5 * cf[i + 1] - dUedx[i + 1] / Ue[i + 1]
					* (2.0 + H[i + 1]) * theta[i + 1]);
			H1[i]         = Math.max(UeThetaH1[i] / (Ue[i] * theta[i]), 3.0 + 0.001);
			H[i]          = myGlobal.head_H1toH(H1[i]);
			delta[i]      = H[i] * theta[i];

			ReTheta[i]    = Ue[i] * theta[i] / freeStream.kin_visc;
			cf[i]         = myGlobal.ludweig_tillman_cf(H[i], ReTheta[i]);

			transition    = "occurred";
		}
	}

	// Drag

	double[] tauW = new double[cf.length];

	for(int i = 0; i < cf.length; i++){
		tauW[i] = 0.5 * freeStream.rho * Math.pow(freeStream.v, 2) * cf[i];
	}

	double[][] temp = new double[2][tauW.length];
	double[][] dF_visc = new double[2][tauW.length];

	for(int i = 0; i < tauW.length; i++){
		temp[0][i] = tauW[i] * lenW[i] * vTiW[i] / Ue[i];
		temp[1][i] = tauW[i] * lenW[i] * vTiW[i] / Ue[i];
	}

	for(int i = 0; i <tauW.length; i++){
		dF_visc[0][i] = temp[0][i] * tversW[0][i];
		dF_visc[1][i] = temp[1][i] * tversW[1][i];
	}
	double[] F_visc = new double[2];

	F_visc[0] = DoubleStream.of(dF_visc[0]).sum() / airfoilW.chord;
	F_visc[1] = DoubleStream.of(dF_visc[1]).sum() / airfoilW.chord;

	double L_visc = -F_visc[0] * Math.sin(freeStream.alpha) + F_visc[1] * Math.cos(freeStream.alpha);
	double D_visc = F_visc[0] * Math.cos(freeStream.alpha) + F_visc[1] * Math.sin(freeStream.alpha);

	double[][] temp0 = new double[2][tauW.length];
	double[][] dF_Pres = new double[2][tauW.length];

	for(int i = 0; i<tauW.length;i++){
		temp0[0][i] = -P[i] * lenW[i];
		temp0[1][i] = -P[i] * lenW[i];
	}
	for(int i = 0; i<tauW.length;i++){
		dF_Pres[0][i] = temp0[0][i] * nversW[0][i];
		dF_Pres[1][i] = temp0[1][i] * nversW[1][i];
	}
	double[] F_pres = new double[2];

	F_pres[0] = DoubleStream.of(dF_Pres[0]).sum() / airfoilW.chord;
	F_pres[1] = DoubleStream.of(dF_Pres[1]).sum() / airfoilW.chord;

	double L_pres = -F_pres[0] * Math.sin(freeStream.alpha) + F_pres[1] * Math.cos(freeStream.alpha);
	double D_pres = F_pres[0] * Math.cos(freeStream.alpha) + F_pres[1] * Math.sin(freeStream.alpha);

	double L = L_pres + L_visc;
	double D = D_pres + D_visc;

	double cL = L / (0.5 * freeStream.rho * Math.pow(freeStream.v, 2));
	double cD = D / (0.5 * freeStream.rho * Math.pow(freeStream.v, 2));


	CSVWriter writer = null;
	try {
		writer = new CSVWriter(new FileWriter(csv, true));

		List<String[]> data = new ArrayList<String[]>();
		data.add(new String[]{nacaid, Double.toString(aoa), Double.toString(clW), Double.toString(CmW), Double.toString(cL), Double.toString(cD)});
		writer.writeAll(data); // data is adding to csv
		writer.close();
	} catch (IOException e) {
		e.printStackTrace();
	}
}

private class MyTask extends AsyncTask<String, Integer, String>{

	private ProgressDialog progressDialog;

	@Override
	protected void onPreExecute(){
		super.onPreExecute();

		progressDialog = new ProgressDialog(GenerateDataset.this, R.style.AlertDialogStyle);
		progressDialog.setIndeterminate(false);
		progressDialog.setProgressStyle(ProgressDialog.STYLE_HORIZONTAL);
		progressDialog.setTitle(R.string.generdata);
		progressDialog.setMessage(getResources().getString(R.string.waiting));

		progressDialog.show();

		//Display a progress bar
	}
	// Background thread
	@Override
	protected String doInBackground(String... params){
		// get the values from params, which is an array

		String csv = params[0];

		if (naca4digits.isChecked()){

			int loop1 = Integer.parseInt(maxfirstdig4) - Integer.parseInt(minfirstdig4) + 1;
			int loop2 = Integer.parseInt(maxsecdig4) - Integer.parseInt(minsecdig4) + 1;
			int loop3 = Integer.parseInt(maxtrqtdig4) - Integer.parseInt(mintrqtdig4) + 1;
			double loop4 = (maxaoa4 - minaoa4) * 10 + 1;

			int total = loop1 * loop2 * loop3 * (int)loop4;
			int count = 1;

			for (int dg1 = Integer.parseInt(minfirstdig4); dg1 <= Integer.parseInt(maxfirstdig4); dg1++) {
				for (int dg2 = Integer.parseInt(minsecdig4); dg2 <= Integer.parseInt(maxsecdig4); dg2++) {
					for (int dg34 = Integer.parseInt(mintrqtdig4); dg34 <= Integer.parseInt(maxtrqtdig4); dg34++) {
						for (double aoa = minaoa4; aoa <= maxaoa4+0.0001; aoa += 0.1) {
							if (dg34 < 10) {
								nacaid = "NACA" + dg1 + dg2 + "0" + dg34;
							} else {
								nacaid = "NACA" + dg1 + dg2 + dg34;
							}
							generateDataset(nacaid, aoa, speed4, csv);
							count++;
							publishProgress((int) ((count / (float) total) * 100));

							if (isCancelled()) {
								break;
							}
						}
					}
				}
			}
		}
		else if (naca5digits.isChecked()){
			int loop1 = Integer.parseInt(maxfirstdig5) - Integer.parseInt(minfirstdig5) + 1;
			int loop2 = (Integer.parseInt(maxsectrdig5) - Integer.parseInt(minsectrdig5)) / 10 + 1;
			int loop3 = Integer.parseInt(maxqtcin5) - Integer.parseInt(minqtcin5) + 1;
			double loop4 = (maxaoa5 - minaoa5) * 10 + 1;

			int total = loop1 * loop2 * loop3 * (int)loop4;
			int count = 1;

			for (int dg1 = Integer.parseInt(minfirstdig5); dg1 <= Integer.parseInt(maxfirstdig5); dg1++) {
				for (int dg23 = Integer.parseInt(minsectrdig5); dg23 <= Integer.parseInt(maxsectrdig5); dg23+=10) {
					for (int dg45 = Integer.parseInt(minqtcin5); dg45 <= Integer.parseInt(maxqtcin5); dg45++) {
						for (double aoa = minaoa5; aoa <= maxaoa5+0.0001; aoa += 0.1) {

							if (dg45 < 10) {
								nacaid = "NACA" + dg1 + dg23 + "0" + dg45;
							} else {
								nacaid = "NACA" + dg1 + dg23 + dg45;
							}
							generateDataset(nacaid, aoa, speed5, csv);
							count++;
							publishProgress((int) ((count / (float) total) * 100));

							if (isCancelled()) {
								break;
							}
						}
					}
				}
			}
		}

		return csv;
	}
	//Called from background thread but runs in UI
	@Override
	protected void onProgressUpdate(Integer... values){
		super.onProgressUpdate();
		progressDialog.setProgress(values[0]);

		//update progress bar
	}
	//This runs in UI when background thread finishes
	@Override
	protected void onPostExecute(String result){
		super.onPostExecute(result);

		//Hide progress bar

		progressDialog.dismiss();
		Toast.makeText(GenerateDataset.this, getResources().getString(R.string.generated) + " " + result, Toast.LENGTH_LONG).show();
	}

}

}

