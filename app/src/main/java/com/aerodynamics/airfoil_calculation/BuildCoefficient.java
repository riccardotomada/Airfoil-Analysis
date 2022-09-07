package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//This class contains all the Java code about the "CP, CL, CD" screen   					      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import static java.lang.Math.abs;

import androidx.appcompat.app.AppCompatActivity;
import androidx.constraintlayout.widget.ConstraintLayout;

import android.graphics.Color;
import android.os.AsyncTask;
import android.os.Bundle;
import android.view.View;
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

public class BuildCoefficient extends AppCompatActivity {

	String          nacaDW, nacaDT;
	double          chLnW, chLnT, staggW, staggT, grndH, wTX, wTY, pr, den, aOA, spd;
	int             nPan;
	boolean         isGround, isWT;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_build_coefficient);
		FirebaseAnalytics mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

			// Receiving the values typed by the user in the MainActivity

			nacaDW   = "NACA";
			nacaDW   = nacaDW + getIntent().getStringExtra("NACA_DW");
			chLnW    = getIntent().getDoubleExtra("ChordW", 1.0);
			staggW   = getIntent().getDoubleExtra("StaggerW", 0.0);
			nPan     = getIntent().getIntExtra("nPan", 20);
			nacaDT   = "NACA";
			nacaDT   = nacaDT + getIntent().getStringExtra("NACA_DT");
			chLnT    = getIntent().getDoubleExtra("ChordT", 1.0);
			staggT   = getIntent().getDoubleExtra("StaggerT", 0.0);
			isGround = getIntent().getBooleanExtra("checkGround", false);
			grndH    = getIntent().getDoubleExtra("GroundH", 1.0);
			isWT     = getIntent().getBooleanExtra("checkTail", false);
			wTX      = getIntent().getDoubleExtra("wTX",1.5);
			wTY      = getIntent().getDoubleExtra("wTY",0.2);
			spd      = getIntent().getDoubleExtra("Speed", 30);
			aOA      = getIntent().getDoubleExtra("AngleOA", 0.0);

			pr  = 103125; //Pa
			den = 1.225; //kg/m3.

			double[] refPoint = {0,0};
			if(isGround)
				refPoint[1] = grndH;

			double[] refPointT = {wTX, wTY};

			if(isGround)
				refPointT[1] = grndH + wTY;

			// Now a FreeStream type object, freeStream, is initialized, using the values typed by the user
			// in the MainActivity screen. The same is done with an Airfoil type object, named airfoil.
			FreeStream freeStream = new FreeStream(pr, den, spd, Math.toRadians(aOA));

			Airfoil airfoilW = new Airfoil(1, nacaDW, chLnW, -Math.toRadians(staggW),
					0.25,refPoint, nPan);

			Airfoil airfoilT = new Airfoil(2, nacaDT, chLnT, -Math.toRadians(staggT),
					0.25,refPointT, nPan);

			// Since many functions are meant to be used many times and shared with other activities,
			// I decided to put all of them inside a new class, named GlobalFunctions. In order to use
			// them a GlobalFunctions object, myGlobal, initialized as follows, will be called every time
			// it will be necessary.

			GlobalFunctions myGlobal = new GlobalFunctions();

			double[][] rrW         = myGlobal.getCoordinates(airfoilW, false);
			int[][]    eeW         = myGlobal.connectivityMatrix(rrW[0]);
			int        nelems      = Array.getLength(eeW[0]);
			Elems[]    elemsW      = myGlobal.elemsMethod(airfoilW,rrW,eeW,nelems);
			Elems[]    elemsT      = elemsW;
			Elems[]    elemsW_mirr = elemsW;
			Elems[]    elemsT_mirr = elemsW;

			double[][] rrcW        = new double[2][nelems];
			double[]   lenW        = new double[nelems];
			double[][] tversW      = new double[2][nelems];
			double[][] nversW      = new double[2][nelems];
			double[][] rrcT        = new double[2][nelems];
			double[]   lenT        = new double[nelems];
			double[][] tversT      = new double[2][nelems];
			double[][] nversT      = new double[2][nelems];


			for(int i = 0; i<nelems; i++){
				rrcW[0][i]   = elemsW[i].cen[0];
				rrcW[1][i]   = elemsW[i].cen[1];
				lenW[i]      = elemsW[i].len;
				tversW[0][i] = elemsW[i].tver[0];
				tversW[1][i] = elemsW[i].tver[1];
				nversW[0][i] = elemsW[i].nver[0];
				nversW[1][i] = elemsW[i].nver[1];
			}

			if(isWT){
				double[][] rrT    = myGlobal.getCoordinates(airfoilT, false);
				int[][]    eeT    = myGlobal.connectivityMatrix(rrT[0]);
				elemsT = myGlobal.elemsMethod(airfoilT,rrT,eeT,nelems);

				rrcT   = new double[2][nelems];
				lenT   = new double[nelems];
				tversT = new double[2][nelems];
				nversT = new double[2][nelems];

				for(int i = 0; i<nelems; i++){
					rrcT[0][i]   = elemsT[i].cen[0];
					rrcT[1][i]   = elemsT[i].cen[1];
					lenT[i]      = elemsT[i].len;
					tversT[0][i] = elemsT[i].tver[0];
					tversT[1][i] = elemsT[i].tver[1];
					nversT[0][i] = elemsT[i].nver[0];
					nversT[1][i] = elemsT[i].nver[1];
				}
			}

			if(isGround){
				double[][] rrW_mirr         = myGlobal.getCoordinates(airfoilW, true);
				int[][]    eeW_mirr         = myGlobal.connectivityMatrix(rrW_mirr[0]);
				elemsW_mirr      = myGlobal.elemsMethod(airfoilW,rrW_mirr,eeW_mirr,nelems);

				if(isWT){
					double[][] rrT_mirr         = myGlobal.getCoordinates(airfoilT, true);
					int[][]    eeT_mirr         = myGlobal.connectivityMatrix(rrT_mirr[0]);
					elemsT_mirr      = myGlobal.elemsMethod(airfoilT,rrT_mirr,eeT_mirr,nelems);
				}

			}

			double[][]  A, Au, Av;
			double[]    b;

			if(!isGround && !isWT){
				A  = myGlobal.coefficientMatrix(elemsW, nPan);
				b  = myGlobal.coefficientArray(freeStream, elemsW, nPan);
				Au = myGlobal.onBodyUMatrix(elemsW, nPan);
				Av = myGlobal.onBodyVMatrix(elemsW, nPan);
			}else if(!isGround && isWT){
				A  = myGlobal.coefficientMatrixWT(elemsW, elemsT, nPan);
				b  = myGlobal.coefficientArrayWT(freeStream, elemsW, elemsT, nPan);
				Au = myGlobal.onBodyUMatrixWT(elemsW, elemsT, nPan);
				Av = myGlobal.onBodyVMatrixWT(elemsW, elemsT, nPan);
			}else if(isGround && !isWT){
				A  = myGlobal.coefficientMatrixG(elemsW, elemsW_mirr, nPan);
				b  = myGlobal.coefficientArray(freeStream, elemsW, nPan);
				Au = myGlobal.onBodyUMatrixG(elemsW, elemsW_mirr, nPan);
				Av = myGlobal.onBodyVMatrixG(elemsW, elemsW_mirr, nPan);
			}else{
				A  = myGlobal.coefficientMatrixGWT(elemsW, elemsT, elemsW_mirr, elemsT_mirr, nPan);
				b  = myGlobal.coefficientArrayWT(freeStream, elemsW, elemsT, nPan);
				Au = myGlobal.onBodyUMatrixGWT(elemsW, elemsT, elemsW_mirr, elemsT_mirr, nPan);
				Av = myGlobal.onBodyVMatrixGWT(elemsW, elemsT, elemsW_mirr, elemsT_mirr, nPan);
			}
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

			RealMatrix uvec, uvecTi, uvecNi;

			if(isWT) {
				uvec = new Array2DRowRealMatrix(2 * nelems, 2);
			}else{
				uvec = new Array2DRowRealMatrix(nelems, 2);
			}

			u = AuCoeff.operate(x);
			v = AvCoeff.operate(x);
			u = u.mapAdd(freeStream.v*Math.cos(freeStream.alpha));
			v = v.mapAdd(freeStream.v*Math.sin(freeStream.alpha));

			uvec.setColumnVector(0,u);
			uvec.setColumnVector(1,v);

			uvecTi = uvec.transpose().copy();
			uvecNi = uvec.transpose().copy();

			double[][] uvec_vTiW = new double[2][nelems];
			double[][] uvec_vNiW = new double[2][nelems];

			double[] vTiW = new double[nelems];
			double[] vNiW = new double[nelems];
			double[] vNi_checkW = new double[nelems];


			for (int j = 0; j < nelems; j++) {
				for (int i = 0; i< 2; i++){
					uvecTi.multiplyEntry(i, j, tversW[i][j]);
					uvecNi.multiplyEntry(i, j, nversW[i][j]);
					uvec_vTiW[i][j] = uvecTi.getEntry(i, j);
					uvec_vNiW[i][j] = uvecNi.getEntry(i, j);
				}
				vTiW[j] = uvec_vTiW[0][j] + uvec_vTiW[1][j];
				vNiW[j] = uvec_vNiW[0][j] + uvec_vNiW[1][j];
			}


			for(int i = 0; i<nelems; i++){
				if(vNiW[i]<0)
					vNi_checkW[i] = -vNiW[i];
				else
					vNi_checkW[i] = vNiW[i];

				if(vNi_checkW[i] > Math.pow(10,-6))
					Toast.makeText(BuildCoefficient.this, "Something is wrong", Toast.LENGTH_LONG).show();
			}

			double[] cPW = new double[nelems];
			double[] cPT = new double[nelems];
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

			double cl = -Math.sin(freeStream.alpha)*cfx + Math.cos(freeStream.alpha)*cfy;
			double clT = cl;

			double[] vTiT = vTiW;

			if(isWT){
				double[][] uvec_vTiT = new double[2][nelems];
				double[][] uvec_vNiT = new double[2][nelems];

				vTiT       = new double[nelems];
				double[] vNiT       = new double[nelems];
				double[] vNi_checkT = new double[nelems];


				for (int j = nelems; j < 2*nelems; j++) {
					for (int i = 0; i< 2; i++){
						uvecTi.multiplyEntry(i, j, tversT[i][j-nelems]);
						uvecNi.multiplyEntry(i, j, nversT[i][j-nelems]);
						uvec_vTiT[i][j-nelems] = uvecTi.getEntry(i, j);
						uvec_vNiT[i][j-nelems] = uvecNi.getEntry(i, j);
					}
					vTiT[j-nelems] = uvec_vTiT[0][j-nelems] + uvec_vTiT[1][j-nelems];
					vNiT[j-nelems] = uvec_vNiT[0][j-nelems] + uvec_vNiT[1][j-nelems];
				}

				for(int i = 0; i<nelems; i++){
					if(vNiT[i]<0)
						vNi_checkT[i] = -vNiT[i];
					else
						vNi_checkT[i] = vNiT[i];

					if(vNi_checkT[i] > Math.pow(10,-6))
						Toast.makeText(BuildCoefficient.this, "Something is wrong", Toast.LENGTH_LONG).show();
				}

				for(int i = 0; i<nelems;i++){
					cPT[i] = 1.0 - vTiT[i]*vTiT[i]/(freeStream.v*freeStream.v);
				}

				double[][] cp_lenT = new double[2][cPT.length];

				for(int i = 0;  i<cPT.length; i++){
					cp_lenT[0][i] = cPT[i]*lenT[i];
					cp_lenT[1][i] = cPT[i]*lenT[i];
				}

				double[][] cp_len_nversT = new double[2][cPT.length];

				for(int i = 0;  i<cPT.length; i++){
					cp_len_nversT[0][i] = cp_lenT[0][i]*nversT[0][i];
					cp_len_nversT[1][i] = cp_lenT[1][i]*nversT[1][i];
				}

				double cfxT = -DoubleStream.of(cp_len_nversT[0]).sum()/airfoilT.chord;
				double cfyT = -DoubleStream.of(cp_len_nversT[1]).sum()/airfoilT.chord;

				clT  = -Math.sin(freeStream.alpha)*cfxT + Math.cos(freeStream.alpha)*cfyT;
			}

			double sumW = DoubleStream.of(lenW).sum();
			double GammaW;
			if(isWT) {
				GammaW =x.getEntry(2*nelems) * sumW;
			}
			else{
				GammaW =x.getEntry(nelems) * sumW;
			}

			double liftKJW = -freeStream.rho*freeStream.v*GammaW;
			double liftKJT = 0;

			if(isWT){
				double sumT = DoubleStream.of(lenT).sum();
				double GammaT = x.getEntry(2*nelems+1)*sumT;
				liftKJT = -freeStream.rho*freeStream.v*GammaT;
			}


			double CmW = 0;
			double dxW,dyW;
			double x_acW = 0.25*airfoilW.chord;

			for(int i = 0; i<cPW.length;i++){

				dxW = elemsW[i].ver2[0]-elemsW[i].ver1[0];
				dyW = elemsW[i].ver2[1]-elemsW[i].ver1[1];
				CmW = CmW + cPW[i]*(dxW*(elemsW[i].cen[0]-0.25)+dyW*(elemsW[i].cen[1]-refPoint[1]));
			}

			double CmT = 0;
			if(isWT){

				double dxT,dyT;
				double x_acT = 0.25*airfoilT.chord;

				for(int i = 0; i<cPT.length;i++){

					dxT = elemsT[i].ver2[0]-elemsT[i].ver1[0];
					dyT = elemsT[i].ver2[1]-elemsT[i].ver1[1];
					CmT = CmT + cPT[i]*(dxT*(elemsT[i].cen[0]-0.25)+dyT*(elemsT[i].cen[1]-refPoint[1]));
				}
			}

			//CL calculus using the Kutta-Joukowski theorem.
			double cLKJW = liftKJW / (0.5*freeStream.rho*Math.pow(freeStream.v,2)*airfoilW.chord);
			double clKJT = 0;

			if(isWT){
				clKJT = liftKJT / (0.5*freeStream.rho*Math.pow(freeStream.v,2)*airfoilT.chord);
			}


			double[] rrc_scaledxW = new double[nelems];
			for(int i = 0; i <nelems; i++){
				rrc_scaledxW[i] = rrcW[0][i]/chLnW;
			}

			double[] rrc_scaledyW = new double[nelems];
			for(int i = 0; i <nelems; i++){
				rrc_scaledyW[i] = rrcW[1][i]/chLnW;
			}

			double[] rrc_scaledxT = new double[nelems];
			double[] rrc_scaledyT = new double[nelems];
			double[] rrc_scaledx1W = new double[nelems/2];
			double[] rrc_scaledx2W = new double[nelems/2];
			double[] rrc_scaledy1W = new double[nelems/2];
			double[] rrc_scaledy2W = new double[nelems/2];
			double[] cP1W = new double[nelems/2];
			double[] cP2W = new double[nelems/2];
			double[] rrc_scaledx1T = new double[nelems/2];
			double[] rrc_scaledx2T = new double[nelems/2];
			double[] rrc_scaledy1T = new double[nelems/2];
			double[] rrc_scaledy2T = new double[nelems/2];
			double[] cP1T = new double[nelems/2];
			double[] cP2T = new double[nelems/2];
			System.arraycopy(rrc_scaledxW, 0, rrc_scaledx1W, 0, nelems / 2);
			if (nelems - nelems / 2 >= 0)
				System.arraycopy(rrc_scaledxW, nelems / 2, rrc_scaledx2W, 0, nelems - nelems / 2);
			if (nelems / 2 >= 0) System.arraycopy(rrc_scaledyW, 0, rrc_scaledy1W, 0, nelems / 2);
			if (nelems - nelems / 2 >= 0)
				System.arraycopy(rrc_scaledyW, nelems / 2, rrc_scaledy2W, 0, nelems - nelems / 2);
			if (nelems / 2 >= 0) System.arraycopy(cPW, 0, cP1W, 0, nelems / 2);
			if (nelems - nelems / 2 >= 0)
				System.arraycopy(cPW, nelems / 2, cP2W, 0, nelems - nelems / 2);

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
					if(rrc_scaledx1W[i]>rrc_scaledx1W[i+1]){
						temp1 = rrc_scaledx1W[i];
						temp2 = rrc_scaledy1W[i];
						temp3 = cP1W[i];
						rrc_scaledx1W[i] = rrc_scaledx1W[i+1];
						rrc_scaledy1W[i] = rrc_scaledy1W[i+1];
						cP1W[i] = cP1W[i+1];
						rrc_scaledx1W[i+1] = temp1;
						rrc_scaledy1W[i+1] = temp2;
						cP1W[i+1] = temp3;
						swapped = true;
					}
					if(rrc_scaledx2W[i]>rrc_scaledx2W[i+1]){
						temp1 = rrc_scaledx2W[i];
						temp2 = rrc_scaledy2W[i];
						temp3 = cP2W[i];
						rrc_scaledx2W[i] = rrc_scaledx2W[i+1];
						rrc_scaledy2W[i] = rrc_scaledy2W[i+1];
						cP2W[i] = cP2W[i+1];
						rrc_scaledx2W[i+1] = temp1;
						rrc_scaledy2W[i+1] = temp2;
						cP2W[i+1] = temp3;
						swapped = true;
					}
				}
			}

			if(isWT){
				for(int i = 0; i <nelems; i++){
					rrc_scaledxT[i] = rrcT[0][i]/chLnT;
				}

				for(int i = 0; i <nelems; i++){
					rrc_scaledyT[i] = rrcT[1][i]/chLnT;
				}

				cP1T = new double[nelems/2];
				cP2T = new double[nelems/2];
				System.arraycopy(rrc_scaledxT, 0, rrc_scaledx1T, 0, nelems / 2);
				if (nelems - nelems / 2 >= 0)
					System.arraycopy(rrc_scaledxT, nelems / 2, rrc_scaledx2T, 0, nelems - nelems / 2);
				if (nelems / 2 >= 0)
					System.arraycopy(rrc_scaledyT, 0, rrc_scaledy1T, 0, nelems / 2);
				if (nelems - nelems / 2 >= 0)
					System.arraycopy(rrc_scaledyT, nelems / 2, rrc_scaledy2T, 0, nelems - nelems / 2);
				if (nelems / 2 >= 0) System.arraycopy(cPT, 0, cP1T, 0, nelems / 2);
				if (nelems - nelems / 2 >= 0)
					System.arraycopy(cPT, nelems / 2, cP2T, 0, nelems - nelems / 2);

				swapped = true;
				j = 0;

				while (swapped){
					swapped = false;
					j++;
					for(int i = 0; i<nelems/2 - j; i++){
						if(rrc_scaledx1T[i]>rrc_scaledx1T[i+1]){
							temp1 = rrc_scaledx1T[i];
							temp2 = rrc_scaledy1T[i];
							temp3 = cP1T[i];
							rrc_scaledx1T[i] = rrc_scaledx1T[i+1];
							rrc_scaledy1T[i] = rrc_scaledy1T[i+1];
							cP1T[i] = cP1T[i+1];
							rrc_scaledx1T[i+1] = temp1;
							rrc_scaledy1T[i+1] = temp2;
							cP1T[i+1] = temp3;
							swapped = true;
						}
						if(rrc_scaledx2T[i]>rrc_scaledx2T[i+1]){
							temp1 = rrc_scaledx2T[i];
							temp2 = rrc_scaledy2T[i];
							temp3 = cP2T[i];
							rrc_scaledx2T[i] = rrc_scaledx2T[i+1];
							rrc_scaledy2T[i] = rrc_scaledy2T[i+1];
							cP2T[i] = cP2T[i+1];
							rrc_scaledx2T[i+1] = temp1;
							rrc_scaledy2T[i+1] = temp2;
							cP2T[i+1] = temp3;
							swapped = true;
						}
					}
				}
			}

			// The following lines of code refer to the airfoil and CP graph construction.

			DataPoint[] values = new DataPoint[nelems/2];
			DataPoint[] values1 = new DataPoint[nelems/2];

			DataPoint[] values2 = new DataPoint[nelems/2];
			DataPoint[] values3 = new DataPoint[nelems/2];

			for(int i =0; i<nelems/2; i++){

				DataPoint value, value1;

				if(isGround) {
					value = new DataPoint(rrc_scaledx1W[i], rrc_scaledy1W[i] - grndH);
					value1 = new DataPoint(rrc_scaledx2W[i], rrc_scaledy2W[i] - grndH);
				}else {
					value = new DataPoint(rrc_scaledx1W[i], rrc_scaledy1W[i]);
					value1 = new DataPoint(rrc_scaledx2W[i], rrc_scaledy2W[i]);
				}

				values[i] = value;
				values1[i] = value1;

				DataPoint cp = new DataPoint(rrc_scaledx1W[i], -cP1W[i]);
				DataPoint cp1 = new DataPoint(rrc_scaledx2W[i], -cP2W[i]);

				values2[i] = cp;
				values3[i] = cp1;
			}

			final GraphView graph = findViewById(R.id.graph);

			LineGraphSeries<DataPoint> series = new LineGraphSeries<>(values2);
			LineGraphSeries<DataPoint> series1 = new LineGraphSeries<>(values3);
			LineGraphSeries<DataPoint> series2 = new LineGraphSeries<>(values);
			LineGraphSeries<DataPoint> series3 = new LineGraphSeries<>(values1);
			graph.addSeries(series);
			graph.addSeries(series1);
			double max_x = 1.0*chLnW;
			double max_y, min_y;
			if(isGround) {
				max_y = Arrays.stream(rrc_scaledyW).max().getAsDouble() - grndH;
				min_y = Arrays.stream(rrc_scaledyW).min().getAsDouble() - grndH;
			}else{
				max_y = Arrays.stream(rrc_scaledyW).max().getAsDouble();
				min_y = Arrays.stream(rrc_scaledyW).min().getAsDouble();
			}
			graph.getViewport().setMaxX(max_x);
			series.setDrawDataPoints(true);
			series.setDataPointsRadius(6);
			series1.setDrawDataPoints(true);
			series1.setDataPointsRadius(6);
			graph.getViewport().setXAxisBoundsManual(true);
			graph.getViewport().setYAxisBoundsManual(true);
			graph.getViewport().setScalableY(true);
			graph.setTitle("-CP");
			graph.setTitleTextSize(40);
			graph.getSecondScale().addSeries(series2);
			graph.getSecondScale().addSeries(series3);
			graph.getSecondScale().setMinY(-5*abs(min_y));
			graph.getSecondScale().setMaxY(5*max_y);

			series2.setColor(Color.rgb(0,145,234));
			series3.setColor(Color.rgb(0,145,234));
			series.setColor(Color.rgb(230,81,0));
			series1.setColor(Color.rgb(126,87,194));
			graph.getGridLabelRenderer().setVerticalLabelsSecondScaleColor(Color.WHITE);

			// The code below refers to 2 tables, which show the CP behaviour as a function of the position
			// respectively on the lower and upper surface of the airfoil.

			TableLayout table = findViewById(R.id.table);
			TableLayout table1 = findViewById(R.id.table1);

			TableRow row_header = new TableRow(BuildCoefficient.this);
			TableRow row_header1 = new TableRow(BuildCoefficient.this);
			TextView x_c = new TextView(BuildCoefficient.this);
			x_c.setText(R.string.x_c);
			TextView x_c1 = new TextView(BuildCoefficient.this);
			x_c1.setText(R.string.x_c);
			TextView cP_text = new TextView(BuildCoefficient.this);
			cP_text.setText(R.string.cp);
			TextView cP_text1 = new TextView(BuildCoefficient.this);
			cP_text1.setText(R.string.cp);
			row_header.addView(x_c);
			row_header.addView(cP_text);
			row_header1.addView(x_c1);
			row_header1.addView(cP_text1);
			table.addView(row_header);
			table1.addView(row_header1);

			for(int i=0;i<nelems/2;i++){
				TableRow row = new TableRow(BuildCoefficient.this);
				TableRow row1 = new TableRow(BuildCoefficient.this);
				TextView lowerSurface = new TextView(BuildCoefficient.this);
				TextView upperSurface = new TextView(BuildCoefficient.this);
				lowerSurface.setText(String.format("%.3f", rrc_scaledxW[nelems/2-1-i]));
				lowerSurface.append("          ");
				upperSurface.setText(String.format("%.3f", rrc_scaledxW[nelems/2+i]));
				upperSurface.append("          ");
				TextView cpLow = new TextView(this);
				TextView cpUp = new TextView(this);
				cpLow.setText(String.format("%.3f", cPW[nelems/2-1-i]));
				cpUp.setText(String.format("%.3f", cPW[nelems/2+i]));
				row.addView(lowerSurface);
				row.addView(cpLow);
				row1.addView(upperSurface);
				row1.addView(cpUp);
				table.addView(row);
				table1.addView(row1);
			}

			if(isWT){
				DataPoint[] values4 = new DataPoint[nelems/2];
				DataPoint[] values5 = new DataPoint[nelems/2];

				DataPoint[] values6 = new DataPoint[nelems/2];
				DataPoint[] values7 = new DataPoint[nelems/2];

				for(int i =0; i<nelems/2; i++){

					DataPoint value4, value5;

					if(isGround) {
						value4 = new DataPoint(rrc_scaledx1T[i]*chLnT - refPointT[0], rrc_scaledy1T[i]*chLnT - wTY - grndH);
						value5 = new DataPoint(rrc_scaledx2T[i]*chLnT - refPointT[0], rrc_scaledy2T[i]*chLnT - wTY - grndH);
					}else{
						value4 = new DataPoint(rrc_scaledx1T[i]*chLnT - refPointT[0], rrc_scaledy1T[i]*chLnT - wTY);
						value5 = new DataPoint(rrc_scaledx2T[i]*chLnT - refPointT[0], rrc_scaledy2T[i]*chLnT - wTY);
					}

					values4[i] = value4;
					values5[i] = value5;

					DataPoint cp2 = new DataPoint(rrc_scaledx1T[i]*chLnT-refPointT[0], -cP1T[i]);
					DataPoint cp3 = new DataPoint(rrc_scaledx2T[i]*chLnT-refPointT[0], -cP2T[i]);

					values6[i] = cp2;
					values7[i] = cp3;
				}

				GraphView graphT = findViewById(R.id.graph_tail);
				graphT.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout = findViewById(R.id.tail);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.wing);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all13);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all14);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all15);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all16);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all17);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all18);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all19);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all20);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all21);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.all22);
				constraintLayout.setVisibility(View.VISIBLE);
				TextView textView = findViewById(R.id.titolo);
				textView.setVisibility(View.VISIBLE);
				textView = findViewById(R.id.titolo1);
				textView.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.s2);
				constraintLayout.setVisibility(View.VISIBLE);
				constraintLayout = findViewById(R.id.d2);
				constraintLayout.setVisibility(View.VISIBLE);

				LineGraphSeries<DataPoint> series4 = new LineGraphSeries<>(values6);
				LineGraphSeries<DataPoint> series5 = new LineGraphSeries<>(values7);
				LineGraphSeries<DataPoint> series6 = new LineGraphSeries<>(values4);
				LineGraphSeries<DataPoint> series7 = new LineGraphSeries<>(values5);
				graphT.addSeries(series4);
				graphT.addSeries(series5);
				max_x = 1.0*chLnT;
				if (isGround) {
					max_y = Arrays.stream(rrc_scaledyT).max().getAsDouble()-grndH;
					min_y = Arrays.stream(rrc_scaledyT).min().getAsDouble()-grndH;
				}else {
					max_y = Arrays.stream(rrc_scaledyT).max().getAsDouble();
					min_y = Arrays.stream(rrc_scaledyT).min().getAsDouble();
				}
				graphT.getViewport().setMaxX(max_x);
				series4.setDrawDataPoints(true);
				series4.setDataPointsRadius(6);
				series5.setDrawDataPoints(true);
				series5.setDataPointsRadius(6);
				graphT.getViewport().setXAxisBoundsManual(true);
				graphT.getViewport().setYAxisBoundsManual(true);
				graphT.getViewport().setScalableY(true);
				graphT.setTitle("-CP");
				graphT.setTitleTextSize(40);
				graphT.getSecondScale().addSeries(series6);
				graphT.getSecondScale().addSeries(series7);
				graphT.getSecondScale().setMinY(-1.5*abs(min_y));
				graphT.getSecondScale().setMaxY(1.5*max_y);

				series6.setColor(Color.rgb(0,145,234));
				series7.setColor(Color.rgb(0,145,234));
				series4.setColor(Color.rgb(230,81,0));
				series5.setColor(Color.rgb(126,87,194));
				graphT.getGridLabelRenderer().setVerticalLabelsSecondScaleColor(Color.WHITE);

				// The code below refers to 2 tables, which show the CP behaviour as a function of the position
				// respectively on the lower and upper surface of the airfoil.

				table = findViewById(R.id.table_tail);
				table1 = findViewById(R.id.table1_tail);

				row_header  = new TableRow(this);
				row_header1 = new TableRow(this);
				TextView x_cT         = new TextView(this);
				x_cT.setText(R.string.x_c);
				TextView x_c1T        = new TextView(this);
				x_c1T.setText(R.string.x_c);
				TextView cP_textT     = new TextView(this);
				cP_textT.setText(R.string.cp);
				TextView cP_text1T    = new TextView(this);
				cP_text1T.setText(R.string.cp);
				row_header.addView(x_cT);
				row_header.addView(cP_textT);
				row_header1.addView(x_c1T);
				row_header1.addView(cP_text1T);
				table.addView(row_header);
				table1.addView(row_header1);

				for(int i=0;i<nelems/2;i++){
					TableRow row = new TableRow(this);
					TableRow row1 = new TableRow(this);
					TextView lowerSurface = new TextView(this);
					TextView upperSurface = new TextView(this);
					lowerSurface.setText(String.format("%.3f", rrc_scaledxT[nelems/2-1-i]-refPointT[0]));
					lowerSurface.append("          ");
					upperSurface.setText(String.format("%.3f", rrc_scaledxT[nelems/2+i]-refPointT[0]));
					upperSurface.append("          ");
					TextView cpLow = new TextView(this);
					TextView cpUp = new TextView(this);
					cpLow.setText(String.format("%.3f", cPT[nelems/2-1-i]));
					cpUp.setText(String.format("%.3f", cPT[nelems/2+i]));
					row.addView(lowerSurface);
					row.addView(cpLow);
					row1.addView(upperSurface);
					row1.addView(cpUp);
					table.addView(row);
					table1.addView(row1);
				}
			}

			//The following lines of code contain the implementation of the Thwaites method for the
			//viscous flow corrections.

			int n_stg = 0;
			double csi_stg = 0.0;
			int i_stg = 0;

			for(int i = 0; i<nelems-1; i++){
				if(vTiW[i]*vTiW[i+1]<=0) {
					n_stg++;
					i_stg = i;
					csi_stg = -vTiW[i]/(vTiW[i+1] - vTiW[i]);
				}
			}

			double[] s = new double[nelems];

			s[i_stg+1] = csi_stg*(lenW[i_stg] + lenW[i_stg+1]) * 0.5;
			s[i_stg] = (1.0-csi_stg) * (lenW[i_stg] + lenW[i_stg+1]) * 0.5;

			for(int i = i_stg+2; i < nelems; i++){
				s[i] = s[i-1] + 0.5*(lenW[i] + lenW[i-1]);
			}
			for(int i = i_stg-1; i>-1; i--){
				s[i] = s[i+1] + 0.5*(lenW[i] + lenW[i+1]);
			}

			double th = Math.acos(nversW[0][i_stg]*nversW[0][i_stg+1]+nversW[1][i_stg]*nversW[1][i_stg+1]);
			double r0 = lenW[1]*Math.sin(th) + (lenW[0] + lenW[1]*Math.cos(th))/Math.tan(th);
			double k = freeStream.v/r0;
			double theta0 = Math.sqrt(0.075*freeStream.kin_visc/k);  //never used

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
			double[] UeThetaH1 			= 		new double[nelems];
			double[] H1					= 		new double[nelems];

			double[] Res				=		new double[Ue.length];
			double[] P					=		new double[cPW.length];
			for(int i = 0; i<Ue.length;i++){
				Res[i] = Ue[i]*s[i] / freeStream.kin_visc;
				P[i] = cPW[i]*(0.5*freeStream.rho*freeStream.v*freeStream.v);
			}

			double[] dPdx			= 		new double[nelems];
			double[] dUedx			= 		new double[nelems];

			double Ptot = freeStream.P + 0.5*freeStream.rho * freeStream.v * freeStream.v;

			dPdx[i_stg+1]	 = 	(P[i_stg+1] - Ptot) / (0.5*csi_stg*(lenW[i_stg]+ lenW[i_stg+1]));
			dPdx[i_stg]		 = 	(P[i_stg] - Ptot) / (0.5*(1.0-csi_stg)*(lenW[i_stg]+ lenW[i_stg+1]));
			dUedx[i_stg+1]	 = 	(Ue[i_stg+1] - 0.0) / (0.5*csi_stg*(lenW[i_stg]+ lenW[i_stg+1]));
			dUedx[i_stg]	 = 	(Ue[i_stg] -   0.0) / (0.5*(1.0-csi_stg)*(lenW[i_stg]+ lenW[i_stg+1]));

			theta2Ue6[i_stg+1] = 0.45*freeStream.kin_visc*Math.pow(Ue[i_stg+1],5.0)*
					(0.25*csi_stg*(lenW[i_stg]+lenW[i_stg+1]));
			theta2Ue6[i_stg]   = 0.45*freeStream.kin_visc*Math.pow(Ue[i_stg],5.0)*
					(0.25*(1.0-csi_stg)*(lenW[i_stg]+lenW[i_stg+1]));
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
					dPdx[i]  = (P[i+1] - P[i-1]) / (0.5 * (lenW[i+1] + lenW[i-1]) + lenW[i]);
					dUedx[i] = ( Ue[i+1] - Ue[i-1] ) / ( 0.5 * ( lenW[i+1] + lenW[i-1] ) + lenW[i]);
				}
				else{
					dPdx[i]  = (P[i] - P[i-1]) / (0.5 * (lenW[i] + lenW[i-1]));
					dUedx[i] = ( Ue[i] - Ue[i-1] ) / ( 0.5 * ( lenW[i] + lenW[i-1] ));
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

			for(int i = i_stg-1;i>-1;i--){
				if(i>0){
					dPdx[i]  = (P[i-1]-P[i+1])/(0.5*(lenW[i+1]+lenW[i-1])+lenW[i]);
					dUedx[i] = (Ue[i-1]-Ue[i+1])/(0.5*(lenW[i+1]+lenW[i-1])+lenW[i]);
				}
				else{
					dPdx[i]  = (P[i]-P[i+1])/(0.5*(lenW[i]+lenW[i+1]));
					dUedx[i] = (Ue[i]-Ue[i+1])/(0.5*(lenW[i]+lenW[i+1]));
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

			double[] tauW = new double[cf.length];

			for(int i = 0; i<cf.length;i++){
				tauW[i] = 0.5*freeStream.rho*Math.pow(freeStream.v,2)*cf[i];
			}

			double[][] temp = new double[2][tauW.length];
			double[][] dF_visc = new double[2][tauW.length];

			for(int i = 0; i<tauW.length;i++){
				temp[0][i] = tauW[i]*lenW[i]*vTiW[i]/Ue[i];
				temp[1][i] = tauW[i]*lenW[i]*vTiW[i]/Ue[i];
			}

			for(int i = 0; i<tauW.length;i++){
				dF_visc[0][i] = temp[0][i]*tversW[0][i];
				dF_visc[1][i] = temp[1][i]*tversW[1][i];
			}
			double[] F_visc = new double[2];

			F_visc[0] = DoubleStream.of(dF_visc[0]).sum()/airfoilW.chord;
			F_visc[1] = DoubleStream.of(dF_visc[1]).sum()/airfoilW.chord;

			double L_visc = -F_visc[0]*Math.sin(freeStream.alpha)+F_visc[1]*Math.cos(freeStream.alpha);
			double D_visc = F_visc[0]*Math.cos(freeStream.alpha) + F_visc[1]*Math.sin(freeStream.alpha);

			double[][] temp0 = new double[2][tauW.length];
			double[][] dF_Pres = new double[2][tauW.length];

			for(int i = 0; i<tauW.length;i++){
				temp0[0][i] = -P[i]*lenW[i];
				temp0[1][i] = -P[i]*lenW[i];
			}
			for(int i = 0; i<tauW.length;i++){
				dF_Pres[0][i] = temp0[0][i]*nversW[0][i];
				dF_Pres[1][i] = temp0[1][i]*nversW[1][i];
			}
			double[] F_pres = new double[2];

			F_pres[0] = DoubleStream.of(dF_Pres[0]).sum()/airfoilW.chord;
			F_pres[1] = DoubleStream.of(dF_Pres[1]).sum()/airfoilW.chord;

			double L_pres = -F_pres[0]*Math.sin(freeStream.alpha)+F_pres[1]*Math.cos(freeStream.alpha);
			double D_pres = F_pres[0]*Math.cos(freeStream.alpha) + F_pres[1]*Math.sin(freeStream.alpha);

			double L = L_pres + L_visc;
			double D = D_pres + D_visc;

			double cL = L/(0.5*freeStream.rho*freeStream.v*freeStream.v);
			double cD = D/(0.5*freeStream.rho*freeStream.v*freeStream.v);

			if(isWT){
				n_stg = 0;
				csi_stg = 0.0;
				i_stg = 0;

				for(int i = 0; i<nelems-1; i++){
					if(vTiT[i]*vTiT[i+1]<=0) {
						n_stg++;
						i_stg = i;
						csi_stg = -vTiT[i]/(vTiT[i+1] - vTiT[i]);
					}
				}

				s = new double[nelems];

				s[i_stg+1] = csi_stg*(lenT[i_stg] + lenT[i_stg+1]) * 0.5;
				s[i_stg] = (1.0-csi_stg) * (lenT[i_stg] + lenT[i_stg+1]) * 0.5;

				for(int i = i_stg+2; i < nelems; i++){
					s[i] = s[i-1] + 0.5*(lenT[i] + lenT[i-1]);
				}
				for(int i = i_stg-1; i>-1; i--){
					s[i] = s[i+1] + 0.5*(lenT[i] + lenT[i+1]);
				}

				th = Math.acos(nversT[0][i_stg]*nversT[0][i_stg+1]+nversT[1][i_stg]*nversT[1][i_stg+1]);
				r0 = lenT[1]*Math.sin(th) + (lenT[0] + lenT[1]*Math.cos(th))/Math.tan(th);
				k = freeStream.v/r0;
				theta0 = Math.sqrt(0.075*freeStream.kin_visc/k);  //never used

				Ue = new double[vTiT.length];
				for(int i = 0; i<vTiT.length; i++){
					if(vTiT[i] > 0)
						Ue[i] = vTiT[i];
					else
						Ue[i] = -vTiT[i];
				}

				theta2Ue6 			= 		new double[nelems];
				theta 				= 		new double[nelems];
				ReTheta			    =		new double[nelems];
				lambda 			    =		new double[nelems];
				ell 				= 		new double[nelems];
				H					=		new double[nelems];
				delta				= 		new double[nelems];
				cf					= 		new double[nelems];
				UeThetaH1 			= 		new double[nelems];
				H1					= 		new double[nelems];

				Res				    =		new double[Ue.length];
				P					=		new double[cPT.length];

				for(int i = 0; i<Ue.length;i++){
					Res[i] = Ue[i]*s[i] / freeStream.kin_visc;
					P[i] = cPT[i]*(0.5*freeStream.rho*freeStream.v*freeStream.v);
				}

				dPdx			= 		new double[nelems];
				dUedx			= 		new double[nelems];

				Ptot = freeStream.P + 0.5*freeStream.rho * freeStream.v * freeStream.v;

				dPdx[i_stg+1]	 = 	(P[i_stg+1] - Ptot) / (0.5*csi_stg*(lenT[i_stg]+ lenT[i_stg+1]));
				dPdx[i_stg]		 = 	(P[i_stg] - Ptot)   / (0.5*(1.0-csi_stg)*(lenT[i_stg]+ lenT[i_stg+1]));
				dUedx[i_stg+1]	 = 	(Ue[i_stg+1] - 0.0) / (0.5*csi_stg*(lenT[i_stg]+ lenT[i_stg+1]));
				dUedx[i_stg]	 = 	(Ue[i_stg] -   0.0) / (0.5*(1.0-csi_stg)*(lenT[i_stg]+ lenT[i_stg+1]));

				theta2Ue6[i_stg+1] = 0.45*freeStream.kin_visc*Math.pow(Ue[i_stg+1],5.0)*
						(0.25*csi_stg*(lenT[i_stg]+lenT[i_stg+1]));
				theta2Ue6[i_stg]   = 0.45*freeStream.kin_visc*Math.pow(Ue[i_stg],5.0)*
						(0.25*(1.0-csi_stg)*(lenT[i_stg]+lenT[i_stg+1]));
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


				regime = "laminar";
				transition = "no";

				for(int i = i_stg+2; i<nelems; i++){

					if(i<nelems-1){
						dPdx[i]  = (P[i+1] - P[i-1]) / (0.5 * (lenT[i+1] + lenT[i-1]) + lenT[i]);
						dUedx[i] = ( Ue[i+1] - Ue[i-1] ) / ( 0.5 * ( lenT[i+1] + lenT[i-1] ) + lenT[i]);
					}
					else{
						dPdx[i]  = (P[i] - P[i-1]) / (0.5 * (lenT[i] + lenT[i-1]));
						dUedx[i] = ( Ue[i] - Ue[i-1] ) / ( 0.5 * ( lenT[i] + lenT[i-1] ));
					}

					if(regime.equals("laminar")) {
						theta2Ue6[i] = theta2Ue6[i - 1] + 0.45 * freeStream.kin_visc * (Math.pow(Ue[i], 5.0) +
								Math.pow(Ue[i - 1], 5.0)) * (0.25 * (lenT[i] + lenT[i - 1]));

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
						}

					} else {

						double dx1 = (0.5 * (lenT[i] + lenT[i - 1]));
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

				for(int i = i_stg-1;i>-1;i--){
					if(i>0){
						dPdx[i]  = (P[i-1]-P[i+1])/(0.5*(lenT[i+1]+lenT[i-1])+lenT[i]);
						dUedx[i] = (Ue[i-1]-Ue[i+1])/(0.5*(lenT[i+1]+lenT[i-1])+lenT[i]);
					}
					else{
						dPdx[i]  = (P[i]-P[i+1])/(0.5*(lenT[i]+lenT[i+1]));
						dUedx[i] = (Ue[i]-Ue[i+1])/(0.5*(lenT[i]+lenT[i+1]));
					}

					if(regime.equals("laminar")){
						theta2Ue6[i] = theta2Ue6[i + 1] + 0.45 * freeStream.kin_visc * (Math.pow(Ue[i], 5.0) +
								Math.pow(Ue[i + 1], 5.0)) * (0.25 * (lenT[i] + lenT[i + 1]));
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
						}

					}  else {

						double dx1 = (0.5 * (lenT[i] + lenT[i + 1]));
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

				double[] tauT = new double[cf.length];

				for(int i = 0; i<cf.length;i++){
					tauT[i] = 0.5*freeStream.rho*Math.pow(freeStream.v,2)*cf[i];
				}

				temp = new double[2][tauT.length];
				dF_visc = new double[2][tauT.length];

				for(int i = 0; i<tauT.length;i++){
					temp[0][i] = tauT[i]*lenT[i]*vTiT[i]/Ue[i];
					temp[1][i] = tauT[i]*lenT[i]*vTiT[i]/Ue[i];
				}

				for(int i = 0; i<tauT.length;i++){
					dF_visc[0][i] = temp[0][i]*tversT[0][i];
					dF_visc[1][i] = temp[1][i]*tversT[1][i];
				}
				F_visc = new double[2];

				F_visc[0] = DoubleStream.of(dF_visc[0]).sum()/airfoilT.chord;
				F_visc[1] = DoubleStream.of(dF_visc[1]).sum()/airfoilT.chord;

				double L_viscT = -F_visc[0]*Math.sin(freeStream.alpha)+F_visc[1]*Math.cos(freeStream.alpha);
				double D_viscT = F_visc[0]*Math.cos(freeStream.alpha) + F_visc[1]*Math.sin(freeStream.alpha);

				temp0 = new double[2][tauT.length];
				dF_Pres = new double[2][tauT.length];

				for(int i = 0; i<tauT.length;i++){
					temp0[0][i] = -P[i]*lenT[i];
					temp0[1][i] = -P[i]*lenT[i];
				}
				for(int i = 0; i<tauT.length;i++){
					dF_Pres[0][i] = temp0[0][i]*nversT[0][i];
					dF_Pres[1][i] = temp0[1][i]*nversT[1][i];
				}
				F_pres = new double[2];

				F_pres[0] = DoubleStream.of(dF_Pres[0]).sum()/airfoilT.chord;
				F_pres[1] = DoubleStream.of(dF_Pres[1]).sum()/airfoilT.chord;

				double L_presT = -F_pres[0]*Math.sin(freeStream.alpha)+F_pres[1]*Math.cos(freeStream.alpha);
				double D_presT = F_pres[0]*Math.cos(freeStream.alpha) + F_pres[1]*Math.sin(freeStream.alpha);

				double LT = L_presT + L_viscT;
				double DT = D_presT + D_viscT;

				double cLT = LT/(0.5*freeStream.rho*freeStream.v*freeStream.v);
				double cDT = DT/(0.5*freeStream.rho*freeStream.v*freeStream.v);

				TextView cltextthwaites = findViewById(R.id.textCltwhaitesvalueTail);
				cltextthwaites.setText(String.format("%.5f", cLT)+ " [-]");

				TextView cdtextthwaites = findViewById(R.id.textCdtwhaitesvalueTail);
				cdtextthwaites.setText(String.format("%.5f", cDT)+ " [-]");
			}




			//Print the following results: CL, CLKJ, Cm CL_thwaites, CD_thwaites

			TextView cltext = findViewById(R.id.textClvalue);
			cltext.setText(String.format("%.5f", cl) + " [-]");

			TextView cltextzukowski = findViewById(R.id.textClKJvalue);
			cltextzukowski.setText(String.format("%.5f", cLKJW) + " [-]");

			TextView cm = findViewById(R.id.textCmvalue);
			cm.setText(String.format("%.5f", CmW) + " [-]");

			TextView cltextthwaites = findViewById(R.id.textCltwhaitesvalue);
			cltextthwaites.setText(String.format("%.5f", cL)+ " [-]");

			TextView cdtextthwaites = findViewById(R.id.textCdtwhaitesvalue);
			cdtextthwaites.setText(String.format("%.5f", cD)+ " [-]");

			if(isWT){
				cltext = findViewById(R.id.textClvalueTail);
				cltext.setText(String.format("%.5f", clT) + " [-]");

				cltextzukowski = findViewById(R.id.textClKJvalueTail);
				cltextzukowski.setText(String.format("%.5f", clKJT) + " [-]");

				cm = findViewById(R.id.textCmvalueTail);
				cm.setText(String.format("%.5f", CmT) + " [-]");
			}
		}
	}
