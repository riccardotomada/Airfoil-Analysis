package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//This class contains all the functions shared between the activities. 						      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import java.lang.reflect.Array;

public class GlobalFunctions {

	// All the following methods, except for the ones which are explicitly commented, are the copy
	// of the MATLAB functions we used in the fluid dynamics laboratory.

	double[][] getCoordinates(Airfoil airfoil, boolean checkGround) {

		String naca = "NACA";
		for (int i = 0; i < 4; i++) {
			if (airfoil.airfoil_str.charAt(i) != naca.charAt(i)) {
				return new double[][]{{0.0}, {0.0}};
			}
		}

		double MC, POS, TH;
		char m, p;
		char[] th, mm;
		int coeff;

		if (airfoil.airfoil_str.length() == 8) {
			m = airfoil.airfoil_str.charAt(4);
			p = airfoil.airfoil_str.charAt(5);
			th = new char[]{airfoil.airfoil_str.charAt(6), airfoil.airfoil_str.charAt(7)};
			MC = Double.parseDouble(String.valueOf(m));
			POS = Double.parseDouble(String.valueOf(p));
			TH = Double.parseDouble(String.valueOf(th));

			if (checkGround) {
				return setAirfoil4(MC, POS, TH, airfoil.chord, airfoil.nChordPanels,
						airfoil.refPoint[0], -airfoil.refPoint[1], airfoil.xcRefPoint, airfoil.theta, true); }
			else{
				return setAirfoil4(MC, POS, TH, airfoil.chord, airfoil.nChordPanels,
						airfoil.refPoint[0], airfoil.refPoint[1], airfoil.xcRefPoint, airfoil.theta, false);
			}
		} else {
			mm = new char[]{airfoil.airfoil_str.charAt(4), airfoil.airfoil_str.charAt(5),
					airfoil.airfoil_str.charAt(6)};
			th = new char[]{airfoil.airfoil_str.charAt(7), airfoil.airfoil_str.charAt(8)};
			coeff = Integer.parseInt(String.valueOf(mm));
			TH = Double.parseDouble(String.valueOf(th));
			if (checkGround) {
				return setAirfoil5(coeff, TH, airfoil.chord, airfoil.nChordPanels,
						airfoil.refPoint[0], airfoil.refPoint[1], airfoil.xcRefPoint, airfoil.theta, true); }
			else{
				return setAirfoil5(coeff, TH, airfoil.chord, airfoil.nChordPanels,
						airfoil.refPoint[0], airfoil.refPoint[1], airfoil.xcRefPoint, airfoil.theta, false);
			}
		}
	}

	Elems[] elemsMethod(Airfoil airfoil, double[][]rr, int[][]ee, int nelems){
		Elems[] elems = new Elems[nelems];

		for(int i = 0; i<nelems;i++){
			elems[i]           = new Elems();
			elems[i].airfoilId = airfoil.id;
			elems[i].id        = i;
			elems[i].ver1      = new double[] {rr[0][ee[0][i]-1], rr[1][ee[0][i]-1]};
			elems[i].ver2      = new double[] {rr[0][ee[1][i]-1], rr[1][ee[1][i]-1]};
			elems[i].cen       = midpoint(elems[i].ver1, elems[i].ver2);
			elems[i].len       = len_norm(elems[i].ver2, elems[i].ver1);
			elems[i].tver      = tver(elems[i].ver2, elems[i].ver1, elems[i].len);
			elems[i].nver      = new double[]{-elems[i].tver[1], elems[i].tver[0]};
		}
		return elems;
	}

	double[][] naca4digit(double MC, double POS, double TH, double c, int n) {

		double mc = MC / 100;
		double pos = POS / 10;
		double th = TH / 100;

		double[] xVector = linspace(0, c, n + 1);
		for (int i = 0; i < Array.getLength(xVector); i++)
			xVector[i] = c * (1 - 0.5 * (1 + Math.cos(Math.PI * xVector[i] / c)));

		final double a0 = 0.2969;
		final double a1 = -0.126;
		final double a2 = -0.3516;
		final double a3 = 0.2843;
		final double a4 = -0.1015;

		double[] yth = new double[Array.getLength(xVector)];
		for (int i = 0; i < Array.getLength(yth); i++) {
			yth[i] = th / 0.2 * c * (a0 * Math.pow(xVector[i] / c, 0.5) + a1 * xVector[i] / c +
					a2 * Math.pow(xVector[i] / c, 2) + a3 * Math.pow(xVector[i] / c, 3) +
					a4 * Math.pow(xVector[i] / c, 4));
		}

		double[] ycamber = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			if (xVector[i] < pos*c)
				ycamber[i] = c*(mc / Math.pow(pos, 2) * (xVector[i]/c) * (2 * pos - (xVector[i]/c)));
			else
				ycamber[i] = c*(mc / Math.pow((1 - pos), 2) * (1 + (2*pos - (xVector[i]/c)) * xVector[i]/c -2*pos));
		}

		double[] dycamber = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			if (xVector[i] < pos*c)
				dycamber[i] = 2 * mc / Math.pow(pos, 2) * (pos - xVector[i]/c);
			else
				dycamber[i] = 2 * mc / Math.pow((1 - pos), 2) * (pos - xVector[i]/c);
		}

		double[] theta = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			theta[i] = Math.atan2(dycamber[i], 1);
		}

		double[] xU = new double[Array.getLength((xVector))];
		double[] yU = new double[Array.getLength((xVector))];
		double[] xL = new double[Array.getLength((xVector))];
		double[] yL = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			xU[i] = xVector[i] - yth[i] * Math.sin(theta[i]);
			yU[i] = ycamber[i] + yth[i] * Math.cos(theta[i]);
			xL[i] = xVector[i] + yth[i] * Math.sin(theta[i]);
			yL[i] = ycamber[i] - yth[i] * Math.cos(theta[i]);
		}

		double[] x = new double[2 * n + 1];
		double[] y = new double[2 * n + 1];

		for (int i = 0; i < n; i++) {
			x[i] = xL[n - i];
			y[i] = yL[n - i];
		}
		for (int i = n; i < (2 * n + 1); i++) {
			x[i] = xU[i - n];
			y[i] = yU[i - n];
		}

		double[][] xy = new double[2][Array.getLength(x)];

		for (int i = 0; i < Array.getLength(x); i++) {
			xy[0][i] = x[i];
			xy[1][i] = y[i];
		}
		return xy;

	}

	double[][] naca5digit(int coeff, double TH, double c, int n) {

		double t = TH / 100;

		double[] xVector = linspace(0, c, n + 1);
		for (int i = 0; i < Array.getLength(xVector); i++)
			xVector[i] = c / 2 * (1 - Math.cos(Math.PI * xVector[i] / c));

		double q, k;

		switch (coeff) {
			case 210:
				q = 0.0580;
				k = 361.4;
				break;
			case 220:
				q = 0.1260;
				k = 51.64;
				break;
			case 230:
				q = 0.2025;
				k = 15.957;
				break;
			case 240:
				q = 0.2900;
				k = 6.643;
				break;
			case 250:
				q = 0.3910;
				k = 3.230;
				break;
			default:
				q = 0;
				k = 0;
		}
		final double a0 = 0.2969;
		final double a1 = -0.126;
		final double a2 = -0.3516;
		final double a3 = 0.2843;
		final double a4 = -0.1015;

		double[] yth = new double[Array.getLength(xVector)];
		for (int i = 0; i < Array.getLength(yth); i++) {
			yth[i] = t / 0.2 * c * (a0 * Math.pow(xVector[i] / c, 0.5) + a1 * xVector[i] / c +
					a2 * Math.pow(xVector[i] / c, 2) + a3 * Math.pow(xVector[i] / c, 3) +
					a4 * Math.pow(xVector[i] / c, 4));
		}

		double[] ycamber = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			if (xVector[i] < q * c)
				ycamber[i] = c * (k / 6 * (Math.pow(xVector[i] / c, 3) - 3 * q * Math.pow(xVector[i] / c, 2)
						+ q * q * (3 - q) * (xVector[i] / c)));
			else
				ycamber[i] = c * k / 6 * q * q * q * (1 - (xVector[i] / c));
		}

		double[] dycamber = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			if (xVector[i] < q * c)
				dycamber[i] = k / 6 * (3 * Math.pow(xVector[i] / c, 2) - 6 * q * xVector[i] / c + q * q * (3 - q));
			else
				dycamber[i] = -k / 6 * q * q * q;
		}
		double[] theta = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			theta[i] = Math.atan2(dycamber[i], 1);
		}

		double[] xU = new double[Array.getLength((xVector))];
		double[] yU = new double[Array.getLength((xVector))];
		double[] xL = new double[Array.getLength((xVector))];
		double[] yL = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			xU[i] = xVector[i] - yth[i] * Math.sin(theta[i]);
			yU[i] = ycamber[i] + yth[i] * Math.cos(theta[i]);
			xL[i] = xVector[i] + yth[i] * Math.sin(theta[i]);
			yL[i] = ycamber[i] - yth[i] * Math.cos(theta[i]);
		}

		double[] x = new double[2 * n + 1];
		double[] y = new double[2 * n + 1];

		for (int i = 0; i < n; i++) {
			x[i] = xL[n - i];
			y[i] = yL[n - i];
		}
		for (int i = n; i < (2 * n + 1); i++) {
			x[i] = xU[i - n];
			y[i] = yU[i - n];
		}

		double[][] xy = new double[2][Array.getLength(x)];

		for (int i = 0; i < Array.getLength(x); i++) {
			xy[0][i] = x[i];
			xy[1][i] = y[i];
		}
		return xy;
	}

	double[][] setAirfoil4(double MC, double POS, double TH, double c, int n, double x0, double y0,
						   double ca, double theta, boolean checkGround) {
		double[][] ab;

		ab = naca4digit(MC, POS, TH, c, n);

		if(checkGround) {
			for (int i = 0; i < 2 * n + 1; i++) {
				ab[1][i] = -ab[1][i];
			}
			ab[0] = reverseArray(ab[0],ab[0].length);
			ab[1] = reverseArray(ab[1],ab[1].length);
		}

		double c1 = ca * c;

		double[] x = new double[2 * n + 1];
		double[] y = new double[2 * n + 1];

		for (int i = 0; i < Array.getLength(x); i++) {

			if(checkGround){
				x[i] = x0 + c1 + Math.cos(-theta) * (ab[0][i] - c1) - Math.sin(-theta) * ab[1][i];
				y[i] = y0 + Math.sin(-theta) * (ab[0][i] - c1) + Math.cos(-theta) * ab[1][i];
			}
			else {
				x[i] = x0 + c1 + Math.cos(-theta) * (ab[0][i] - c1) + Math.sin(-theta) * ab[1][i];
				y[i] = y0 - Math.sin(-theta) * (ab[0][i] - c1) + Math.cos(-theta) * ab[1][i];
			}
		}

		double[][] xy = new double[2][2 * n + 1];
		for (int i = 0; i < Array.getLength(x); i++) {
			xy[0][i] = x[i];
			xy[1][i] = y[i];
		}

		return xy;
	}

	double[][] setAirfoil5(int coeff, double TH, double c, int n, double x0, double y0,
						   double ca, double theta, boolean checkGround) {
		double[][] ab;

		ab = naca5digit(coeff, TH, c, n);

		if(checkGround)
			for(int i = 0; i<2*n+1; i++){
				ab[1][i] = -ab[1][i];
				ab[0][i] = ab[0][2*n-i];
				ab[1][i] = ab[1][2*n-i];
			}

		double c1 = ca * c;

		double[] x = new double[2 * n + 1];
		double[] y = new double[2 * n + 1];

		for (int i = 0; i < Array.getLength(x); i++) {
			if(checkGround){
				x[i] = x0 + c1 + Math.cos(-theta) * (ab[0][i] - c1) - Math.sin(-theta) * ab[1][i];
				y[i] = y0 + Math.sin(-theta) * (ab[0][i] - c1) + Math.cos(-theta) * ab[1][i];
			}
			else {
				x[i] = x0 + c1 + Math.cos(-theta) * (ab[0][i] - c1) + Math.sin(-theta) * ab[1][i];
				y[i] = y0 - Math.sin(-theta) * (ab[0][i] - c1) + Math.cos(-theta) * ab[1][i];
			}
		}

		double[][] xy = new double[2][2 * n + 1];
		for (int i = 0; i < Array.getLength(x); i++) {
			xy[0][i] = x[i];
			xy[1][i] = y[i];
		}

		return xy;
	}

	int[][] connectivityMatrix(double[] x) {
		int[][] conMat = new int[2][Array.getLength(x) - 1];
		for (int i = 1; i < Array.getLength(x); i++) {
			conMat[0][i - 1] = i;
			conMat[1][i - 1] = i + 1;
		}
		return conMat;
	}

	double[] computeVelocityVortex(Elems elem_j, Elems elem_i, boolean diffAirfoil) {
		double[] r1 = new double[2];
		double[] r2 = new double[2];
		double[] vstar, v;
		double sij, bij, vcross, sinbij, cosbij;

		r1[0] = elem_i.cen[0] - elem_j.ver1[0];
		r1[1] = elem_i.cen[1] - elem_j.ver1[1];
		r2[0] = elem_i.cen[0] - elem_j.ver2[0];
		r2[1] = elem_i.cen[1] - elem_j.ver2[1];

		if (elem_i.id == elem_j.id && !diffAirfoil) {
			sij = 1.0;
			bij = Math.PI;
		} else {
			sij = Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2)) /
					Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2));
			vcross = cross(r1, r2);
			sinbij = vcross / Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2)) /
					Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2));
			cosbij = scalar2(r1, r2) / Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2)) /
					Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2));
			bij = Math.atan2(sinbij, cosbij);
		}
		vstar = new double[]{-bij / (2 * Math.PI), -Math.log(sij) / (2 * Math.PI)};

		v = new double[]{elem_j.tver[0] * vstar[0] + elem_j.nver[0] * vstar[1],
				elem_j.tver[1] * vstar[0] + elem_j.nver[1] * vstar[1]};
		return v;
	}

	double[] computeVelocitySource(Elems elem_j, Elems elem_i, boolean diffAirfoil) {
		double[] r1 = new double[2];
		double[] r2 = new double[2];
		double[] vstar, v;
		double sij, bij, vcross, sinbij, cosbij;

		r1[0] = elem_i.cen[0] - elem_j.ver1[0];
		r1[1] = elem_i.cen[1] - elem_j.ver1[1];
		r2[0] = elem_i.cen[0] - elem_j.ver2[0];
		r2[1] = elem_i.cen[1] - elem_j.ver2[1];

		if (elem_i.id == elem_j.id && !diffAirfoil) {
			sij = 1.0;
			bij = Math.PI;
		} else {
			sij = Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2)) /
					Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2));
			vcross = cross(r1, r2);
			sinbij = vcross / Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2)) /
					Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2));
			cosbij = scalar2(r1, r2) / Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2)) /
					Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2));
			bij = Math.atan2(sinbij, cosbij);
		}


		vstar = new double[]{-Math.log(sij) / (2 * Math.PI), bij / (2 * Math.PI)};
		v = new double[]{elem_j.tver[0] * vstar[0] + elem_j.nver[0] * vstar[1],
				elem_j.tver[1] * vstar[0] + elem_j.nver[1] * vstar[1]};
		return v;
	}

	// The midpoint method returns an array whose elements are the average of the correspondent
	// elements of the given arrays.
	double[] midpoint(double[] a, double[] b) {
		double[] c = new double[Array.getLength(a)];
		for (int i = 0; i < Array.getLength(a); ++i) {
			c[i] = a[i] + b[i];
			c[i] = c[i] / 2;
		}
		return c;
	}

	// The len_norm method compute the norm of a 2 dimension array which is obtained by the difference
	// of the given arrays.
	double len_norm(double[] a, double[] b) {
		double[] c = new double[2];
		for (int i = 0; i < 2; i++) {
			c[i] = a[i] - b[i];
		}
		return Math.sqrt(c[0] * c[0] + c[1] * c[1]);
	}

	double[] tver(double[] a, double[] b, double norm) {
		double[] c = new double[2];
		for (int i = 0; i < 2; i++) {
			c[i] = (a[i] - b[i]) / norm;
		}
		return c;
	}

	// The following method is the correspondent of a part of the Matlab function build_lin_sys
	// and returns the coefficient matrix A.

	double[][] coefficientMatrix(Elems[] elems, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] A = new double[nelems + n_te][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = computeVelocitySource(elems[j], elems[i], false);
				vv_ij = computeVelocityVortex(elems[j], elems[i], false);

				A[i][j]      = scalar2(elems[i].nver, vs_ij);
				A[i][nelems] = A[i][nelems] + scalar2(elems[i].nver, vv_ij);
				if (i == 0 || i == nelems - 1) {
					A[nelems][j]      = A[nelems][j] + scalar2(elems[i].tver, vs_ij);
					A[nelems][nelems] = A[nelems][nelems] + scalar2(elems[i].tver, vv_ij);
				}

			}
		}
		return A;
	}

	double[][] coefficientMatrixWT(Elems[] elemsW, Elems[] elemsT, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 2;

		double[][] A = new double[2*nelems + n_te][2*nelems + n_te];
		double[] vs_ij11, vv_ij11, vs_ij12, vv_ij12, vs_ij21, vv_ij21, vs_ij22, vv_ij22;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij11 = computeVelocitySource(elemsW[j], elemsW[i], false);
				vv_ij11 = computeVelocityVortex(elemsW[j], elemsW[i], false);

				vs_ij12 = computeVelocitySource(elemsT[j], elemsW[i], true);
				vv_ij12 = computeVelocityVortex(elemsT[j], elemsW[i], true);

				vs_ij21 = computeVelocitySource(elemsW[j], elemsT[i], true);
				vv_ij21 = computeVelocityVortex(elemsW[j], elemsT[i], true);

				vs_ij22 = computeVelocitySource(elemsT[j], elemsT[i], false);
				vv_ij22 = computeVelocityVortex(elemsT[j], elemsT[i], false);

				A[i][j]               = scalar2(elemsW[i].nver, vs_ij11);
				A[i][nelems+j]        = scalar2(elemsW[i].nver, vs_ij12);
				A[nelems+i][j]        = scalar2(elemsT[i].nver, vs_ij21);
				A[nelems+i][nelems+j] = scalar2(elemsT[i].nver, vs_ij22);

				A[i][2*nelems]            = A[i][2*nelems] + scalar2(elemsW[i].nver, vv_ij11);
				A[i][2*nelems+1]          = A[i][2*nelems+1] + scalar2(elemsW[i].nver, vv_ij12);
				A[nelems+i][2*nelems]     = A[nelems+i][2*nelems] + scalar2(elemsT[i].nver, vv_ij21);
				A[nelems+i][2*nelems+1]   = A[nelems+i][2*nelems+1] + scalar2(elemsT[i].nver, vv_ij22);
				if (i == 0 || i == nelems - 1) {
					A[2*nelems][j]          = A[2*nelems][j] + scalar2(elemsW[i].tver, vs_ij11);
					A[2*nelems][2*nelems]   = A[2*nelems][2*nelems] + scalar2(elemsW[i].tver, vv_ij11);
					A[2*nelems][nelems+j]   = A[2*nelems][nelems+j] + scalar2(elemsW[i].tver, vs_ij12);
					A[2*nelems][2*nelems+1] =  A[2*nelems][2*nelems+1] + scalar2(elemsW[i].tver, vv_ij12);

					A[2*nelems+1][j]          = A[2*nelems+1][j] + scalar2(elemsT[i].tver, vs_ij21);
					A[2*nelems+1][2*nelems]   = A[2*nelems+1][2*nelems] + scalar2(elemsT[i].tver, vv_ij21);
					A[2*nelems+1][nelems+j]   = A[2*nelems+1][nelems+j] + scalar2(elemsT[i].tver, vs_ij22);
					A[2*nelems+1][2*nelems+1] =  A[2*nelems+1][2*nelems+1] + scalar2(elemsT[i].tver, vv_ij22);
				}

			}
		}
		return A;
	}

	double[][] coefficientMatrixG(Elems[] elems, Elems[] elems_mirr, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] A = new double[nelems + n_te][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = mySum(computeVelocitySource(elems[j], elems[i], false),
				              computeVelocitySource(elems_mirr[nelems-1-j],elems[i], true));
				vv_ij = mySub(computeVelocityVortex(elems[j], elems[i], false),
							  computeVelocityVortex(elems_mirr[nelems-1-j],elems[i],true));

				A[i][j]      = scalar2(elems[i].nver, vs_ij);
				A[i][nelems] = A[i][nelems] + scalar2(elems[i].nver, vv_ij);
				if (i == 0 || i == nelems - 1) {
					A[nelems][j]      = A[nelems][j] + scalar2(elems[i].tver, vs_ij);
					A[nelems][nelems] = A[nelems][nelems] + scalar2(elems[i].tver, vv_ij);
				}

			}
		}
		return A;
	}

	double[][] coefficientMatrixGWT(Elems[] elemsW, Elems[] elemsT, Elems[] elemsW_mirr, Elems[] elemsT_mirr, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 2;

		double[][] A = new double[2*nelems + n_te][2*nelems + n_te];
		double[] vs_ij11, vv_ij11, vs_ij12, vv_ij12, vs_ij21, vv_ij21, vs_ij22, vv_ij22;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij11 = mySum(computeVelocitySource(elemsW[j], elemsW[i], false),
				                computeVelocitySource(elemsW_mirr[nelems-1-j],elemsW[i],true));
				vv_ij11 = mySub(computeVelocityVortex(elemsW[j], elemsW[i], false),
						        computeVelocityVortex(elemsW_mirr[nelems-1-j],elemsW[i],true));

				vs_ij12 = mySum(computeVelocitySource(elemsT[j], elemsW[i], true),
						        computeVelocitySource(elemsT_mirr[nelems-1-j],elemsW[i],true));
				vv_ij12 = mySub(computeVelocityVortex(elemsT[j], elemsW[i], true),
						        computeVelocityVortex(elemsT_mirr[nelems-1-j],elemsW[i],true));

				vs_ij21 = mySum(computeVelocitySource(elemsW[j], elemsT[i], true),
						        computeVelocitySource(elemsW_mirr[nelems-1-j],elemsT[i],true));
				vv_ij21 = mySub(computeVelocityVortex(elemsW[j], elemsT[i], true),
						        computeVelocityVortex(elemsW_mirr[nelems-1-j],elemsT[i],true));

				vs_ij22 = mySum(computeVelocitySource(elemsT[j], elemsT[i], false),
						        computeVelocitySource(elemsT_mirr[nelems-1-j],elemsT[i],true));
				vv_ij22 = mySub(computeVelocityVortex(elemsT[j], elemsT[i], false),
						        computeVelocityVortex(elemsT_mirr[nelems-1-j],elemsT[i],true));

				A[i][j]               = scalar2(elemsW[i].nver, vs_ij11);
				A[i][nelems+j]        = scalar2(elemsW[i].nver, vs_ij12);
				A[nelems+i][j]        = scalar2(elemsT[i].nver, vs_ij21);
				A[nelems+i][nelems+j] = scalar2(elemsT[i].nver, vs_ij22);

				A[i][2*nelems]            = A[i][2*nelems] + scalar2(elemsW[i].nver, vv_ij11);
				A[i][2*nelems+1]          = A[i][2*nelems+1] + scalar2(elemsW[i].nver, vv_ij12);
				A[nelems+i][2*nelems]     = A[nelems+i][2*nelems] + scalar2(elemsT[i].nver, vv_ij21);
				A[nelems+i][2*nelems+1]   = A[nelems+i][2*nelems+1] + scalar2(elemsT[i].nver, vv_ij22);
				if (i == 0 || i == nelems - 1) {
					A[2*nelems][j]          = A[2*nelems][j] + scalar2(elemsW[i].tver, vs_ij11);
					A[2*nelems][2*nelems]   = A[2*nelems][2*nelems] + scalar2(elemsW[i].tver, vv_ij11);
					A[2*nelems][nelems+j]   = A[2*nelems][nelems+j] + scalar2(elemsW[i].tver, vs_ij12);
					A[2*nelems][2*nelems+1] =  A[2*nelems][2*nelems+1] + scalar2(elemsW[i].tver, vv_ij12);

					A[2*nelems+1][j]          = A[2*nelems+1][j] + scalar2(elemsT[i].tver, vs_ij21);
					A[2*nelems+1][2*nelems]   = A[2*nelems+1][2*nelems] + scalar2(elemsT[i].tver, vv_ij21);
					A[2*nelems+1][nelems+j]   = A[2*nelems+1][nelems+j] + scalar2(elemsT[i].tver, vs_ij22);
					A[2*nelems+1][2*nelems+1] =  A[2*nelems+1][2*nelems+1] + scalar2(elemsT[i].tver, vv_ij22);
				}

			}
		}
		return A;
	}

	// The following method is the correspondent of a part of the Matlab function build_lin_sys
	// and returns the array of the known terms.

	double[] coefficientArray(FreeStream freeStream, Elems[] elems, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[] b = new double[nelems + n_te];

		for (int i = 0; i < nelems; i++) {
			b[i] = -scalar2(elems[i].nver, freeStream.vvec);

			if (i == 0 || i == (nelems - 1)) {
				b[nelems] = b[nelems] - scalar2(elems[i].tver, freeStream.vvec);
			}
		}

		return b;
	}

	double[] coefficientArrayWT(FreeStream freeStream, Elems[] elems, Elems[] elemsT, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 2;

		double[] b = new double[2*nelems + n_te];

		for (int i = 0; i < nelems; i++) {
			b[i]        = -scalar2(elems[i].nver, freeStream.vvec);
			b[nelems+i] = -scalar2(elemsT[i].nver, freeStream.vvec);

			if (i == 0 || i == (nelems - 1)) {
				b[2*nelems]   = b[2*nelems] - scalar2(elems[i].tver, freeStream.vvec);
				b[2*nelems+1] = b[2*nelems+1] - scalar2(elemsT[i].tver, freeStream.vvec);
			}
		}

		return b;
	}

	// The following method is the correspondent of a part of the Matlab function build_lin_sys
	// and returns the auxiliary matrix Av.

	double[][] onBodyVMatrix(Elems[] elems, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] Av = new double[nelems][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = computeVelocitySource(elems[j], elems[i],false);
				vv_ij = computeVelocityVortex(elems[j], elems[i],false);

				Av[i][j] = vs_ij[1];
				Av[i][nelems] = Av[i][nelems] + vv_ij[1];
			}
		}
		return Av;
	}

	double[][] onBodyVMatrixWT(Elems[] elemsW, Elems[] elemsT, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 2;

		double[][] Av = new double[2*nelems][2*nelems + n_te];
		double[] vs_ij11, vv_ij11, vs_ij12, vv_ij12, vs_ij21, vv_ij21, vs_ij22, vv_ij22;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij11 = computeVelocitySource(elemsW[j], elemsW[i], false);
				vv_ij11 = computeVelocityVortex(elemsW[j], elemsW[i], false);

				vs_ij12 = computeVelocitySource(elemsT[j], elemsW[i], true);
				vv_ij12 = computeVelocityVortex(elemsT[j], elemsW[i], true);

				vs_ij21 = computeVelocitySource(elemsW[j], elemsT[i], true);
				vv_ij21 = computeVelocityVortex(elemsW[j], elemsT[i], true);

				vs_ij22 = computeVelocitySource(elemsT[j], elemsT[i], false);
				vv_ij22 = computeVelocityVortex(elemsT[j], elemsT[i], false);

				Av[i][j]                 = vs_ij11[1];
				Av[i][2*nelems]          = Av[i][2*nelems] + vv_ij11[1];
				Av[i][nelems+j]          = vs_ij12[1];
				Av[i][2*nelems+1]        = Av[i][2*nelems+1] + vv_ij12[1];
				Av[nelems+i][j]          = vs_ij21[1];
				Av[nelems+i][2*nelems]   = Av[nelems+i][2*nelems] + vv_ij21[1];
				Av[nelems+i][nelems+j]   = vs_ij22[1];
				Av[nelems+i][2*nelems+1] = Av[nelems+i][2*nelems+1] + vv_ij22[1];
			}
		}
		return Av;
	}

	double[][] onBodyVMatrixG(Elems[] elems, Elems[] elems_mirr, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] Av = new double[nelems][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = mySum(computeVelocitySource(elems[j], elems[i], false),
						computeVelocitySource(elems_mirr[nelems-1-j],elems[i], true));
				vv_ij = mySub(computeVelocityVortex(elems[j], elems[i], false),
						computeVelocityVortex(elems_mirr[nelems-1-j],elems[i],true));

				Av[i][j] = vs_ij[1];
				Av[i][nelems] = Av[i][nelems] + vv_ij[1];
			}
		}
		return Av;
	}

	double[][] onBodyVMatrixGWT(Elems[] elemsW, Elems[] elemsT, Elems[] elemsW_mirr, Elems[] elemsT_mirr, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 2;

		double[][] Av = new double[2*nelems][2*nelems + n_te];
		double[] vs_ij11, vv_ij11, vs_ij12, vv_ij12, vs_ij21, vv_ij21, vs_ij22, vv_ij22;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij11 = mySum(computeVelocitySource(elemsW[j], elemsW[i], false),
						computeVelocitySource(elemsW_mirr[nelems-1-j],elemsW[i],true));
				vv_ij11 = mySub(computeVelocityVortex(elemsW[j], elemsW[i], false),
						computeVelocityVortex(elemsW_mirr[nelems-1-j],elemsW[i],true));

				vs_ij12 = mySum(computeVelocitySource(elemsT[j], elemsW[i], true),
						computeVelocitySource(elemsT_mirr[nelems-1-j],elemsW[i],true));
				vv_ij12 = mySub(computeVelocityVortex(elemsT[j], elemsW[i], true),
						computeVelocityVortex(elemsT_mirr[nelems-1-j],elemsW[i],true));

				vs_ij21 = mySum(computeVelocitySource(elemsW[j], elemsT[i], true),
						computeVelocitySource(elemsW_mirr[nelems-1-j],elemsT[i],true));
				vv_ij21 = mySub(computeVelocityVortex(elemsW[j], elemsT[i], true),
						computeVelocityVortex(elemsW_mirr[nelems-1-j],elemsT[i],true));

				vs_ij22 = mySum(computeVelocitySource(elemsT[j], elemsT[i], false),
						computeVelocitySource(elemsT_mirr[nelems-1-j],elemsT[i],true));
				vv_ij22 = mySub(computeVelocityVortex(elemsT[j], elemsT[i], false),
						computeVelocityVortex(elemsT_mirr[nelems-1-j],elemsT[i],true));

				Av[i][j]                 = vs_ij11[1];
				Av[i][2*nelems]          = Av[i][2*nelems] + vv_ij11[1];
				Av[i][nelems+j]          = vs_ij12[1];
				Av[i][2*nelems+1]        = Av[i][2*nelems+1] + vv_ij12[1];
				Av[nelems+i][j]          = vs_ij21[1];
				Av[nelems+i][2*nelems]   = Av[nelems+i][2*nelems] + vv_ij21[1];
				Av[nelems+i][nelems+j]   = vs_ij22[1];
				Av[nelems+i][2*nelems+1] = Av[nelems+i][2*nelems+1] + vv_ij22[1];
			}
		}
		return Av;
	}

	// The following method is the correspondent of a part of the Matlab function build_lin_sys
	// and returns the auxiliary matrix Au.

	double[][] onBodyUMatrix(Elems[] elems, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] Au = new double[nelems][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = computeVelocitySource(elems[j], elems[i], false);
				vv_ij = computeVelocityVortex(elems[j], elems[i], false);

				Au[i][j] = vs_ij[0];
				Au[i][nelems] = Au[i][nelems] + vv_ij[0];
			}
		}
		return Au;
	}

	double[][] onBodyUMatrixWT(Elems[] elemsW, Elems[] elemsT, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 2;

		double[][] Au = new double[2*nelems][2*nelems + n_te];
		double[] vs_ij11, vv_ij11, vs_ij12, vv_ij12, vs_ij21, vv_ij21, vs_ij22, vv_ij22;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij11 = computeVelocitySource(elemsW[j], elemsW[i], false);
				vv_ij11 = computeVelocityVortex(elemsW[j], elemsW[i], false);

				vs_ij12 = computeVelocitySource(elemsT[j], elemsW[i], true);
				vv_ij12 = computeVelocityVortex(elemsT[j], elemsW[i], true);

				vs_ij21 = computeVelocitySource(elemsW[j], elemsT[i], true);
				vv_ij21 = computeVelocityVortex(elemsW[j], elemsT[i], true);

				vs_ij22 = computeVelocitySource(elemsT[j], elemsT[i], false);
				vv_ij22 = computeVelocityVortex(elemsT[j], elemsT[i], false);

				Au[i][j]                 = vs_ij11[0];
				Au[i][2*nelems]          = Au[i][2*nelems] + vv_ij11[0];
				Au[i][nelems+j]          = vs_ij12[0];
				Au[i][2*nelems+1]        = Au[i][2*nelems+1] + vv_ij12[0];
				Au[nelems+i][j]          = vs_ij21[0];
				Au[nelems+i][2*nelems]   = Au[nelems+i][2*nelems] + vv_ij21[0];
				Au[nelems+i][nelems+j]   = vs_ij22[0];
				Au[nelems+i][2*nelems+1] = Au[nelems+i][2*nelems+1] + vv_ij22[0];
			}
		}
		return Au;
	}

	double[][] onBodyUMatrixG(Elems[] elems, Elems[] elems_mirr, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] Au = new double[nelems][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = mySum(computeVelocitySource(elems[j], elems[i], false),
						computeVelocitySource(elems_mirr[nelems-1-j],elems[i], true));
				vv_ij = mySub(computeVelocityVortex(elems[j], elems[i], false),
						computeVelocityVortex(elems_mirr[nelems-1-j],elems[i],true));

				Au[i][j] = vs_ij[0];
				Au[i][nelems] = Au[i][nelems] + vv_ij[0];
			}
		}
		return Au;
	}

	double[][] onBodyUMatrixGWT(Elems[] elemsW, Elems[] elemsT, Elems[] elemsW_mirr, Elems[] elemsT_mirr, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 2;

		double[][] Au = new double[2*nelems][2*nelems + n_te];
		double[] vs_ij11, vv_ij11, vs_ij12, vv_ij12, vs_ij21, vv_ij21, vs_ij22, vv_ij22;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij11 = mySum(computeVelocitySource(elemsW[j], elemsW[i], false),
						computeVelocitySource(elemsW_mirr[nelems-1-j],elemsW[i],true));
				vv_ij11 = mySub(computeVelocityVortex(elemsW[j], elemsW[i], false),
						computeVelocityVortex(elemsW_mirr[nelems-1-j],elemsW[i],true));

				vs_ij12 = mySum(computeVelocitySource(elemsT[j], elemsW[i], true),
						computeVelocitySource(elemsT_mirr[nelems-1-j],elemsW[i],true));
				vv_ij12 = mySub(computeVelocityVortex(elemsT[j], elemsW[i], true),
						computeVelocityVortex(elemsT_mirr[nelems-1-j],elemsW[i],true));

				vs_ij21 = mySum(computeVelocitySource(elemsW[j], elemsT[i], true),
						computeVelocitySource(elemsW_mirr[nelems-1-j],elemsT[i],true));
				vv_ij21 = mySub(computeVelocityVortex(elemsW[j], elemsT[i], true),
						computeVelocityVortex(elemsW_mirr[nelems-1-j],elemsT[i],true));

				vs_ij22 = mySum(computeVelocitySource(elemsT[j], elemsT[i], false),
						computeVelocitySource(elemsT_mirr[nelems-1-j],elemsT[i],true));
				vv_ij22 = mySub(computeVelocityVortex(elemsT[j], elemsT[i], false),
						computeVelocityVortex(elemsT_mirr[nelems-1-j],elemsT[i],true));

				Au[i][j]                 = vs_ij11[0];
				Au[i][2*nelems]          = Au[i][2*nelems] + vv_ij11[0];
				Au[i][nelems+j]          = vs_ij12[0];
				Au[i][2*nelems+1]        = Au[i][2*nelems+1] + vv_ij12[0];
				Au[nelems+i][j]          = vs_ij21[0];
				Au[nelems+i][2*nelems]   = Au[nelems+i][2*nelems] + vv_ij21[0];
				Au[nelems+i][nelems+j]   = vs_ij22[0];
				Au[nelems+i][2*nelems+1] = Au[nelems+i][2*nelems+1] + vv_ij22[0];
			}
		}
		return Au;
	}



	// This method returns the value of the cross product between the two bidimensional given arrays
	// The result is used as a scalar number, even if mathematically should be threaten as an array
	// which is orthogonal to the plane shared by the two entry arrays.

	double cross(double[] a, double[] b) {
		return a[0] * b[1] - a[1] * b[0];
	}

	// This method returns the value of the scalar product between the two bidimensional given arrays.

	double scalar2(double[] a, double[] b) {
		return a[0] * b[0] + a[1] * b[1];
	}

	double[] mySum(double[] a, double[] b) {
		double[] c = new double[a.length];
		for(int i = 0; i < a.length; i++)
			c[i] = a[i] + b[i];
		return c;
	}

	double[] mySub(double[] a, double[] b) {
		double[] c = new double[a.length];
		for(int i = 0; i < a.length; i++)
			c[i] = a[i] - b[i];
		return c;
	}

	double[] reverseArray(double[] a, int n){
		double[] b = new double[n];
		int j = n;
		for(int i = 0; i<n; i++){
			b[j-1] = a[i];
			j = j-1;
		}
		return b;
	}

	// The linspace method is the equivalent of the linspace function in Matlab.

	public static double[] linspace(double min, double max, int points) {
		double[] d = new double[points];
		for (int i = 0; i < points; i++) {
			d[i] = min + i * (max - min) / (points - 1);
		}
		return d;
	}

	//Numerical Integration for obtaining the alfa0 in the thin wall theory, using a Gauss Legendre
	//method with 6 nodes and 6 weights.
	// For NACA 5 digits only, which need a numerical calculus of the integral.

	double gauss_legendre_integration(double a, double b, double c, double d, double a2,
										double a1, double a0, double b0, int check){
		double[] ybar = new double[]{-0.9324695142031520278123, -0.661209386466264513661,
									 -0.2386191860831969086305, 0.238619186083196908631,
				                     0.661209386466264513661, 0.9324695142031520278123};
		double[] alfabar = new double[]{0.1713244923791703450403, 0.3607615730481386075698,
										0.4679139345726910473899, 0.46791393457269104739,
										0.3607615730481386075698, 0.1713244923791703450403};
		double[] y = new double[6];
		double[] alfa = new double[6];
		double I1 = 0;
		double I2 = 0;
		for(int i = 0; i<6; i++){
			y[i] = (a+b)/2+(b-a)/2*ybar[i];
			alfa[i] = (b-a)/2*alfabar[i];
			I1 = I1 + alfa[i]*f1x(y[i], a2, a1, a0, check);
			y[i] = (c+d)/2+(d-c)/2*ybar[i];
			alfa[i] = (d-c)/2*alfabar[i];
			I2 = I2 + alfa[i]*f2x(y[i], b0, check);
		}
		return I1+I2;
	}

	// Function which need to be integrated in order to obtain alfa0, beta0 and alfa_theodorsen in
	// the thin wall theory. The integration will be computed between -0.5 and q-0.5;
	// For NACA 5 digits only, which need a numerical calculus of the integral.

	double f1x(double x, double a2, double a1, double a0, int check){
		if(check == 0)
			return Math.sqrt((0.5+x)/(0.5-x))*(a2*Math.pow(x,2)+a1*x+a0);
		else if(check == 1)
			return Math.sqrt(0.25-Math.pow(x,2))*(a2*Math.pow(x,2)+a1*x+a0);
		else
			return (a2*Math.pow(x,2)+a1*x+a0)/Math.sqrt(0.25-Math.pow(x,2));
	}

	// Function which need to be integrated in order to obtain alfa0, beta0 and alfa_theodorsen in
	// the thin wall theory. The integration will be computed between q-0.5 and 0.5;
	// For NACA 5 digits only, which need a numerical calculus of the integral.

	double f2x(double x, double b0, int check){
		if(check == 0)
			return Math.sqrt((0.5+x)/(0.5-x))*(b0);
		else if(check == 1)
			return Math.sqrt(0.25-Math.pow(x,2))*(b0);
		else
			return (b0)/Math.sqrt(0.25-Math.pow(x,2));
	}

	// The following methods are used in order to apply the Thwaites correction to the boundary layer.

	double thwaites_H(double lambda){
		double H;
		if(lambda <= -0.10)
			H = (2.088+0.0731/(-0.12+0.14));
		else if(-0.10<lambda && lambda <= 0.0)
			H = (2.088+0.0731/(lambda+0.14));
		else if(0.0<lambda && lambda<=0.25)
			H = (2.61-3.75*lambda-5.24*lambda*lambda);
		else
			H = (2.61-3.75*0.25+5.24*0.25*0.25);

		return H;
	}

	double thwaites_ell(double lambda){
		double ell;
		if(lambda <= -0.10)
			ell = (0.22+1.402*(-0.1)+0.018*(-0.1)/(-0.1+0.107));
		else if(-0.10<lambda && lambda <= 0.0)
			ell = (0.22+1.402*lambda+0.018*lambda/(lambda+0.107));
		else if(0.0<lambda && lambda<=0.25)
			ell = (0.22+1.57*lambda-1.8*lambda*lambda);
		else
			ell = (0.22+1.57*0.25-1.8*0.25*0.25);

		return ell;
	}

	double head_H1toH(double H1){
		double H;
		if(H1 < 3.3)
			H = 3.0;
		else if(3.3<=H1 && H1<5.3)
			H = 0.6778+1.1536*Math.pow((H1-3.3),(-0.326));
		else
			H = 1.1000+0.8600*Math.pow((H1-3.3),(-0.777));
		return H;
	}

	double head_HtoH1(double H){
		double H1;
		if(H<=1.6)
			H1 = 3.3 + 0.8234 * Math.pow((H - 1.1000),-1.287);
		else
			H1 = 3.3 + 1.5501 * Math.pow((H - 0.6778),-3.064);
		return H1;
	}

	double ludweig_tillman_cf(double H, double RetTheta){
		return 0.246*Math.pow(10.0,(-0.678*H))*Math.pow(RetTheta,-0.268);
	}
}

