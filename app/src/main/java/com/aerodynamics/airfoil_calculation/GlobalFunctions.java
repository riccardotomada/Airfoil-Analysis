package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//Questa classe contiene tutte le funzioni condivise tra più activity						      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import java.lang.reflect.Array;

public class GlobalFunctions {

	// Tutti i metodi presenti, eccetto quelli commentati nello specifico, sono la trasposizione
	// in Java delle funzioni scritte in Matlab

	double[][] getCoordinates(Airfoil airfoil){

		String prova = "NACA";
		for(int i = 0; i < 4; i++){
			if (airfoil.airfoil_str.charAt(i) != prova.charAt(i)) {
				return new double[][]  {{0.0},{0.0}};
			}
		}

		double MC, POS, TH;
		char m,p;
		char[] th,mm;
		int coeff;

		if(airfoil.airfoil_str.length() == 8){
			m = airfoil.airfoil_str.charAt(4);
			p = airfoil.airfoil_str.charAt(5);
			th = new char[]{airfoil.airfoil_str.charAt(6), airfoil.airfoil_str.charAt(7)};
			MC = Double.parseDouble(String.valueOf(m));
			POS = Double.parseDouble(String.valueOf(p));
			TH = Double.parseDouble(String.valueOf(th));
			coeff = 0;
		}
		else{
			mm = new char[]{airfoil.airfoil_str.charAt(4), airfoil.airfoil_str.charAt(5),
			         		airfoil.airfoil_str.charAt(6)};
			th = new char[]{airfoil.airfoil_str.charAt(7), airfoil.airfoil_str.charAt(8)};
			coeff = Integer.parseInt(String.valueOf(mm));
			TH = Double.parseDouble(String.valueOf(th));
			MC = 0.0;
			POS = 0.0;
		}


		if(airfoil.airfoil_str.length() == 8)
		return setAirfoil4(MC, POS, TH, airfoil.chord, airfoil.nChordPanels,
				airfoil.refPoint[0], airfoil.refPoint[1], airfoil.xcRefPoint, airfoil.theta);
		else
			return setAirfoil5(coeff, TH, airfoil.chord, airfoil.nChordPanels,
					airfoil.refPoint[0], airfoil.refPoint[1], airfoil.xcRefPoint, airfoil.theta);
	}

	double[][] naca4digit(double MC, double POS, double TH, double c, int n) {

		double mc = MC / 100;
		double pos = POS / 10;
		double th = TH / 100;

		double[] xVector = linspace(0, c, n + 1);
		for (int i = 0; i < Array.getLength(xVector); i++)
			xVector[i] = c * (1 - Math.cos(0.5*Math.PI * xVector[i] / c));

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
			if (xVector[i] < pos)
				ycamber[i] = mc / (Math.pow(pos, 2)) * xVector[i] * (2 * pos - xVector[i]);
			else
				ycamber[i] = mc / (Math.pow((1 - pos), 2)) * (1 - 2 * pos + 2 * pos * xVector[i] - xVector[i] * xVector[i]);
		}

		double[] dycamber = new double[Array.getLength((xVector))];

		for (int i = 0; i < (n + 1); i++) {
			if (xVector[i] < pos)
				dycamber[i] = 2 * mc / Math.pow(pos, 2) * (pos - xVector[i]);
			else
				dycamber[i] = 2 * mc / Math.pow((1 - pos), 2) * (pos - xVector[i]);
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

		double t = TH/100;

		double[] xVector = linspace(0, c, n + 1);
		for (int i = 0; i < Array.getLength(xVector); i++)
			xVector[i] = c/2 * (1 - Math.cos(Math.PI * xVector[i] / c));

		double q,k;

		switch(coeff){
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
				if (xVector[i] <= q * c)
					ycamber[i] = c * (k / 6 * Math.pow(xVector[i] / c, 3) - 3 * q * Math.pow(xVector[i] / c, 2)
							+ q * q * (3 - q) * (xVector[i] / c));
				else
					ycamber[i] = c * k / 6 * q * q * q * (1 - (xVector[i] / c));
			}

			double[] dycamber = new double[Array.getLength((xVector))];

			for (int i = 0; i < (n + 1); i++) {
				if (xVector[i] <= q * c)
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
						   double ca, double theta){
		double[][] ab;

		ab = naca4digit(MC, POS, TH, c, n);

		double c1 = ca*c;

		double[] x = new double[2 * n + 1];
		double[] y = new double[2 * n + 1];

		for(int i = 0; i < Array.getLength(x); i++){
			x[i] = x0 + c1 + Math.cos(-theta)*(ab[0][i]-c1) + Math.sin(-theta)*ab[1][i];
			y[i] = y0 - Math.sin(-theta)*(ab[0][i] - c1) + Math.cos(-theta)*ab[1][i];
		}

		double[][] xy = new double[2][2 * n +1];
		for (int i = 0; i < Array.getLength(x); i++) {
			xy[0][i] = x[i];
			xy[1][i] = y[i];
		}

		return xy;
	}

	double[][] setAirfoil5(int coeff, double TH, double c, int n, double x0, double y0,
						   double ca, double theta){
		double[][] ab;

		ab = naca5digit(coeff, TH, c, n);

		double c1 = ca*c;

		double[] x = new double[2 * n + 1];
		double[] y = new double[2 * n + 1];

		for(int i = 0; i < Array.getLength(x); i++){
			x[i] = x0 + c1 + Math.cos(-theta)*(ab[0][i]-c1) + Math.sin(-theta)*ab[1][i];
			y[i] = y0 - Math.sin(-theta)*(ab[0][i] - c1) + Math.cos(-theta)*ab[1][i];
		}

		double[][] xy = new double[2][2 * n +1];
		for (int i = 0; i < Array.getLength(x); i++) {
			xy[0][i] = x[i];
			xy[1][i] = y[i];
		}

		return xy;
	}

	int[][] connectivityMatrix(double[] x) {
		int[][] conMat = new int[2][Array.getLength(x)-1];
		for(int i = 1; i < Array.getLength(x); i++){
			conMat[0][i-1] = i;
			conMat[1][i-1] = i+1;
		}
		return conMat;
	}

	double[] computeVelocityVortex(Elems elem_j, Elems elem_i){
		double[] r1 = new double[2];
		double[] r2 = new double[2];
		double[] vstar, v;
		double sij, bij, vcross, sinbij, cosbij;

		r1[0] = elem_i.cen[0] - elem_j.ver1[0];
		r1[1] = elem_i.cen[1] - elem_j.ver1[1];
		r2[0] = elem_i.cen[0] - elem_j.ver2[0];
		r2[1] = elem_i.cen[1] - elem_j.ver2[1];

		if(elem_i.id == elem_j.id){
			sij = 1.0;
			bij = Math.PI;
		}
		else{
			sij    = Math.sqrt(Math.pow(r2[0],2)+Math.pow(r2[1],2))/
					 Math.sqrt(Math.pow(r1[0],2)+Math.pow(r1[1],2));
			vcross = cross(r1,r2);
			sinbij = vcross/Math.sqrt(Math.pow(r1[0],2)+Math.pow(r1[1],2))/
					 Math.sqrt(Math.pow(r2[0],2)+Math.pow(r2[1],2));
			cosbij = scalar2(r1,r2)/Math.sqrt(Math.pow(r1[0],2)+Math.pow(r1[1],2))/
					 Math.sqrt(Math.pow(r2[0],2)+Math.pow(r2[1],2));
			bij    = Math.atan2(sinbij,cosbij);
		}
		vstar = new double[]{-bij / (2 * Math.PI), -Math.log(sij) / (2 * Math.PI)};

		v = new double[]{elem_j.tver[0]*vstar[0] + elem_j.nver[0]*vstar[1],
				elem_j.tver[1]*vstar[0] + elem_j.nver[1]*vstar[1]};
		return v;
	}

	double[] computeVelocitySource(Elems elem_j, Elems elem_i){
		double[] r1 = new double[2];
		double[] r2 = new double[2];
		double[] vstar, v;
		double sij, bij, vcross, sinbij, cosbij;

		r1[0] = elem_i.cen[0] - elem_j.ver1[0];
		r1[1] = elem_i.cen[1] - elem_j.ver1[1];
		r2[0] = elem_i.cen[0] - elem_j.ver2[0];
		r2[1] = elem_i.cen[1] - elem_j.ver2[1];

		if(elem_i.id == elem_j.id){
			sij = 1.0;
			bij = Math.PI;
		}
		else {
			sij    = Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2)) /
					 Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2));
			vcross = cross(r1, r2);
			sinbij = vcross / Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2)) /
					 Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2));
			cosbij = scalar2(r1, r2) / Math.sqrt(Math.pow(r1[0], 2) + Math.pow(r1[1], 2)) /
					 Math.sqrt(Math.pow(r2[0], 2) + Math.pow(r2[1], 2));
			bij    = Math.atan2(sinbij, cosbij);
		}


		vstar = new double[]{-Math.log(sij) / (2 * Math.PI), bij / (2 * Math.PI)};
		v     = new double[]{elem_j.tver[0]*vstar[0] + elem_j.nver[0]*vstar[1],
				elem_j.tver[1]*vstar[0] + elem_j.nver[1]*vstar[1]};
		return v;
	}

	// La function midpoint genera un vettore i cui elementi sono la media degli elementi corrispondenti
	// di due vettori in ingresso.
	double[] midpoint(double[] a, double[] b){
		double[] c = new double[Array.getLength(a)];
		for (int i = 0; i < Array.getLength(a); ++i) {
			c[i] = a[i] + b[i];
			c[i] = c[i]/2;
		}
		return c;
	}

	// La function len_norm calcola la norma di un vettore di dimensione 2 dato dalla differenza tra i due in
	// ingresso.
	double len_norm(double[] a, double[] b){
		double[] c = new double[2];
		for (int i = 0; i < 2; i++) {
			c[i] = a[i] - b[i];
		}
		return Math.sqrt(c[0]*c[0]+c[1]*c[1]);
	}

	double[] tver(double[] a, double[] b, double norm){
		double[] c = new double[2];
		for (int i = 0; i < 2; i++) {
			c[i] = (a[i] - b[i])/norm;
		}
		return c;
	}

	//Questo metodo è una parte dell function MatLab build_lin_sys e restituisce la matrice dei
	//coefficienti A

	double[][] coefficientMatrix(Elems[] elems, int nChordPanels){
		int nelems = 2*nChordPanels;
		int n_te = 1;

		double[][] A = new double[nelems+n_te][nelems+n_te];
		double[] vs_ij;
		double[] vv_ij;

		for(int i = 0; i < nelems; i++){
			for(int j = 0; j < nelems; j++) {
				vs_ij = computeVelocitySource(elems[j], elems[i]);
				vv_ij = computeVelocityVortex(elems[j], elems[i]);

				A[i][j] = scalar2(elems[i].nver, vs_ij);
				A[i][nelems] = A[i][nelems] + scalar2(elems[i].nver, vv_ij);
				if(i == 0 || i == nelems-1){
					A[nelems][j] = A[nelems][j] + scalar2(elems[i].tver,vs_ij);
					A[nelems][nelems] = A[nelems][nelems] + scalar2(elems[i].tver,vv_ij);
				}

			}
		}
		return A;
	}

	//Questo metodo è una parte dell function MatLab build_lin_sys e restituisce il vettore dei
	//termini noti b

	double[] coefficientArray(FreeStream freeStream, Elems[] elems, int nChordPanels){
		int nelems = 2*nChordPanels;
		int n_te = 1;

		double[] b = new double[nelems+n_te];

		for(int i = 0; i < nelems; i++){
			b[i] = scalar2(elems[i].nver,freeStream.vvec);
			b[i] = -b[i];

			if(i == 0 || i == (nelems-1)){
				b[nelems] = b[nelems] - scalar2(elems[i].tver, freeStream.vvec);
			}
		}

		return b;
	}

	//Questo metodo è una parte dell function MatLab build_lin_sys e restituisce la matrice
	//ausiliaria Av

	double[][] onBodyVMatrix(Elems[] elems, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] Av = new double[nelems][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = computeVelocitySource(elems[j], elems[i]);
				vv_ij = computeVelocityVortex(elems[j], elems[i]);

				Av[i][j] = vs_ij[1];
				Av[i][nelems] = Av[i][nelems] + vv_ij[1];
			}
		}
		return Av;
	}

	//Questo metodo è una parte dell function MatLab build_lin_sys e restituisce la matrice
	//ausiliaria Au

	double[][] onBodyUMatrix(Elems[] elems, int nChordPanels) {
		int nelems = 2 * nChordPanels;
		int n_te = 1;

		double[][] Au = new double[nelems][nelems + n_te];
		double[] vs_ij;
		double[] vv_ij;

		for (int i = 0; i < nelems; i++) {
			for (int j = 0; j < nelems; j++) {
				vs_ij = computeVelocitySource(elems[j], elems[i]);
				vv_ij = computeVelocityVortex(elems[j], elems[i]);

				Au[i][j] = vs_ij[0];
				Au[i][nelems] = Au[i][nelems] + vv_ij[0];
			}
		}
		return Au;
	}

	//Questo metodo restituisce il valore del prodotto vettoriale tra i due vettori in ingresso di dimensione 2
	//Il risultato viene trattato per quanto serve nel programma come uno scalare, anchè se in realtà
	//dovrebbe essere un vettore di direzione ortogonale al piano condiviso dai due vettori in ingresso

	double cross(double[] a, double[] b){
		return a[0]*b[1]-a[1]*b[0];
	}

	//Questo metodo restituisce il valore del prodotto scalare tra i due vettori in ingresso di dimensione 2

	double scalar2(double[] a, double[] b){
		return a[0]*b[0]+a[1]*b[1];
	}

	//Equivalente della funzione linspace di MatLab

	public static double[] linspace(double min, double max, int points) {
		double[] d = new double[points];
		for (int i = 0; i < points; i++) {
			d[i] = min + i * (max - min) / (points - 1);
		}
		return d;
	}

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

