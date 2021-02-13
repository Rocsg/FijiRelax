package com.vitimage.fijirelax.curvefit;

import lma.ArrayConverter;
import lma.LMAFunction;
import lma.LMAMatrix;
import lma.LMAMultiDimFunction;
import  lma.ArrayConverter.SeparatedData;
import lma.LMAMatrix.InvertException;
import  lma.implementations.JAMAMatrix;
import lma.implementations.LMA;

import java.util.Arrays;


/**
 * A class which implements the <i>Levenberg-Marquardt Algorithm</i>
 * (LMA) fit for non-linear, multidimensional parameter space
 * for any multidimensional fit function.
 * <p>
 * 
 * The algorithm is described in <i>Numerical Recipes in FORTRAN</i>,
 * 2nd edition, p. 676-679, ISBN 0-521-43064X, 1992 and also
 * <a href="http://www.nrbook.com/b/bookfpdf/f15-5.pdf">here</a> as a pdf file.
 * <p>
 * 
 * The matrix (<code>LMAMatrix</code>) class used in the fit is an interface, so you can use your
 * favourite implementation. This package uses <code>Matrix</code> from JAMA-math libraries,
 * but feel free to use anything you want. Note that you have to implement
 * the actual model function and its partial derivates as <code>LMAFunction</code>
 * or <code>LMAMultiDimFunction</code> before making the fit.
 * <p>
 * 
 * Note that there are <i>three</i> different ways to input the data points.
 * Read the documentation for each constructor carefully.
 * 
 * @author Janne Holopainen (jaolho@utu.fi, tojotamies@gmail.com)
 * @version 1.2, 24.04.2007
 * 
 * The algorithm is free for non-commercial use. 
 * 
 */
public class MRLMADual extends LMA{
	/** Set true to print details while fitting. */
	public boolean verbose = false;
	/** 
	 * The model function to be fitted, y = y(x[], a[]),
	 * where <code>x[]</code> the array of x-values and <code>a</code>
	 * is the array of fit parameters.
	 */
	public LMAMultiDimFunction function;
	/** 
	 * The array of fit parameters (a.k.a, the a-vector).
	 */
	public double[] parameters;
	/** 
	 * Measured y-data points for which the model function is to be fitted,
	 * yDataPoints[j] = y(xDataPoints[j], a[]).
	 */
	public double yDataPoints[];
	/** 
	 * Measured x-data point arrays for which the model function is to be fitted,
	 * yDataPoints[j] = y(xDataPoints[j], a[]).
	 * xDataPoints.length must be equal to yDataPoints.length and
	 * xDataPoints[].length must equal to the fit function's dimension.
	 */
	public double xDataPoints[][];	
	/** 
	 * Weights for each data point. The merit function is:
	 * chi2 = Sum[(y_i - y(x_i;a))^2 * w_i].
	 * For gaussian errors in datapoints, set w_i = (sigma_i)^-2.
	 */
	public double[] weights;
	
	public LMAMatrix alpha;
	public double[] beta;
	public double[] da;
	public double lambda = 0.001;
	public double lambdaFactor = 10;
	public double incrementedChi2;
	public double[] incrementedParameters;
	public int iterationCount;
	public double chi2;
	
	// default end conditions
	public double minDeltaChi2 = 1e-30;
	public int maxIterations = 100;
	private float[][] parameterBounds;
	
	
	
	/**
	 * One dimensional convenience constructor for LMAFunction.
	 * You can also implement the same function using LMAMultiDimFunction.
	 * <p>
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of data points, M is the number of fit parameters.
	 * Call <code>fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Must be able to take M input parameters.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in an array, <code>double[0 = x, 1 = y][point index]</code>.
	 * Size must be <code>double[2][N]</code>.
	 */
	public MRLMADual(final LMAFunction function, double[] parameters, double[][] dataPoints, float[][]bounds) {
		super(function, parameters, dataPoints, function.constructWeights(dataPoints));
		this.parameterBounds=bounds;
	}
	
	/**
	 * One dimensional convenience constructor for LMAFunction.
	 * You can also implement the same function using LMAMultiDimFunction.
	 * <p>
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of data points, M is the number of fit parameters.
	 * Call <code>fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Must be able to take M input parameters.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in an array, <code>double[0 = x, 1 = y][point index]</code>.
	 * Size must be <code>double[2][N]</code>.
	 */
	public MRLMADual(final LMAFunction function, double[] parameters, double[][] dataPoints, double[] weights, float[][]bounds) {
		super (
			// convert LMAFunction to LMAMultiDimFunction
			new LMAMultiDimFunction() { 
				private LMAFunction f = function;
				@Override
				public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
					return f.getPartialDerivate(x[0], a, parameterIndex);
				}
				@Override
				public double getY(double[] x, double[] a) {
					return f.getY(x[0], a);
				}
			}, 
			parameters,
			dataPoints[1], // y-data
			ArrayConverter.transpose(dataPoints[0]), // x-data
			weights,
			new JAMAMatrix(parameters.length, parameters.length)
		);		
		this.parameterBounds=bounds;

	}
	
	/**
	 * One dimensional convenience constructor for LMAFunction.
	 * You can also implement the same function using LMAMultiDimFunction.
	 * <p>
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of data points, M is the number of fit parameters.
	 * Call <code>fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Must be able to take M input parameters.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in an array, <code>float[0 = x, 1 = y][point index]</code>.
	 * Size must be <code>float[2][N]</code>.
	 */
	public MRLMADual(final LMAFunction function, float[] parameters, float[][] dataPoints, float[][]bounds) {
		super(
			function,
			ArrayConverter.asDoubleArray(parameters),
			ArrayConverter.asDoubleArray(dataPoints)
		);
		this.parameterBounds=bounds;
	}
	
	/**
	 * One dimensional convenience constructor for LMAFunction.
	 * You can also implement the same function using LMAMultiDimFunction.
	 * <p>
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of data points, M is the number of fit parameters.
	 * Call <code>fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Must be able to take M input parameters.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in an array, <code>float[0 = x, 1 = y][point index]</code>.
	 * Size must be <code>float[2][N]</code>.
	 * @param weights The weights, normally given as: <code>weights[i] = 1 / sigma_i^2</code>.
	 * If you have a bad data point, set its weight to zero.
	 * If the given array is null, a new array is created with all elements set to 1.
	 */
	public MRLMADual(final LMAFunction function, float[] parameters, float[][] dataPoints, float[] weights, float[][]bounds) {
		super(
			function,
			ArrayConverter.asDoubleArray(parameters),
			ArrayConverter.asDoubleArray(dataPoints),
			ArrayConverter.asDoubleArray(weights)
		);
		this.parameterBounds=bounds;

	}
	
	/**
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of y-data points, K is the dimension of the fit function and
	 * M is the number of fit parameters. Call <code>this.fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Input parameter sizes K and M.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in two dimensional array where each array, dataPoints[i],
	 * contains one y-value followed by the corresponding x-array values.
	 * I.e., the arrays should look like this:
	 * <p>
	 * dataPoints[0] = y0 x00 x01 x02 ... x0[K-1]<br>
	 * dataPoints[1] = y1 x10 x11 x12 ... x1[K-1]<br>
	 * ...<br>
	 * dataPoints[N] = yN xN0 xN1 xN2 ... x[N-1][K-1]
	 */
	public MRLMADual(LMAMultiDimFunction function, float[] parameters, float[][] dataPoints, float[][]bounds) {
		super (
			function, 
			ArrayConverter.asDoubleArray(parameters),
			ArrayConverter.asDoubleArray(dataPoints),
			function.constructWeights(ArrayConverter.asDoubleArray(dataPoints)),
			new JAMAMatrix(parameters.length, parameters.length)
		);
		this.parameterBounds=bounds;

	}
	
	/**
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of y-data points, K is the dimension of the fit function and
	 * M is the number of fit parameters. Call <code>this.fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Input parameter sizes K and M.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in two dimensional array where each array, dataPoints[i],
	 * contains one y-value followed by the corresponding x-array values.
	 * I.e., the arrays should look like this:
	 * <p>
	 * dataPoints[0] = y0 x00 x01 x02 ... x0[K-1]<br>
	 * dataPoints[1] = y1 x10 x11 x12 ... x1[K-1]<br>
	 * ...<br>
	 * dataPoints[N] = yN xN0 xN1 xN2 ... x[N-1][K-1]
	 */
	public MRLMADual(LMAMultiDimFunction function, double[] parameters, double[][] dataPoints, float[][]bounds) {
		super (
			function, 
			parameters,
			dataPoints,
			function.constructWeights(dataPoints),
			new JAMAMatrix(parameters.length, parameters.length)
		);
		this.parameterBounds=bounds;

	}
	
	/**
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of y-data points, K is the dimension of the fit function and
	 * M is the number of fit parameters. Call <code>this.fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Input parameter sizes K and M.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param yDataPoints The y-data points in an array.
	 * @param xDataPoints The x-data points for each y data point, double[y-index][x-index]
	 */
	public MRLMADual(LMAMultiDimFunction function, double[] parameters, float[] yDataPoints, float[][] xDataPoints, float[][]bounds) {
		super (
			function, 
			parameters,
			ArrayConverter.asDoubleArray(yDataPoints),
			ArrayConverter.asDoubleArray(xDataPoints),
			function.constructWeights(ArrayConverter.combineMultiDimDataPoints(yDataPoints, xDataPoints)),
			new JAMAMatrix(parameters.length, parameters.length)
		);
		this.parameterBounds=bounds;

	}
	
	/**
	 * Initiates the fit with function constructed weights and a JAMA matrix.
	 * N is the number of y-data points, K is the dimension of the fit function and
	 * M is the number of fit parameters. Call <code>this.fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Input parameter sizes K and M.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param yDataPoints The y-data points in an array.
	 * @param xDataPoints The x-data points for each y data point, double[y-index][x-index]
	 */
	public MRLMADual(LMAMultiDimFunction function, double[] parameters, double[] yDataPoints, double[][] xDataPoints, float[][]bounds) {
		super (
			function, 
			parameters,
			yDataPoints,
			xDataPoints,
			function.constructWeights(ArrayConverter.combineMultiDimDataPoints(yDataPoints, xDataPoints)),
			new JAMAMatrix(parameters.length, parameters.length)
		);
		this.parameterBounds=bounds;

	}
	
	/**
	 * Initiates the fit. N is the number of y-data points, K is the dimension of the fit function and
	 * M is the number of fit parameters. Call <code>this.fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Input parameter sizes K and M.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in two dimensional array where each array, dataPoints[i],
	 * contains one y-value followed by the corresponding x-array values.
	 * I.e., the arrays should look like this:
	 * <p>
	 * dataPoints[0] = y0 x00 x01 x02 ... x0[K-1]<br>
	 * dataPoints[1] = y1 x10 x11 x12 ... x1[K-1]<br>
	 * ...<br>
	 * dataPoints[N] = yN xN0 xN1 xN2 ... x[N-1][K-1]
	 * <p>
	 * @param weights The weights, normally given as: <code>weights[i] = 1 / sigma_i^2</code>.
	 * If you have a bad data point, set its weight to zero.
	 * If the given array is null, a new array is created with all elements set to 1.
	 * @param alpha An LMAMatrix instance. Must be initiated to (M x M) size.
	 */
	public MRLMADual(LMAMultiDimFunction function, float[] parameters, float[][] dataPoints, float[] weights, LMAMatrix alpha, float[][]bounds) {
		super(
			function, 
			ArrayConverter.asDoubleArray(parameters),
			ArrayConverter.asDoubleArray(dataPoints),
			ArrayConverter.asDoubleArray(weights),
			alpha);
		this.parameterBounds=bounds;

	}	
	
	/**
	 * Initiates the fit. N is the number of y-data points, K is the dimension of the fit function and
	 * M is the number of fit parameters. Call <code>this.fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Input parameter sizes K and M.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param dataPoints The data points in two dimensional array where each array, dataPoints[i],
	 * contains one y-value followed by the corresponding x-array values.
	 * I.e., the arrays should look like this:
	 * <p>
	 * dataPoints[0] = y0 x00 x01 x02 ... x0[K-1]<br>
	 * dataPoints[1] = y1 x10 x11 x12 ... x1[K-1]<br>
	 * ...<br>
	 * dataPoints[N] = yN xN0 xN1 xN2 ... x[N-1][K-1]
	 * <p>
	 * @param weights The weights, normally given as: <code>weights[i] = 1 / sigma_i^2</code>.
	 * If you have a bad data point, set its weight to zero.
	 * If the given array is null, a new array is created with all elements set to 1.
	 * @param alpha An LMAMatrix instance. Must be initiated to (M x M) size.
	 */
	public MRLMADual(LMAMultiDimFunction function, double[] parameters, double[][] dataPoints, double[] weights, LMAMatrix alpha, float[][]bounds) {
		super(function,parameters,dataPoints,weights,alpha);
		SeparatedData s = ArrayConverter.separateMultiDimDataToXY(dataPoints);
		this.yDataPoints = s.yDataPoints;
		this.xDataPoints = s.xDataPoints;
		init(function, parameters, yDataPoints, xDataPoints, weights, alpha);
		this.parameterBounds=bounds;
	}
	

	/**
	 * Initiates the fit. N is the number of y-data points, K is the dimension of the fit function and
	 * M is the number of fit parameters. Call <code>this.fit()</code> to start the actual fitting.
	 * 
	 * @param function The model function to be fitted. Must be able to take M input parameters.
	 * @param parameters The initial guess for the fit parameters, length M.
	 * @param yDataPoints The y-data points in an array.
	 * @param xDataPoints The x-data points for each y data point, double[y-index][x-index]
	 * Size must be <code>double[N][K]</code>, where N is the number of measurements
	 * and K is the dimension of the fit function.
	 * @param weights The weights, normally given as: <code>weights[i] = 1 / sigma_i^2</code>.
	 * If you have a bad data point, set its weight to zero. If the given array is null,
	 * a new array is created with all elements set to 1.
	 * @param alpha An LMAMatrix instance. Must be initiated to (M x M) size.
	 */
	public MRLMADual(LMAMultiDimFunction function, double[] parameters, double[] yDataPoints, double[][] xDataPoints, double[] weights, LMAMatrix alpha, float[][]bounds) {
		super(function,parameters,xDataPoints,weights,alpha);
		init(function, parameters, yDataPoints, xDataPoints, weights, alpha);
		this.parameterBounds=bounds;

	}
	

}
