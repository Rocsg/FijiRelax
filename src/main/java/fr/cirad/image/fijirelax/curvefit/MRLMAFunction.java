package fr.cirad.image.fijirelax.curvefit;
import java.util.Arrays;

import lma.ArrayConverter;

/** 
 * Implement this for your fit function.
 * 
 * @author Janne Holopainen (jaolho@utu.fi, tojotamies@gmail.com)
 * @version 1.2, 16.11.2006
 */
public abstract class MRLMAFunction {
	/**
	 * @return The <i>y</i>-value of the function.
	 * @param x The <i>x</i>-value for which the <i>y</i>-value is calculated.
	 * @param a The fitting parameters. 
	 */
	public abstract double getY(double x, double[] a);
	
	public abstract double getY(double Tr,double Te, double[] a);

	/** 
	 * The method which gives the partial derivates used in the LMA fit.
	 * If you can't calculate the derivate, use a small <code>a</code>-step (e.g., <i>da</i> = 1e-20)
	 * and return <i>dy/da</i> at the given <i>x</i> for each fit parameter.
	 * @return The partial derivate of the function with respect to parameter <code>parameterIndex</code> at <i>x</i>.
	 * @param x The <i>x</i>-value for which the partial derivate is calculated.
	 * @param a The fitting parameters.
	 * @param parameterIndex The parameter index for which the partial derivate is calculated. 
	 */
	public abstract double getPartialDerivate(double x, double[] a, int parameterIndex);

	public abstract double getPartialDerivate(double Tr, double Te,double[] a, int parameterIndex);
	
	/** @return Calculated function values with the given x- and parameter-values. */
	public double[] generateData(double[] Tr,double[] Te, double[] a) {
		double[] result = new double[Tr.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = getY(Tr[i],Te[i], a);
		}
		return result;
	}
	
	/** @return Calculated function values with the given x- and parameter-values. */
	public float[] generateData(float[] Tr,float[] Te, double[] a) {
		float[] result = new float[Tr.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = (float) getY(Tr[i],Te[i], a);
		}
		return result;
	}
	
	/** The default weights-array constructor. Override for your purposes. */
	public double[] constructWeights(double[][] dataPoints) {
		double[] result = new double[dataPoints[0].length];
		Arrays.fill(result, 1);
		return result;
	}
	
	public float[] generateData(float[] Tr, float[] Te,float[] a) {
		return ArrayConverter.asFloatArray(generateData(ArrayConverter.asDoubleArray(Tr), ArrayConverter.asDoubleArray(Te),ArrayConverter.asDoubleArray(a)));
	}
	
}