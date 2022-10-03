/*
 * 
 */
package io.github.rocsg.fijirelax.lma.implementations;//Initially joalho.data.lma.implementations, see  https://zenodo.org/record/4281064

import io.github.rocsg.fijirelax.lma.LMAMatrix;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.NumberFormat;

import Jama.Matrix;

// TODO: Auto-generated Javadoc
/**
 * The Class JAMAMatrix.
 */
public class JAMAMatrix extends Matrix implements LMAMatrix {
	
	/** The Constant serialVersionUID. */
	private static final long serialVersionUID = -8925816623803983503L;

	/**
	 * Instantiates a new JAMA matrix.
	 *
	 * @param elements the elements
	 */
	public JAMAMatrix(double[][] elements) {
		super(elements);
	}
	
	/**
	 * Instantiates a new JAMA matrix.
	 *
	 * @param rows the rows
	 * @param cols the cols
	 */
	public JAMAMatrix(int rows, int cols) {
		super(rows, cols);
	}
	
	/**
	 * Invert.
	 *
	 * @throws InvertException the invert exception
	 */
	public void invert() throws LMAMatrix.InvertException {
		try {
			Matrix m = inverse();
			setMatrix(0, this.getRowDimension() - 1, 0, getColumnDimension() - 1, m);
		}
		catch (RuntimeException e) {
			StringWriter s = new StringWriter();
			PrintWriter p = new PrintWriter(s);
			p.println(e.getMessage());
			p.println("Inversion failed for matrix:");
			this.print(p, NumberFormat.getInstance(), 5);
			throw new LMAMatrix.InvertException(s.toString());
		}
	}

	/**
	 * Sets the element.
	 *
	 * @param row the row
	 * @param col the col
	 * @param value the value
	 */
	public void setElement(int row, int col, double value) {
		set(row, col, value);
	}

	/**
	 * Gets the element.
	 *
	 * @param row the row
	 * @param col the col
	 * @return the element
	 */
	public double getElement(int row, int col) {
		return get(row, col);
	}

	/**
	 * Multiply.
	 *
	 * @param vector the vector
	 * @param result the result
	 */
	public void multiply(double[] vector, double[] result) {
		for (int i = 0; i < this.getRowDimension(); i++) {
			result[i] = 0;
			for (int j = 0; j < this.getColumnDimension(); j++) {
				 result[i] += this.getElement(i, j) * vector[j];
			}
		}
	}
	
	/**
	 * The main method.
	 *
	 * @param args the arguments
	 */
	public static void main(String[] args) {
		StringWriter s = new StringWriter();
		PrintWriter out = new PrintWriter(s);
		out.println("jakkajaaa");
		System.out.println(s);
	}
}
