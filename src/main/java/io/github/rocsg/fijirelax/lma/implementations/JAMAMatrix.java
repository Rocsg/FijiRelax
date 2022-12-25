package io.github.rocsg.fijirelax.lma.implementations;//Initially joalho.data.lma.implementations, see  https://zenodo.org/record/4281064

import io.github.rocsg.fijirelax.lma.LMAMatrix;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.NumberFormat;

import Jama.Matrix;

/**
 * The Class JAMAMatrix.
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>
 * @author Janne Holopainen (jaolhoATutu.fi, tojotamiesATgmail.com)
 * @version 1.0, 23.03.2007
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
