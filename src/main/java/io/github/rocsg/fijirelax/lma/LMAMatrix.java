/*
 * 
 */
package io.github.rocsg.fijirelax.lma;//Initially joalho.data.lma, see  https://zenodo.org/record/4281064

// TODO: Auto-generated Javadoc
/**
 * The matrix to be used in LMA.
 * Implement this to make LMA operational if you
 * don't or can't use jama or flanagan math libraries.
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>
 */
public interface LMAMatrix {
	
	/**
	 * The Class InvertException.
	 */
	public static class InvertException extends RuntimeException {
		
		/**
		 * Instantiates a new invert exception.
		 *
		 * @param message the message
		 */
		public InvertException(String message) {
			super(message);
		}
	}
	
	/**
	 * Inverts the matrix for solving linear equations for
	 * parameter increments.
	 *
	 * @throws InvertException the invert exception
	 */
	public void invert() throws InvertException;
	
	/**
	 * Set the value of a matrix element.
	 *
	 * @param row the row
	 * @param col the col
	 * @param value the value
	 */
	public void setElement(int row, int col, double value);
	
	/**
	 * Get the value of a matrix element.
	 *
	 * @param row the row
	 * @param col the col
	 * @return the element
	 */
	public double getElement(int row, int col);
	
	/**
	 * Multiplies this matrix with an array (result = this * vector).
	 * The lengths of the arrays must be equal to the number of rows in the matrix.
	 * @param vector The array to be multiplied with the matrix.
	 * @param result The result of the multiplication will be put here.
	 */
	public void multiply(double[] vector, double[] result);
}
