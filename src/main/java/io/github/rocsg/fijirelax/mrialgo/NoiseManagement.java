/*
 * 
 */
package io.github.rocsg.fijirelax.mrialgo;

// TODO: Auto-generated Javadoc
/** This enum define cases of noise management. NOTHING means no noise is estimated, OFFSET means noise is estimated as an additional BIAS, and RICE means full Rice Noise estimation 
 * See the paper Fernandez et al. 2022 (in prep) for more information
 *
 * 	Copyright (C) 2022  <io.github.rocsg>
 *
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
*/
public enum NoiseManagement {
	
	/** The nothing. */
	NOTHING,
	
	/** The offset. */
	OFFSET,
	
	/** The rice. */
	RICE
}