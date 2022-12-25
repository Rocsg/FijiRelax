/*
 * 
 */
package io.github.rocsg.fijirelax.mrialgo;

/** This enum define different common MRI imaging sequence. Their handling is done by HyperMap
 * 
 * @see io.github.rocsg.fijirelax.mrialgo.HyperMap
 * 	Copyright (C) 2022  io.github.rocsg
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
 *   along with this program.  If not, see https://www.gnu.org/licenses/
*/
public enum MRDataType {
	
	/** The t1seq. */
	T1SEQ,
	
	/** The t2seq. */
	T2SEQ,
	
	/** The t1t2seq. */
	T1T2SEQ,
	
	/** The tr600. */
	TR600,
	
	/** The tr1200. */
	TR1200,
	
	/** The tr2400. */
	TR2400,
	
	/** The tr4800. */
	TR4800,
	
	/** The tr10000. */
	TR10000,
	
	/** The t2seqhr. */
	T2SEQHR,
	
	/** The ge3d. */
	GE3D,
	
	/** The t1map. */
	T1MAP,
	
	/** The t2map. */
	T2MAP,
	
	/** The pdmap. */
	PDMAP,
	
	/** The maskmap. */
	MASKMAP,
	
	/** The other. */
	OTHER
}
