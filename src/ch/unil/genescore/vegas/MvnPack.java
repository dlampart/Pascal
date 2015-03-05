/*******************************************************************************
 * Copyright (c) 2015 IBM Corporation and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     IBM Corporation - initial API and implementation
 *******************************************************************************/
package ch.unil.genescore.vegas;

import com.sun.jna.Library;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

public interface MvnPack extends Library {

		/**
		 * See http://www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f
		 * @param n
		 * @param lower
		 * @param upper
		 * @param infin
		 * @param correl
		 * @param maxpts
		 * @param abseps
		 * @param releps
		 * @param error
		 * @param value
		 * @param inform
		 */
		void mvtdst_(IntByReference n,IntByReference df, double[] lower, double[] upper, int[] infin, double[] correl, double[] delta,
				IntByReference maxpts, DoubleByReference abseps, DoubleByReference releps, 
				DoubleByReference error, DoubleByReference value, IntByReference inform);
	
		/**
		 * See http://www.math.wsu.edu/faculty/genz/software/fort77/mvnexppack.f
		 * @param n
		 * @param lower
		 * @param upper
		 * @param infin
		 * @param correl
		 * @param maxpts
		 * @param abseps
		 * @param releps
		 * @param error
		 * @param value
		 * @param inform
		 */
	//	void mvnexp_(IntByReference n, double[] lower, double[] upper, int[] infin, double[] correl,
		//		IntByReference maxpts, DoubleByReference abseps, DoubleByReference releps,
			//	double[] error, double[] value, IntByReference inform);

		/**
		 * See http://www.math.wsu.edu/faculty/genz/software/fort77/mvnxpppack.f
		 * @param n
		 * @param lower
		 * @param upper
		 * @param infin
		 * @param correl
		 * @param maxpts
		 * @param abseps
		 * @param releps
		 * @param error
		 * @param value
		 * @param inform
		 */
	//	void mvnxpp_(IntByReference n, double[] lower, double[] upper, int[] infin, double[] correl,
		//		IntByReference maxpts, DoubleByReference abseps, DoubleByReference releps,
			//	double[] error, double[] value, IntByReference inform);

}


