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
package ch.unil.genescore.vegas.test;

import static org.junit.Assert.assertEquals;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.Snp;

import org.apache.commons.math3.distribution.NormalDistribution;

public class SnpTest {
	
	
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test test statistic */
	@Test
	public void testComputeChiSquaredStatistics() {
		
		// Computed chi2 statistics for 1e-1, 1e-2, ... using R:
		//    > x <- 10^(-(1:32))
		//    > qchisq(x, 1, lower.tail=FALSE)
		double[] chi2Stat = { 2.705543, 6.634897, 10.827566, 15.136705, 19.511421, 23.928127, 28.373987, 32.841253, 37.324893, 41.821456, 
				46.328476, 50.844128, 55.367025, 59.896088, 64.430464, 68.969461, 73.512517, 78.059165, 82.609014, 87.161733,
				91.717041, 96.274696, 100.834486, 105.396229, 109.959765, 114.524952, 119.091664, 123.659790, 128.229229, 132.799893,
				137.371700, 141.944577 }; 
		
		//Main.println("Comparing chi2 stat computed using:");
		//Main.println("p-value\tR\tMath Commons");

		for (int i=1; i<14; i++) {
			Snp snp = new Snp("snp1", Math.pow(10, -i));
			snp.computeChiSquaredStatistics();
			//Main.println("1e-" + i + "\t" + chi2Stat[i-1] + "\t" + snp.getChi2Stat());
			assertEquals(chi2Stat[i-1], snp.getChi2Stat(), 1e-2);
		}
	}

}
	

