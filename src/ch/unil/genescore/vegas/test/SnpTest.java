/*
Copyright (c) 2013 Daniel Marbach

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper available at:
http://networkinference.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
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
	

