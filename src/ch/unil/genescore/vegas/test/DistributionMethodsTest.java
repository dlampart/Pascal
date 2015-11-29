/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.vegas.test;

import static org.junit.Assert.assertEquals;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.vegas.DistributionMethods;


	
public class DistributionMethodsTest {

	//private static ChiSquaredDistribution chiSquared1df_ = new ChiSquaredDistribution(1);
	//private static NormalDistribution normalDist_ = new NormalDistribution();
	
	// ============================================================================
		// SETUP
		
	@BeforeClass
	public static void testSetup() {
	}

	@AfterClass
	public static void testCleanup() { }
	
	@Test
	public void testChiSquared1dfCumulativeProbabilityUpperTail(){
		
		double val=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(49);
		System.out.println(val);
		assertEquals(2.559619e-12,val, 0.01e-12);
		val=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(60);
		System.out.println(val);
		assertEquals(9.436896e-15,val, 0.3e-15);
		val=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(70);
		System.out.println(val);
		assertEquals(5.930446e-17,val, 0.3e-17);
		val=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(1373.873);
		System.out.println(val);
		assertEquals(1e-300,val, 0.01e-300);
	}
	
	@Test
	public void testChiSquared1dfInverseCumulativeProbabilityUpperTail(){
		
		double val=DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(1e-10);
		System.out.println(val);
		assertEquals(41.82146,val, 0.1);
		val=DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(2e-14);
		System.out.println(val);
		assertEquals(58.53368,val, 0.1);
		val=DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(1e-14);
		System.out.println(val);
		assertEquals(59.8976,val, 0.1);
		val=DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(1e-18);
		assertEquals(78.05916,val, 0.1);
		val=DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(1e-200);
		System.out.println(val);
		assertEquals(913.7627,val, 0.1);
		val=DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(1e-200);
		System.out.println(val);
		assertEquals(913.7627,val, 0.1);
		val=DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(1e-300);
		System.out.println(val);
		assertEquals(1373.873,val, 0.1);
	}
	
	@Test
	public void testNormalCumulativeProbability(){
		
		double val=DistributionMethods.normalCumulativeProbability(3);
		assertEquals( 0.9986501,val, 0.000001);
		val=DistributionMethods.normalCumulativeProbability(-8);
		assertEquals(6.220961e-16,val,0.2e-16);
		val=DistributionMethods.normalCumulativeProbability(-8.1);
		assertEquals(2.747959e-16,val,0.2e-16);
		val=DistributionMethods.normalCumulativeProbability(-20);
		assertEquals(2.753624e-89,val,0.2e-89);
	}
	
	@Test
	public void TestNormalInverseCumulativeProbability(){
		
		double val=DistributionMethods.normalInverseCumulativeProbability(1e-5);
		assertEquals( -4.2648910,val, 0.0001);		
		val=DistributionMethods.normalInverseCumulativeProbability(0.9999);
		assertEquals( 3.719016,val, 0.0001);
		val=DistributionMethods.normalInverseCumulativeProbability(1e-14);
		assertEquals(-7.650628,val, 0.01);
		val=DistributionMethods.normalInverseCumulativeProbability(0.5e-14);
		assertEquals(-7.739256,val, 0.01);
		val=DistributionMethods.normalInverseCumulativeProbability(0.5e-18);
		assertEquals(-8.83511,val, 0.01);
		val=DistributionMethods.normalInverseCumulativeProbability(0.5e-300);
		assertEquals(-37.06579,val, 0.01);
		
	}

}
