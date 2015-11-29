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

import static org.junit.Assert.*;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import ch.unil.genescore.main.Pascal;
import ch.unil.genescore.vegas.AnalyticVegas;
import ch.unil.genescore.vegas.Snp;


public class AnalyticVegasTest {
	
	/** The Pascal instance (initializes Settings) */
	private static Pascal psc = new Pascal();

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Pascal.set.resetToDefaults();
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test test statistic */
	@Test
	public void test() {
		double delta = 0.00001;
		Snp fakeSnp1 = new Snp("fakeId1", 0, 0.3);
		Snp fakeSnp2 = new Snp("fakeId2", 0, 0.8);
		Snp fakeSnp3 = new Snp("fakeId3", 0, 1.4);
		Pascal.set.withZScore_=true;
		ArrayList<Snp> geneSnps = new ArrayList<Snp>();
		geneSnps.add(fakeSnp1);geneSnps.add(fakeSnp2);geneSnps.add(fakeSnp3);
		DenseMatrix ld= new DenseMatrix(3,3);
		//DenseMatrix crossLd= new DenseMatrix(3,2);
		
		ArrayList<Double> snpScores = new ArrayList<Double>(3);
		snpScores.add(fakeSnp1.getZscore());
		snpScores.add(fakeSnp2.getZscore());
		snpScores.add(fakeSnp3.getZscore());
		//ld and crossLd calculated as follows:
		// make 0.9-toeplitz mat of size 5.
		// 3,4,5 for ld-mat
		// 1,2 for crossLd-mat
		ld.set(0,0,1);ld.set(1,1,1);ld.set(2,2,1);
		ld.set(0,1,0.9);ld.set(1,2,0.9);
		ld.set(1,0,0.9);ld.set(2,1,0.9);
		ld.set(0,2,0.81);
		ld.set(2,0,0.81);
		double[] weights = {2,2,2};
		AnalyticVegas myAnalyticObj= null;
		//UpperSymmDenseMatrix myMatToDecompose = null;
		double[] myEigenvals = null;
		double[] emptyWeights = {1,1,1};
		myAnalyticObj= new AnalyticVegas(snpScores, ld, emptyWeights);
		myAnalyticObj.computeScore();
		myEigenvals = myAnalyticObj.computeLambda();
		assertEquals(myEigenvals[0],2.74067, delta);
		assertEquals(myEigenvals[1],0.1900, delta);
		assertEquals(myEigenvals[2],0.06932, delta);
		
		myAnalyticObj= new AnalyticVegas(snpScores, ld, weights);
		myAnalyticObj.computeScore();
		myEigenvals = myAnalyticObj.computeLambda();			
		assertEquals(myEigenvals[0], 5.481348, delta);
		assertEquals(myEigenvals[1],0.380000, delta);
		assertEquals(myEigenvals[2],0.138652, delta);
			
		double[] weights2 = {1,2,0.5};
		
		myAnalyticObj= new AnalyticVegas(snpScores, ld, weights2);
		myAnalyticObj.computeScore();
		myEigenvals = myAnalyticObj.computeLambda();		
		assertEquals(myEigenvals[0],3.27694674, delta);
		assertEquals(myEigenvals[1],0.1492338, delta);
		assertEquals(myEigenvals[2],0.07381938, delta);			
		
	}
	
	public void computeTestStatisticRealTest() {
		
		DenseMatrix ld= new DenseMatrix(3,3);
		ld.set(0,0,1);ld.set(1,1,1);ld.set(2,2,1);
		ld.set(0,1,0.9);ld.set(1,2,0.9);
		ld.set(1,0,0.9);ld.set(2,1,0.9);
		ld.set(0,2,0.81);
		ld.set(2,0,0.81);
		
		//double delta = 0.00001;
		Snp fakeSnp1 = new Snp("fakeId1", 0.7641772, 0);
		Snp fakeSnp2 = new Snp("fakeId2", 0.4237108, 0);
		Snp fakeSnp3 = new Snp("fakeId3", 0.1615133, 0);
		Pascal.set.withZScore_=false;
		ArrayList<Snp> geneSnps = new ArrayList<Snp>();
		geneSnps.add(fakeSnp1);geneSnps.add(fakeSnp2);geneSnps.add(fakeSnp3);
		//
		ChiSquaredDistribution chiSquareDist1 = new ChiSquaredDistribution(1);			
		ArrayList<Double> chiSquaredVals= new ArrayList<Double>();
		chiSquaredVals.add(chiSquareDist1.inverseCumulativeProbability(1-fakeSnp1.getPval()));
		chiSquaredVals.add(chiSquareDist1.inverseCumulativeProbability(1-fakeSnp2.getPval()));
		chiSquaredVals.add(chiSquareDist1.inverseCumulativeProbability(1-fakeSnp3.getPval()));
		//cumulativeProbability(chiStat);
		double[] emptyWeights = {1,1,1};
		AnalyticVegas myAnalyticObj= new AnalyticVegas(chiSquaredVals, ld, emptyWeights);
		myAnalyticObj.computeScore();
		
		
		myAnalyticObj.computeTestStatisticReal();
		double firstRes = myAnalyticObj.getTestStatisticReal();
		
		
		//Snp fakeSnp1_1 = new Snp("fakeId1", 0, 0.3);
		//Snp fakeSnp2_1 = new Snp("fakeId2", 0, 0.8);
		//Snp fakeSnp3_1 = new Snp("fakeId3", 0, 1.4);
		Pascal.set.withZScore_=false;
		ArrayList<Double> zscores= new ArrayList<Double>();
		zscores.add(0.3);
		zscores.add(0.8);
		zscores.add(1.4);
				
		//geneSnps_1.add(fakeSnp1_1);geneSnps_1.add(fakeSnp2_1);geneSnps_1.add(fakeSnp3_1);
		AnalyticVegas myAnalyticObj_1= new AnalyticVegas(zscores, ld, emptyWeights);
		myAnalyticObj.computeScore();
		myAnalyticObj_1.computeTestStatisticReal();
		double secRes = myAnalyticObj_1.getTestStatisticReal();
		assertEquals(firstRes,secRes, 1E-4);
	}
}