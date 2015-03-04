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

import static org.junit.Assert.*;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.AnalyticVegas;
import ch.unil.genescore.vegas.MaxVegas;
import ch.unil.genescore.vegas.MaxVegasWithoutPruning;
import ch.unil.genescore.vegas.ProjectionVegas;
import ch.unil.genescore.vegas.Snp;

/**
 * Unit tests for MaxVegas
 */
public class MaxVegasTest {
	

	
	class TestableMaxVegas extends MaxVegas{

		public TestableMaxVegas(ArrayList<Double> snpScores, DenseMatrix ld, double pruningCutoff) {
			super(snpScores, ld, pruningCutoff);			
		}
		public TestableMaxVegas(ArrayList<Double> snpScores, DenseMatrix ld, double[] weights) {
			super(snpScores, ld, weights);			
		}
		public void setCorrel(UpperSymmDenseMatrix mat){correl_= mat;}
		public void setUpper(double[] upper){upper_= upper;}
		public void setVariances(double[] variances){variances_=variances;}
		public int[] getPruningIndices(){return pruningIndices_;}
		public double[] getUpper(){return upper_;}
		public UpperSymmDenseMatrix getCorrel(){return correl_;}
		
		
		
	} 
	
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
	
	@Test
	public void deWeight_test(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(4,4);
		testMat.set(0,0,1);testMat.set(1,1,.81);testMat.set(2,2,1);testMat.set(3,3,1);
		testMat.set(0,1,.81);testMat.set(1,0,.81);
		testMat.set(1,2,.81);testMat.set(2,1,.81);
		testMat.set(0,2,0.81);testMat.set(2,0,0.81);		
		/** 1.81.81  0
		 **.81.81.81 0
		 **.81.81 1  0
		 ** 0  0  0  1 		 
		 * */		
		geneScores.add(1.2);geneScores.add((double) 1);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,0.9,1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.deweight();
		UpperSymmDenseMatrix myMat = testable.getCorrel();
		double[] myUpper=testable.getUpper();
		double[] myUpperTarget={1.2,1.111,1,1};	
		for (int  i=0;i<3;i++)
			assertEquals(myUpper[i],myUpperTarget[i], 1e-3);	
		assertEquals(myMat.get(0,1),0.9, 1e-3);
		assertEquals(myMat.get(0,2),0.81,1e-3);
		assertEquals(myMat.get(1,2),0.9,1e-3);					
	}
	@Test
	public void stretchToVarOne_test(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(4,4);
		testMat.set(0,0,1);testMat.set(1,1,.81);testMat.set(2,2,1);testMat.set(3,3,1);
		testMat.set(0,1,.81);testMat.set(1,0,.81);
		testMat.set(1,2,.81);testMat.set(2,1,.81);
		testMat.set(0,2,0.81);testMat.set(2,0,0.81);		
		/** 1.81.81  0
		 **.81.81.81 0
		 **.81.81 1  0
		 ** 0  0  0  1 		 
		 * */		
		geneScores.add(1.2);geneScores.add((double) 0.9);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,1,1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.deweight();
		testable.stretchToVarOne();
		UpperSymmDenseMatrix myMat = testable.getCorrel();
		double[] myUpper=testable.getUpper();
		double[] myUpperTarget={1.2,1,1,1};	
		for (int  i=0;i<3;i++)
			assertEquals(myUpper[i],myUpperTarget[i], 1e-3);	
		assertEquals(myMat.get(0,1),0.9, 1e-3);
		assertEquals(myMat.get(0,2),0.81,1e-3);
		assertEquals(myMat.get(1,2),0.9,1e-3);			

	}
	@Test
	public void calculateWeightedMaximum_test(){
		
		Settings.withZScore_=true;		
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		geneScores.add(1.2);geneScores.add((double) 0.9);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,1,1,1};
		double[] upper={1.2,0.9,1,1};
		double[] variances={1,1,1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,new DenseMatrix(1,1),weights);
		testable.setUpper(upper);
		testable.setVariances(variances);		
		testable.calculateWeightedMaximum();
		double[] myUpper = testable.getUpper();
		double[] myUpperTarget = {1.2,1.2,1.2,1.2};
		for (int  i=0;i<3;i++)
			assertEquals(myUpper[i],myUpperTarget[i], 1e-3);	
	}
	
	@Test
	public void calculateWeightedMaximum_withDifferentVariances_test(){
		
		Settings.withZScore_=true;		
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		geneScores.add(1.2);geneScores.add((double) 0.9);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,1,1,1};
		double[] upper={1.0,1,1,1};
		double[] variances={0.49,1,0.49,0.49};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,new DenseMatrix(1,1),weights);
		testable.setUpper(upper);
		testable.setVariances(variances);		
		testable.calculateWeightedMaximum();
		double[] myUpper = testable.getUpper();
		double[] myUpperTarget = {1.428,1,1.428,1.428};
		for (int  i=0;i<3;i++)
			assertEquals(myUpper[i],myUpperTarget[i], 1e-2);	
	}
	
	@Test
	public void prune_withoutWeights_test(){
		
		UpperSymmDenseMatrix testMat= new UpperSymmDenseMatrix(4);		
		testMat.set(0,0,1);testMat.set(1,1,1);testMat.set(2,2,1);testMat.set(3,3,1);
		testMat.set(0,1,.99);
		//testMat.set(1,0,.99);
		//testMat.set(3,2,.89);
		testMat.set(2,3,.89);
		
		/** 1 .99 0  0
		 **.99 1  0  0
		 ** 0  0  1 .89
		 ** 0  0 .89 1 		 
		 * */		
		Settings.withZScore_=true;		
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		geneScores.add(1.2);geneScores.add((double) 0.9);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,1,1,1};
		double[] upper={1,1,2,2};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,new DenseMatrix(1,1),weights);
		testable.setPruningCutoff(0.9);
		testable.setCorrel(testMat);
		testable.setUpper(upper);
		testable.prune();
		double[] myUpper = testable.getUpper();
		UpperSymmDenseMatrix myCorrel  = testable.getCorrel();
		double[] myUpperTarget = {1,2,2};
		for (int  i=0;i<3;i++)
			assertEquals(myUpper[i],myUpperTarget[i], 1e-3);	
		assertEquals(myCorrel.get(0, 1), 0, 1e-4);
		assertEquals(myCorrel.get(0, 2), 0, 1e-4);
		assertEquals(myCorrel.get(1, 2), 0.89, 1e-4);
	}
	
	@Test
	public void prune_withWeights_test(){
		
		UpperSymmDenseMatrix testMat= new UpperSymmDenseMatrix(4);		
		testMat.set(0,0,1);testMat.set(1,1,1);testMat.set(2,2,1);testMat.set(3,3,1);
		testMat.set(0,1,.99);
		//testMat.set(1,0,.99);
		//testMat.set(3,2,.89);
		testMat.set(2,3,.89);
		
		/** 1 .99 0  0
		 **.99 1  0  0
		 ** 0  0  1 .89
		 ** 0  0 .89 1 		 
		 * */		
		Settings.withZScore_=true;		
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		geneScores.add(1.2);geneScores.add((double) 0.9);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,2,1,1};
		double[] upper={1,0.86,1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,new DenseMatrix(1,1),weights);
		testable.setPruningCutoff(0.9);
		testable.setCorrel(testMat);
		testable.setUpper(upper);
		testable.prune();
		double[] myUpper = testable.getUpper();
		UpperSymmDenseMatrix myCorrel  = testable.getCorrel();
		double[] myUpperTarget = {0.86,1,1};
		for (int  i=0;i<3;i++)
			assertEquals(myUpper[i],myUpperTarget[i], 1e-3);	
		assertEquals(myCorrel.get(0, 1), 0, 1e-4);
		assertEquals(myCorrel.get(0, 2), 0, 1e-4);
		assertEquals(myCorrel.get(1, 2), 0.89, 1e-4);
	}
	
	@Test
	public void compareToMaxVegasWithoutPruningWithoutWeighting(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(4,4);
		testMat.set(0,0,1);testMat.set(1,1,.81);testMat.set(2,2,1);testMat.set(3,3,1);
		testMat.set(0,1,.81);testMat.set(1,0,.81);
		testMat.set(1,2,.81);testMat.set(2,1,.81);
		testMat.set(0,2,0.81);testMat.set(2,0,0.81);		
		/** 1.81.81  0
		 **.81.81.81 0
		 **.81.81 1  0
		 ** 0  0  0  1 		 
		 * */		
		geneScores.add(1.2);geneScores.add((double) 0.9);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,1,1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.setPruningCutoff(1.1);
		testable.computeScore();
		double[] firstScore = testable.getScore();
		MaxVegasWithoutPruning testable2 = new MaxVegasWithoutPruning(geneScores,testMat);
		testable2.computeScore();
		double[] secScore = testable2.getScore();
		System.out.println("adsf");
		assertEquals(firstScore[0],secScore[0],0.001);
	}
	@Test
	public void compareToMaxVegasWithoutPruningWithWeighting(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(4,4);
		testMat.set(0,0,1);testMat.set(1,1,.81);testMat.set(2,2,1);testMat.set(3,3,4);
		testMat.set(0,1,.81);testMat.set(1,0,.81);
		testMat.set(1,2,.81);testMat.set(2,1,.81);
		testMat.set(0,2,0.81);testMat.set(2,0,0.81);		
		/** 1.81.81  0
		 **.81.81.81 0
		 **.81.81 1  0
		 ** 0  0  0  1 		 
		 * */		
		geneScores.add(1.2);geneScores.add((double) 0.9);geneScores.add(1.0);geneScores.add(2.0);		
		double[] weights={1,1,1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.setPruningCutoff(1.1);
		testable.computeScore();
		double[] firstScore = testable.getScore();
		MaxVegasWithoutPruning testable2 = new MaxVegasWithoutPruning(geneScores,testMat);
		testable2.computeScore();
		double[] secScore = testable2.getScore();
		assertEquals(firstScore[0],secScore[0],0.001);
		
	}
	
	@Test
	public void compareToMaxVegasWithoutPruningWithWeighting2(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(4);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(4,4);
		testMat.set(0,0,1);testMat.set(1,1,.81);testMat.set(2,2,1);testMat.set(3,3,1);
		testMat.set(0,1,.81);testMat.set(1,0,.81);
		testMat.set(1,2,.81);testMat.set(2,1,.81);
		testMat.set(0,2,0.81);testMat.set(2,0,0.81);		
		/** 1.81.81  0
		 **.81.81.81 0
		 **.81.81 1  0
		 ** 0  0  0  1 		 
		 * */		
		geneScores.add(1.2);geneScores.add((double) 1.5);geneScores.add(1.0);geneScores.add(1.0);		
		double[] weights={1,1,1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.setPruningCutoff(1.1);
		testable.computeScore();
		double[] firstScore = testable.getScore();
		MaxVegasWithoutPruning testable2 = new MaxVegasWithoutPruning(geneScores,testMat);
		testable2.computeScore();
		double[] secScore = testable2.getScore();
		
		assertEquals(firstScore[0],secScore[0],0.001);
	}

	@Test
	public void testManual(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(2);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(2,2);
		testMat.set(0,0,1);testMat.set(1,1,.04);
		/** 1 0
		 * 0 0.04
		 * */		
		geneScores.add((double) 1);geneScores.add((double) 0.5);
		double[] weights={1,1};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.setPruningCutoff(1.1);
		testable.computeScore();
		double[] firstScore = testable.getScore();
		MaxVegasWithoutPruning testable2 = new MaxVegasWithoutPruning(geneScores,testMat);
		testable2.computeScore();
		double[] secScore = testable2.getScore();
		
		assertEquals(firstScore[0],secScore[0],0.00001);
		assertEquals(firstScore[0],0.317,0.001);
	}
	
	@Test
	public void testManual2(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(2);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(2,2);
		testMat.set(0,0,1);testMat.set(1,1,1);
		/** 1 0
		 * 0 1
		 * */		
		geneScores.add((double) 1);geneScores.add((double) 1);
		double[] weights={1,10};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.setPruningCutoff(1.1);
		testable.computeScore();
		double[] firstScore = testable.getScore();
		MaxVegasWithoutPruning testable2 = new MaxVegasWithoutPruning(geneScores,testMat);
		testable2.computeScore();
		double[] secScore = testable2.getScore();
		
		assertEquals(firstScore[0],0.5335,0.001);
		assertEquals(secScore[0],0.5335,0.001);
	}
	
	@Test
	public void testManual3(){
		
		Settings.withZScore_=true;
		ArrayList<Double>  geneScores = new ArrayList<Double>(2);
		//TestableMaxVegas testable = new TestableMaxVegas(new ArrayList<Double>(),new DenseMatrix(4, 4));
		DenseMatrix testMat= new DenseMatrix(2,2);
		testMat.set(0,0,1);testMat.set(1,1,36);
		/** 1 0
		 * 0 1
		 * */		
		geneScores.add((double) 1);geneScores.add((double) 6);
		double[] weights={1,6};
		TestableMaxVegas testable = new TestableMaxVegas(geneScores,testMat,weights);
		testable.setPruningCutoff(1.1);
		testable.computeScore();
		double[] firstScore = testable.getScore();
		MaxVegasWithoutPruning testable2 = new MaxVegasWithoutPruning(geneScores,testMat);
		testable2.computeScore();
		double[] secScore = testable2.getScore();
		
		assertEquals(firstScore[0],secScore[0],0.00001);
		assertEquals(firstScore[0],0.317,0.001);
	}
}

