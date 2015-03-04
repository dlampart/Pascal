package ch.unil.genescore.eqtlOverlap.test;

import java.util.ArrayList;

import static org.junit.Assert.*;
import no.uib.cipr.matrix.DenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.eqtlOverlap.EqtlEvaluator;
import ch.unil.genescore.eqtlOverlap.GeneDataWithEqtl;
import ch.unil.genescore.main.Settings;


public class EqtlEvaluatorTest{
	class TestableEqtlEvaluator extends EqtlEvaluator {

		public TestableEqtlEvaluator(GeneDataWithEqtl geneData) {
			super(geneData);
		}
		public TestableEqtlEvaluator(DenseMatrix LdMat, double[] eqtlZscores, double[] gwasZscores) {
			LdMat_=LdMat;
			eqtlZscores_=eqtlZscores;
			gwasZscores_=gwasZscores;
		}

	}
	// SETUP
	
		@BeforeClass
		public static void testSetup() {
		//	Settings.loadSettings();
		}

		@AfterClass
		public static void testCleanup() { }
		
		
		// ============================================================================
		// TESTS
		
		@Test
		public void computeScore_test(){
						
			
			
			DenseMatrix ld= new DenseMatrix(3,3);			
			ld.set(0,0,1);ld.set(1,1,1);ld.set(2,2,1);
			ld.set(0,1,0.9);ld.set(1,2,0.9);
			ld.set(1,0,0.9);ld.set(2,1,0.9);
			ld.set(0,2,0.81);
			ld.set(2,0,0.81);
			double[] beta = { 1,0.8,1};
			double[] x = {2.394,2.4466,3.582};
			
			TestableEqtlEvaluator testObj = new TestableEqtlEvaluator(ld, beta,x);
			testObj.computeScore();
			double[] score = testObj.getScore();
			assertEquals(score[0], 0.000806, 0.00001);								
		}
}