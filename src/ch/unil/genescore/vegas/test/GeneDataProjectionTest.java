package ch.unil.genescore.vegas.test;

import java.util.ArrayList;
import java.util.LinkedList;

import no.uib.cipr.matrix.DenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GeneDataProjection;
import ch.unil.genescore.vegas.LdMat;
import ch.unil.genescore.vegas.ProjectionVegas;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpWeightPairs;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class GeneDataProjectionTest {
	
	
	class stubLdMat extends LdMat {
	
		private double[] rowWeights_=null;
		public stubLdMat(SnpWeightPairs rowSnpWeightPairs, ArrayList<Snp> columnSnpList) {
			super(rowSnpWeightPairs, columnSnpList);
		}			
		public void setRowWeights(double[] rowWeights){
			rowWeights_=rowWeights;
		}
		public double[] getRowWeights(){
			return rowWeights_;
		}
	}
	
	class  TestableGeneDataProject extends GeneDataProjection {

		public DenseMatrix testMatToDecompose_=null;
	
		public TestableGeneDataProject(ArrayList<Snp> snpList,
			LinkedList<Snp> currentSnpsWithGenotypes,
			SnpWeightPairs SnpsAndWeights) {
			super(snpList, currentSnpsWithGenotypes, SnpsAndWeights);
		}
		public void processData(){
			calculateSnpwiseVar();			
		}	
		protected DenseMatrix getMatToDecompose(){
			return testMatToDecompose_;
		}	
		protected void setMatToDecompose(DenseMatrix testMatToDecompose){
			testMatToDecompose_=testMatToDecompose;
		}	
		protected void setWeightedCrossCorrMat(LdMat weightedCrossCorrMat){
			weightedCrossCorrMat_=weightedCrossCorrMat;			
		}
		protected void setSnpWiseVar(double[] snpWiseVar){
			snpWiseVar_=snpWiseVar;
		}
	}
		
		
		
	
	
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
	public void testCalculateSnpwiseVar() {
		
		DenseMatrix myMat = new DenseMatrix(3,3);
		myMat.set(0, 0, 1);myMat.set(1, 1, 1.5);myMat.set(2, 2, 0.9);
		
		TestableGeneDataProject testable = new TestableGeneDataProject(null,null,null);
		testable.setMatToDecompose(myMat);
		
		double[] rowWeights={2,1,1};
		stubLdMat myLdMat = new stubLdMat(null,null);
		myLdMat.setRowWeights(rowWeights);
		testable.setWeightedCrossCorrMat(myLdMat);
		testable.processData();
		double[] res = testable.getSnpWiseVar();
		double[] target = {0.25,1.5,0.9};
		assertTrue(res[0]==target[0] && res[1]==target[1] && res[2]==target[2]);
		
		
	}
	@Test
	public void testGetCorr() {
		
		
		
		DenseMatrix myMat = new DenseMatrix(3,3);
		myMat.set(0, 0, 1);myMat.set(1, 1, 1.5);myMat.set(2, 2, 0.9);
		myMat.set(0, 0, 1);myMat.set(1, 1, 1.5);myMat.set(2, 2, 0.9);
		TestableGeneDataProject testable = new TestableGeneDataProject(null,null,null);
		testable.setMatToDecompose(myMat);
		
		
		double[] mySnpWiseVars={0.5,0.3,0.9};
		testable.setSnpWiseVar(mySnpWiseVars);
		testable.setSnpWiseVarCutoff(0.3);//stuff just on border gets cut
		DenseMatrix res = testable.getCorr();
		assertTrue( res.numColumns()==2 && res.numRows()==2);
		
		testable.setSnpWiseVarCutoff(1);//everything empty
		DenseMatrix res2 = testable.getCorr();
		assertTrue( res2.numColumns()==0 && res2.numRows()==0);
		
	}
	@Test
	public void testCaluculateWeights() {
		
		
		double[] rowWeights= {10,1,2, 3};
		stubLdMat weightedStub = new stubLdMat(null,null);
		weightedStub.setRowWeights(rowWeights);
		
		TestableGeneDataProject testable = new TestableGeneDataProject(null,null,null);
		testable.setWeightedCrossCorrMat(weightedStub);
		double[] snpWiseVar ={1,1,0.5, 1};
		testable.setSnpWiseVar(snpWiseVar);
		testable.setSnpWiseVarCutoff(0.99);
		testable.calculateWeights();
		
		double[] myWeights = testable.returnWeights(false);
		double[] weightsTarget = {10, 1, 3};
		for (int i = 0; i< myWeights.length; i++)
			assertTrue(myWeights[i]==weightsTarget[i]);
		
		
	}
}