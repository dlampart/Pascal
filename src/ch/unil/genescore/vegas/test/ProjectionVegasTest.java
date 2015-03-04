package ch.unil.genescore.vegas.test;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.GeneDataProjection;
import ch.unil.genescore.vegas.ProjectionVegas;
import ch.unil.genescore.vegas.Snp;


/**
 * Unit tests for ProjectionVegas.
 * 
 */
public class ProjectionVegasTest {
	
	
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
	public void testConstructRegularizedInverseCovariance1(){
		Settings.fractionToBeExplained_=1;
		Settings.conditionFraction_=0.01;
		Settings.withZScore_ = true;
		
		Snp fakeSnp1 = new Snp("fakeId1", 0, 0.3);
		Snp fakeSnp2 = new Snp("fakeId2", 0, 0.8);
		Snp fakeSnp3 = new Snp("fakeId3", 0, 1.4);
		ArrayList<Snp> geneSnps = new ArrayList<Snp>();
		geneSnps.add(fakeSnp1);geneSnps.add(fakeSnp2);geneSnps.add(fakeSnp3);
		DenseMatrix ld= new DenseMatrix(3,3);
		DenseMatrix crossLd= new DenseMatrix(2,3);
		//DenseMatrix crossLd= new DenseMatrix(3,2);
		
		
		//ld and crossLd calculated as follows:
		// make 0.9-toeplitz mat of size 5.
		// 3,4,5 for ld-mat
		// 1,2 for crossLd-mat
		ld.set(0,0,1);ld.set(1,1,1);ld.set(2,2,1);
		ld.set(0,1,0.9);ld.set(1,2,0.9);
		ld.set(1,0,0.9);ld.set(2,1,0.9);
		ld.set(0,2,0.81);
		ld.set(2,0,0.81);
				
		crossLd.set(0,0,0.81);crossLd.set(1,0,0.9);
		crossLd.set(0,1,0.729);crossLd.set(1,1,0.81);
		crossLd.set(0,2,0.6561);crossLd.set(1,2,0.729);
		
		ArrayList<Double> snpScores = new ArrayList<Double>();
		snpScores.add(0.3);
		snpScores.add(0.8);
		snpScores.add(1.4);		
		ProjectionVegas myProjectionObj= new ProjectionVegas(snpScores, ld,crossLd,0); 
		myProjectionObj.processData();
		DenseMatrix regularizedInvCovTarget=new DenseMatrix(3,3);
		regularizedInvCovTarget.set(0,0,4.8143905);
		regularizedInvCovTarget.set(0,1,-4.124656);
		regularizedInvCovTarget.set(0,2,-0.1856095);
		regularizedInvCovTarget.set(1,0,-4.1246563);
		regularizedInvCovTarget.set(1,1,8.340972);
		regularizedInvCovTarget.set(1,2,-4.1246563);
		regularizedInvCovTarget.set(2,0,-0.1856095);
		regularizedInvCovTarget.set(2,1,-4.1246563);
		regularizedInvCovTarget.set(2,2,4.8143905);
		
		DenseMatrix regularizedInvCovResult=myProjectionObj.getRegularizedInvCovariance();		
		double[] dataTarget = regularizedInvCovTarget.getData();
		double[] dataResult = regularizedInvCovResult.getData();		
		for (int i=0; i < dataTarget.length;i++)
			assertEquals(dataTarget[i],dataResult[i], 1E-4);
		
		
		
		DenseMatrix matToDecomposeTarget = new DenseMatrix(2,2);
		matToDecomposeTarget.set(0,0,0.6498549);		
		matToDecomposeTarget.set(0,1,0.722061);
		matToDecomposeTarget.set(1,0,0.722061);
		matToDecomposeTarget.set(1,1,0.802290);
		
		DenseMatrix matToDecomposeResult=myProjectionObj.getMatToDecomposeMTJ();		
		dataTarget = matToDecomposeTarget.getData();
		dataResult = matToDecomposeResult.getData();		
		for (int i=0; i < dataTarget.length;i++)
			assertEquals(dataTarget[i],dataResult[i], 1E-4);
		
		DenseVector transformedZscoresTarget=new DenseVector(2);		
		transformedZscoresTarget=new DenseVector(2);
		transformedZscoresTarget.set(0, 0.2601336);
		transformedZscoresTarget.set(1, 0.2890374);
		
		DenseVector transformedZscoresResult=myProjectionObj.getTransformedZscores();
		dataTarget = transformedZscoresTarget.getData();
		dataResult = transformedZscoresTarget.getData();
		for (int i=0; i < dataTarget.length;i++)
			assertEquals(dataTarget[i],dataResult[i], 1E-4);				
	}
		
	/** Test ProjectionVegas other settings*/
	@Test
	public void testConstructRegularizedInverseCovariance2(){
		Settings.fractionToBeExplained_=0.6;
		Settings.conditionFraction_=0;
		Settings.withZScore_ = true;
		
		Snp fakeSnp1 = new Snp("fakeId1", 0, 0.3);
		Snp fakeSnp2 = new Snp("fakeId2", 0, 0.8);
		Snp fakeSnp3 = new Snp("fakeId3", 0, 1.4);
		ArrayList<Snp> geneSnps = new ArrayList<Snp>();
		geneSnps.add(fakeSnp1);geneSnps.add(fakeSnp2);geneSnps.add(fakeSnp3);
		DenseMatrix ld= new DenseMatrix(3,3);
		DenseMatrix crossLd= new DenseMatrix(2,3);
		//DenseMatrix crossLd= new DenseMatrix(3,2);
		
		
		//ld and crossLd calculated as follows:
		// make 0.9-toeplitz mat of size 5.
		// 3,4,5 for ld-mat
		// 1,2 for crossLd-mat
		ld.set(0,0,1);ld.set(1,1,1);ld.set(2,2,1);
		ld.set(0,1,0.9);ld.set(1,2,0.9);
		ld.set(1,0,0.9);ld.set(2,1,0.9);
		ld.set(0,2,0.81);
		ld.set(2,0,0.81);
				
		crossLd.set(0,0,0.81);crossLd.set(1,0,0.9);
		crossLd.set(0,1,0.729);crossLd.set(1,1,0.81);
		crossLd.set(0,2,0.6561);crossLd.set(1,2,0.729);
		
		ArrayList<Double> snpScores = new ArrayList<Double>();
		snpScores.add(0.3);
		snpScores.add(0.8);
		snpScores.add(1.4);		
		ProjectionVegas myProjectionObj= new ProjectionVegas(snpScores, ld,crossLd, 0); 
		myProjectionObj.processData();
//		DenseMatrix regularizedInvCovTarget=new DenseMatrix(3,3);
//		regularizedInvCovTarget.set(0,0,4.8143905);
//		regularizedInvCovTarget.set(0,1,-4.124656);
//		regularizedInvCovTarget.set(0,2,-0.1856095);
//		regularizedInvCovTarget.set(1,0,-4.1246563);
//		regularizedInvCovTarget.set(1,1,8.340972);
//		regularizedInvCovTarget.set(1,2,-4.1246563);
//		regularizedInvCovTarget.set(2,0,-0.1856095);
//		regularizedInvCovTarget.set(2,1,-4.1246563);
//		regularizedInvCovTarget.set(2,2,4.8143905);
//		
		DenseMatrix regularizedInvCovResult=myProjectionObj.getRegularizedInvCovariance();		
//		double[] dataTarget = regularizedInvCovTarget.getData();
//		double[] dataResult = regularizedInvCovResult.getData();		
//		for (int i=0; i < dataTarget.length;i++)
//			assertEquals(dataTarget[i],dataResult[i], 1E-4);
//		
//		
		
		DenseMatrix matToDecomposeTarget = new DenseMatrix(2,2);
		matToDecomposeTarget.set(0,0,0.5858472);		
		matToDecomposeTarget.set(0,1,0.6509414);
		matToDecomposeTarget.set(1,0,0.6509414);
		matToDecomposeTarget.set(1,1,0.7232682);
		
		DenseMatrix matToDecomposeResult=myProjectionObj.getMatToDecomposeMTJ();		
		double[] dataTarget = matToDecomposeTarget.getData();
		double[] dataResult = matToDecomposeResult.getData();		
		for (int i=0; i < dataTarget.length;i++)
			assertEquals(dataTarget[i],dataResult[i], 1E-4);
		
		DenseVector transformedZscoresTarget=new DenseVector(2);		
		transformedZscoresTarget=new DenseVector(2);
		transformedZscoresTarget.set(0,0.6669494);
		transformedZscoresTarget.set(1,0.7410549);
		
		DenseVector transformedZscoresResult=myProjectionObj.getTransformedZscores();
		dataTarget = transformedZscoresTarget.getData();
		dataResult = transformedZscoresTarget.getData();
		for (int i=0; i < dataTarget.length;i++)
			assertEquals(dataTarget[i],dataResult[i], 1E-4);				
	}
			
	
}

