package ch.unil.genescore.eqtlOverlap.test;

import java.util.ArrayList;

import static org.junit.Assert.*;
import no.uib.cipr.matrix.DenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.eqtlOverlap.EqtlEvaluator;
import ch.unil.genescore.eqtlOverlap.EqtlGenePair;
import ch.unil.genescore.eqtlOverlap.GeneDataWithEqtl;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.Snp;

public class EqtlGenePairTest {
		
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
	//	Settings.loadSettings();
	}

	@AfterClass
	public static void testCleanup() { 
		
	}
	
	
	// ============================================================================
	// TESTS
	
	@Test
	public void overlapSnps_test(){
		ArrayList<Snp> eqtlSnpList = new ArrayList<Snp>();
		Snp snp = null;
		snp = new Snp("id1",1,1.5);snp.setMajorAllele('a');snp.setMinorAllele('e');eqtlSnpList.add(snp);
		snp = new Snp("id2",2,2.5);snp.setMajorAllele('b');snp.setMinorAllele('f');eqtlSnpList.add(snp);
		snp = new Snp("id3",3,3.5);snp.setMajorAllele('c');snp.setMinorAllele('g');eqtlSnpList.add(snp);
		snp = new Snp("id4",4,4.5);snp.setMajorAllele('d');snp.setMinorAllele('h');eqtlSnpList.add(snp);
		
		ArrayList<Snp> gwasSnpList = new ArrayList<Snp>();				
		snp = new Snp("id3",2,2.6);snp.setMajorAllele('b');snp.setMinorAllele('f');gwasSnpList.add(snp);
		snp = new Snp("id1",1,1.6);snp.setMajorAllele('a');snp.setMinorAllele('e');gwasSnpList.add(snp);
		snp = new Snp("id5",3,3.6);snp.setMajorAllele('c');snp.setMinorAllele('g');gwasSnpList.add(snp);
		snp = new Snp("id0",4,4.6);snp.setMajorAllele('d');snp.setMinorAllele('h');gwasSnpList.add(snp);
		snp = new Snp("id4",4,4.6);snp.setMajorAllele('h');snp.setMinorAllele('d');gwasSnpList.add(snp);		

		EqtlGenePair testObj = new EqtlGenePair("geneName", eqtlSnpList);
		
		testObj.overlapSnps(gwasSnpList);
		double[] scores = testObj.getOverlappedEqtlZscores();
		ArrayList<Snp> outList = testObj.getOverlappedSnpList();
		assertEquals(scores[0],1.5, 0);assertEquals(scores[1],-4.5, 0);
		snp = outList.get(0);
		assertTrue(snp.id_=="id1");assertEquals(snp.getZscore(),1.6, 0);
		snp = outList.get(1);
		assertTrue(snp.id_=="id4");assertEquals(snp.getZscore(),4.6, 0);
	}
	}
