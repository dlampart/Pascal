package ch.unil.genescore.vegas.test;




import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Set;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpWeightMap;
import ch.unil.genescore.vegas.SnpWeightPairs;

/**
 * Unit tests for GenewiseSnpWeightsTest
 */
public class GenewiseSnpWeightsTest {

	
	static double delta;
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
		delta=1E-14;
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test PathwayMain.run() */
	@Test
	public void testRun() {
		Settings.geneWiseSnpWeightsFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/Chr22CorrelatedStates.txt";
		Settings.useMafCutoffForProjection_= 0.1;
		
		GenewiseSnpWeights testInstance = new GenewiseSnpWeights();
		testInstance.loadSnpWeightsForEachGeneFromFile(null);
		Set<String> allIds =  testInstance.returnAllSnpIds();
		//Snps on chr22:
		assertTrue(allIds.contains("rs143740623"));
		assertTrue(allIds.contains("rs55718707"));
		assertTrue(allIds.contains("rs117058002"));
		assertTrue(allIds.contains("rs147432878"));
		assertTrue(allIds.contains("esv2656764"));		
		//Snps not on chr22:
		assertFalse(allIds.contains("rs186126309"));
		assertFalse(allIds.contains("rs183949376"));
		// Test overlap function
		String geneId="ENSG00000184319";
		// setup snplist
		LinkedList<Snp> mySnpList = new LinkedList<Snp>();
		Snp snpIn1 = new Snp("rs1053764");
		byte[] fakeGenotype = {0,0,1,1,1};
		snpIn1.setGenotypes(fakeGenotype);
		 snpIn1.computeAlleleStats();
		Snp snpIn2 = new Snp("rs200182983");
		snpIn2.setGenotypes(fakeGenotype);
		 snpIn2.computeAlleleStats();		
		Snp snpIn3 = new Snp("rs143893967");
		snpIn3.setGenotypes(fakeGenotype);
		 snpIn3.computeAlleleStats();				
		Snp snpNotIn1 = new Snp("rs201651988");
		snpNotIn1.setGenotypes(fakeGenotype);
		 snpNotIn1.computeAlleleStats();				
		Snp snpNotIn2 = new Snp("rs140138610");
		snpNotIn2.setGenotypes(fakeGenotype);
		 snpNotIn2.computeAlleleStats();						
		mySnpList.add(snpIn1);
		mySnpList.add(snpNotIn1);
		mySnpList.add(snpIn2);
		mySnpList.add(snpNotIn2);
		mySnpList.add(snpIn3);
		//
		SnpWeightPairs myWeightPairs = testInstance.getOverlappedSnpsWithWeights(geneId, mySnpList);
				
		double[] weights = myWeightPairs.getWeightsNotList();
		for (double weight : weights)
			assertEquals(weight,1, 1e-10);
		ArrayList<Snp> snps = myWeightPairs.getSnps();
		assertEquals(snps.size(), 3, 1e-10);
		assertEquals(snps.get(0).id_,"rs1053764");
		assertEquals(snps.get(1).id_,"rs200182983");
		assertEquals(snps.get(2).id_,"rs143893967");
					
	}
	@Test
	public void testUpdateWeights(){
		String[] inputString1 = {"id1;3", "id2;1", "id3;4", "id4;4"};
		String[] inputString2 = {"id2;1", "id1;2", "id5;20", "id4;2"};
		String[] inputString3 = {"id2;1", "id1;2", "id5;20", "id4;2"};
		String[] inputString4 = {"id2;1", "id1;2", "id5;20", "id4;2"};
		String[] inputString5 = {"id2;1", "id1;2", "id5;20", "id4;2"};
		SnpWeightMap myMapMainGene1 = new SnpWeightMap();
		SnpWeightMap myMapAuxGene1 = new SnpWeightMap();
		SnpWeightMap myMapMainGene2 = new SnpWeightMap();
		SnpWeightMap myMapAuxGene2 = new SnpWeightMap();
		SnpWeightMap myMapAuxGene3 = new SnpWeightMap();
		
		myMapMainGene1.readMapFromSplitLine(inputString1);
		myMapAuxGene1.readMapFromSplitLine(inputString2);
		myMapMainGene2.readMapFromSplitLine(inputString3);
		myMapAuxGene2.readMapFromSplitLine(inputString4);		
		myMapAuxGene3.readMapFromSplitLine(inputString5);
		
		String geneName1 = "g1";
		String geneName2 = "g2";
		String geneName3 = "g3";
		
		HashMap<String, SnpWeightMap> main = new HashMap<String, SnpWeightMap>();
		main.put(geneName1, myMapMainGene1);
		main.put(geneName2, myMapMainGene2);
		
		HashMap<String, SnpWeightMap> aux = new HashMap<String, SnpWeightMap>();
		aux.put(geneName1, myMapAuxGene1);
		aux.put(geneName2, myMapAuxGene2);
		aux.put(geneName3, myMapAuxGene3);
		
		GenewiseSnpWeights mainWeights = new  GenewiseSnpWeights(main);
		GenewiseSnpWeights auxWeights = new  GenewiseSnpWeights(aux);
				
		mainWeights.updateWeights(auxWeights);		
		
		HashMap<String, SnpWeightMap> myHash = mainWeights.getHash();
		assertTrue(myHash.get(geneName1).getWeight("id1")==6);
		assertTrue(myHash.get(geneName1).getWeight("id4")==8);
		assertTrue(myHash.get(geneName2).getWeight("id1")==4);
		assertTrue(!myHash.containsKey(geneName3));		
	}
	
}