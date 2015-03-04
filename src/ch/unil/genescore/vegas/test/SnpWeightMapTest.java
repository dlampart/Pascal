package ch.unil.genescore.vegas.test;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Genome;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpWeightMap;
public class SnpWeightMapTest {
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	@Test
	public void testReadMapFromSplitLine(){ 
		
		Settings.useMafCutoffForProjection_= 0.3;
		
		SnpWeightMap myMap = new SnpWeightMap();
		String[] inputString = {"id1;3", "id2;1", "id3;4", "id4;4"};
		ArrayList<Snp> snpList = new ArrayList<Snp>();
		
		Snp snp1 = new Snp("id2");		
		
		snpList.add(snp1);
		Snp.setGenotypeIsPhased(true);
		snp1.setGenotypeIsPhased2(true);
		Snp.setGenotypeLength(5);
		byte[] fakeGenotype = {0,0,1,1,1};
		snp1.setGenotypes(fakeGenotype);
		snp1.computeAlleleStats();		 
		 
		Snp snp2 = new Snp("id5");
			
		snpList.add(snp2);
		snp2.setGenotypeIsPhased2(true);
		byte[] fakeGenotype2 = {0,0,1,1,1};
		snp2.setGenotypes(fakeGenotype2);
		snp2.computeAlleleStats();
		
		Snp snp3 = new Snp("id3");
		
		snpList.add(snp3);
		snp3.setGenotypeIsPhased2(true);
		byte[] fakeGenotype3 = {0,0,0,0,1};
		snp3.setGenotypes(fakeGenotype3);
		snp3.computeAlleleStats();
		
		Snp snp4 = new Snp("id4");
		snp3.setGenotypeIsPhased2(true);
		snpList.add(snp4);
		byte[] fakeGenotype4 = {0,0,0,1,1,1};
		snp4.setGenotypes(fakeGenotype4);
		snp4.computeAlleleStats();
		
		
		myMap.readMapFromSplitLine(inputString);
		myMap.findOverlappedSnps(snpList);
		double[] weightResults=myMap.getOverlappedWeightsNotList();
		assertEquals(weightResults[0],1, 1e-10);
		assertEquals(weightResults[1],4, 1e-10);
		
		}
	@Test
	public void testUpdateWeights(){
		String[] inputString1 = {"id1;3", "id2;1", "id3;4", "id4;4"};
		String[] inputString2 = {"id2;1", "id1;2", "id5;20", "id4;2"};
		SnpWeightMap myMapMain = new SnpWeightMap();
		SnpWeightMap myMapAux = new SnpWeightMap();
		myMapMain.readMapFromSplitLine(inputString1);
		myMapAux.readMapFromSplitLine(inputString2);
		myMapMain.updateWeights(myMapAux);		
		assertTrue(myMapMain.getWeight("id1")==6);
		assertTrue(myMapMain.getWeight("id2")==1);
		assertTrue(myMapMain.getWeight("id3")==4);
		assertTrue(myMapMain.getWeight("id4")==8);
		assertTrue(!myMapMain.checkExistence("id5"));
	}
	
}


