package ch.unil.genescore.vegas.test;




import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.vegas.AllOverlappedElements;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.OverlappedCollectionStream;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpWeightPairs;



/**
 * Unit tests for GenewiseSnpWeightsTest
 */
public class AllOverlappedTest {

	
	static double delta;
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {	
		delta=1E-14;
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test PathwayMain.run() */
//	@Test
	public void test1() {
		String chr1="chr19";
		GenomicElement e1 = new GenomicElement("e1");
		e1.setPosition(chr1, 1100, 1200, true);
		ArrayList<OverlappedGenomicElement> eArray = new ArrayList<OverlappedGenomicElement>(2);		
		eArray.add(new OverlappedGenomicElement(e1));
		Snp s1 = new Snp("s1");
		s1.setPosition(chr1, 1100, 1101, true);		
		ArrayList<OverlappedGenomicElement> sArray = new ArrayList<OverlappedGenomicElement>(6);
		sArray.add(new OverlappedGenomicElement(s1));
		OverlappedCollectionStream eStream = new OverlappedCollectionStream(eArray);
		OverlappedCollectionStream sStream = new OverlappedCollectionStream(sArray);		
		AllOverlappedElements firstObj = new AllOverlappedElements(sStream, eStream);				
		GenewiseSnpWeights myObjView= firstObj.getAsGenewiseSnpWeightDataStructure(0.4, 0);	
	    TreeSet<OverlappedGenomicElement> myLargeList = firstObj.getMyLargeList();
		assertTrue(myLargeList.size()==1);
		OverlappedGenomicElement el = myLargeList.first();
		Set<String> myIds = myObjView.returnAllSnpIds();
		assertTrue(myIds.size()==1);
		for (String myId : myIds)
			assertTrue(myId.equals("s1"));
		 System.out.println("blbu");
		 
		 
	}
	
		
//	AllOverlappedElements firstObj = new AllOverlappedElements(sStream, eStream);
	/** Test PathwayMain.run() */
	@Test
	public void test2() {
		String chr1="chr19";
		Gene g1 = new Gene("g1");
		Gene g2 = new Gene("g2");
		g2.setPosition(chr1, 1000, 3000, true);
		g1.setPosition(chr1, 1000, 4000, true);
		ArrayList<OverlappedGenomicElement> gArray = new ArrayList<OverlappedGenomicElement>(2);
		gArray.add(new OverlappedGenomicElement(g1));
		gArray.add(new OverlappedGenomicElement(g2));
		Collections.sort(gArray);
		GenomicElement e2 = new GenomicElement("e2");
		e2.setPosition(chr1, 1905, 1990, true);	
		ArrayList<OverlappedGenomicElement> eArray = new ArrayList<OverlappedGenomicElement>(1);		
		eArray.add(new OverlappedGenomicElement(e2));
		Snp s1 = new Snp("s1");
		//s1.setPosition(chr1, 1900, 2000, true);
		s1.setPosition(chr1, 1910, 1980, true);
		s1.setMaf(0.5);
		ArrayList<OverlappedGenomicElement> sArray = new ArrayList<OverlappedGenomicElement>(1);
		sArray.add(new OverlappedGenomicElement(s1));
		OverlappedCollectionStream gStream = new OverlappedCollectionStream(gArray);
		OverlappedCollectionStream eStream = new OverlappedCollectionStream(eArray);
		OverlappedCollectionStream sStream = new OverlappedCollectionStream(sArray);
		
		AllOverlappedElements firstObj = new AllOverlappedElements(sStream, eStream);
		firstObj.fillTreeSet();
		 GenewiseSnpWeights myObjView= firstObj.getAsGenewiseSnpWeightDataStructure(0.4, 0);
		 TreeSet<OverlappedGenomicElement> myOverlappedGenomicRegions = firstObj.getMyLargeList();		
		 OverlappedCollectionStream myOverlappedGenomicRegionStream = new OverlappedCollectionStream(myOverlappedGenomicRegions);
		 AllOverlappedElements secObj = new AllOverlappedElements(myOverlappedGenomicRegionStream, gStream);
		 secObj.fillTreeSet();
			GenewiseSnpWeights secObjView= secObj.getAsGenewiseSnpWeightDataStructure(0.4, 1);
			LinkedList<Snp> SnpList = new LinkedList<Snp>();
			 SnpList.add(s1);
			 SnpWeightPairs overlappedList =  secObjView.getOverlappedSnpsWithWeights("g1", SnpList);
			 ArrayList<Snp>  finalList = overlappedList.getSnps();
			 assertTrue(finalList.contains(s1));
			 overlappedList =  secObjView.getOverlappedSnpsWithWeights("g2", SnpList);
			 finalList = overlappedList.getSnps();
			 assertTrue(finalList.contains(s1));
			 	 
	}
	@Test
	public void test3(){
		String chr2="chr22";
		Gene g3 = new Gene("g3");
		g3.setPosition(chr2, 1000, 2000, true);
		ArrayList<OverlappedGenomicElement> gArray = new ArrayList<OverlappedGenomicElement>(3);		
		gArray.add(new OverlappedGenomicElement(g3));
		GenomicElement e3 = new GenomicElement("e3");
		e3.setPosition(chr2, 1000, 1200, true);	
		ArrayList<OverlappedGenomicElement> eArray = new ArrayList<OverlappedGenomicElement>(4);
		eArray.add(new OverlappedGenomicElement(e3));
		Snp s3 = new Snp("s3");
		Snp s5 = new Snp("s5");
		s3.setPosition(chr2, 1000, 1200, true);		
		s5.setPosition(chr2, 1000, 1200, true);		
		s3.setMaf(0.5);
		s5.setMaf(0.5);
		ArrayList<OverlappedGenomicElement> sArray = new ArrayList<OverlappedGenomicElement>(6);
		sArray.add(new OverlappedGenomicElement(s3));
		sArray.add(new OverlappedGenomicElement(s5));
		Collections.sort(sArray);		
		OverlappedCollectionStream gStream = new OverlappedCollectionStream(gArray);
		OverlappedCollectionStream eStream = new OverlappedCollectionStream(eArray);
		OverlappedCollectionStream sStream = new OverlappedCollectionStream(sArray);
		AllOverlappedElements firstObj = new AllOverlappedElements(sStream, eStream);				
		firstObj.fillTreeSet();
		GenewiseSnpWeights myObjView= firstObj.getAsGenewiseSnpWeightDataStructure(0.4, 0);
		TreeSet<OverlappedGenomicElement> myOverlappedGenomicRegions = firstObj.getMyLargeList();		
		OverlappedCollectionStream myOverlappedGenomicRegionStream = new OverlappedCollectionStream(myOverlappedGenomicRegions);		
		AllOverlappedElements secObj = new AllOverlappedElements(myOverlappedGenomicRegionStream, gStream);		
		secObj.fillTreeSet();
		GenewiseSnpWeights secObjView= secObj.getAsGenewiseSnpWeightDataStructure(0.4, 1);	
		 LinkedList<Snp> SnpList = new LinkedList<Snp>();		
		 SnpList.add(s3);		 
		 SnpList.add(s5);	
		SnpWeightPairs overlappedList =  secObjView.getOverlappedSnpsWithWeights("g3", SnpList);
		ArrayList<Snp>  finalList = overlappedList.getSnps();;
		 finalList = overlappedList.getSnps();
		 assertTrue(finalList.contains(s5));		 
		 assertTrue(finalList.contains(s3));
		
			
		 
	}
	
	/** Test PathwayMain.run() */
	@Test
	public void testRun() {
		String chr1="chr19";
		String chr2="chr22";
		Gene g1 = new Gene("g1");
		Gene g2 = new Gene("g2");
		Gene g3 = new Gene("g3");
		g1.setPosition(chr1, 1000, 2000, true);
		g2.setPosition(chr1, 1900, 3000, true);
		g3.setPosition(chr2, 1000, 2000, true);
		
		
		ArrayList<OverlappedGenomicElement> gArray = new ArrayList<OverlappedGenomicElement>(3);
		gArray.add(new OverlappedGenomicElement(g1));
		gArray.add(new OverlappedGenomicElement(g2));
		gArray.add(new OverlappedGenomicElement(g3));
		Collections.sort(gArray);
		
		GenomicElement e1 = new GenomicElement("e1");
		GenomicElement e2 = new GenomicElement("e2");
		GenomicElement e3 = new GenomicElement("e3");
		GenomicElement e4 = new GenomicElement("e4");
		
		e1.setPosition(chr1, 1100, 1200, true);
			e3.setPosition(chr2, 1000, 1200, true);		
			e2.setPosition(chr1, 1900, 2000, true);
		e4.setPosition(chr2, 3000, 3100, true);
		
		ArrayList<OverlappedGenomicElement> eArray = new ArrayList<OverlappedGenomicElement>(4);
		eArray.add(new OverlappedGenomicElement(e1));
		eArray.add(new OverlappedGenomicElement(e2));
		eArray.add(new OverlappedGenomicElement(e3));
		eArray.add(new OverlappedGenomicElement(e4));
		
		Collections.sort(eArray);
		
		Snp s1 = new Snp("s1");
		Snp s2 = new Snp("s2");
		Snp s3 = new Snp("s3");
		Snp s4 = new Snp("s4");
		Snp s5 = new Snp("s5");
		Snp s6 = new Snp("s6");
		
		
		s1.setPosition(chr1, 1100, 1101, true);		
		s2.setPosition(chr1, 1900, 2000, true);
		s3.setPosition(chr2, 1000, 1200, true);		
		s4.setPosition(chr2, 3000, 3100, true);
		s5.setPosition(chr2, 1000, 1200, true);		
		s6.setPosition(chr2, 3000, 3100, true);
		s1.setMaf(0.5);
		s2.setMaf(0.5);
		s3.setMaf(0.5);
		s4.setMaf(0.5);
		s5.setMaf(0.5);
		s6.setMaf(0.5);
		
		
		ArrayList<OverlappedGenomicElement> sArray = new ArrayList<OverlappedGenomicElement>(6);
		sArray.add(new OverlappedGenomicElement(s1));
		sArray.add(new OverlappedGenomicElement(s2));
		sArray.add(new OverlappedGenomicElement(s3));
		sArray.add(new OverlappedGenomicElement(s4));
		sArray.add(new OverlappedGenomicElement(s5));
		sArray.add(new OverlappedGenomicElement(s6));
		
		Collections.sort(sArray);
		
		OverlappedCollectionStream gStream = new OverlappedCollectionStream(gArray);
		OverlappedCollectionStream eStream = new OverlappedCollectionStream(eArray);
		OverlappedCollectionStream sStream = new OverlappedCollectionStream(sArray);
		
		AllOverlappedElements firstObj = new AllOverlappedElements(sStream, eStream);
		firstObj.fillTreeSet();
		 GenewiseSnpWeights myObjView= firstObj.getAsGenewiseSnpWeightDataStructure(0.4, 0);
		 
		 HashSet<String> myIds = myObjView.returnAllSnpIds();
		 assertTrue(myIds.contains("s1"));
		 assertTrue(myIds.contains("s1"));
		 			
		 TreeSet<OverlappedGenomicElement> myOverlappedGenomicRegions = firstObj.getMyLargeList();		
		OverlappedCollectionStream myOverlappedGenomicRegionStream = new OverlappedCollectionStream(myOverlappedGenomicRegions);
		AllOverlappedElements secObj = new AllOverlappedElements(myOverlappedGenomicRegionStream, gStream);
		secObj.fillTreeSet();
		GenewiseSnpWeights secObjView= secObj.getAsGenewiseSnpWeightDataStructure(0.4, 1);
		 HashSet<String> myIds2 = myObjView.returnAllSnpIds();
		 
//		 assertTrue(myIds2.contains("s1"));
//		 assertTrue(myIds2.contains("s2"));
//		 assertTrue(myIds2.contains("s3"));
//		 assertTrue(!myIds2.contains("s4"));
//		 assertTrue(myIds2.contains("s5"));
//		 assertTrue(!myIds2.contains("s6"));
		 
		 Snp newSnp = new Snp("idnotIn");
		 newSnp.setMaf(0.5);
		 LinkedList<Snp> SnpList = new LinkedList<Snp>();
		 SnpList.add(newSnp);
		 SnpList.add(s1);
		 SnpList.add(s2);
		 SnpList.add(s3);
		 SnpList.add(s4);
		 SnpList.add(s5);		 
		 SnpList.add(s6);
		 		 
		 SnpWeightPairs overlappedList =  secObjView.getOverlappedSnpsWithWeights("g1", SnpList);
		 ArrayList<Snp>  finalList = overlappedList.getSnps();
		 double[] finalWeights = overlappedList.getWeightsNotList();
		 for (double weight : finalWeights){
			 assertTrue(weight==0.4);
		 }
		 assertTrue(finalList.contains(s1));
		 assertTrue(finalList.contains(s2));
		 
		 assertTrue(!finalList.contains(s5));		 
		 assertTrue(!finalList.contains(s3));
		 assertTrue(!finalList.contains(s4));
		 assertTrue(!finalList.contains(s6));
		 assertTrue(!finalList.contains(newSnp));
		 
		 overlappedList =  secObjView.getOverlappedSnpsWithWeights("g2", SnpList);
		 finalList = overlappedList.getSnps();
		 assertTrue(!finalList.contains(s1));
		 assertTrue(finalList.contains(s2));
		 assertTrue(!finalList.contains(s5));		 
		 assertTrue(!finalList.contains(s3));
		 assertTrue(!finalList.contains(s4));
		 assertTrue(!finalList.contains(s6));
		assertTrue(!finalList.contains(newSnp));
		
		overlappedList =  secObjView.getOverlappedSnpsWithWeights("g3", SnpList);
		 finalList = overlappedList.getSnps();
		 assertTrue(!finalList.contains(s1));
		 assertTrue(!finalList.contains(s2));
		 assertTrue(finalList.contains(s5));		 
		 assertTrue(finalList.contains(s3));
		 assertTrue(!finalList.contains(s4));
		 assertTrue(!finalList.contains(s6));
		assertTrue(!finalList.contains(newSnp));
		 
		 
		 
		 
	}
}