package ch.unil.genescore.prioritization.test;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.collect.TreeMultimap;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.test.ParserStub;
import ch.unil.genescore.prioritization.SparseNet;

public class SparseNetTest {

	@BeforeClass
	public static void testSetup() {

	}

	@AfterClass
	public static void testCleanup() {
	}

	// ============================================================================
	// TESTS

	@Test
	public void loadNetworkDataTestLoadingFails() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t1");
		myDat.add("g2\tg1");
		myDat.add("g2\tg1\t1");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet net = new SparseNet();
		net.setParser(parser);
		boolean loadFail = false;
		try {
			net.loadNetworkData();
		} catch (RuntimeException e) {
			loadFail = e.getMessage().equals("Exception during parsing of network: not 3 columns in network file");

		}
		assertTrue(loadFail);
	}

	@Test
	public void loadNetworkDataTestLoadingFails2() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t1");
		myDat.add("g2\tg1\t0.3");
		myDat.add("g2\tg1\t2");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet net = new SparseNet();
		net.setParser(parser);
		boolean loadFail = false;
		try {
			net.loadNetworkData();
		} catch (RuntimeException e) {
			loadFail = e
					.getMessage()
					.equals("Exception during parsing of network: Value is larger than one. One is the allowed maximum. Scale network externally before trying again.");

		}
		assertTrue(loadFail);
	}

	@Test
	public void loadNetworkDataTestLoadingWorks() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t1");
		myDat.add("g2\tg1\t0.4");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet net = new SparseNet();
		net.setParser(parser);
		boolean loadWorks = false;
		// try{
		net.loadNetworkData();
		System.out.print("asdfs");
		// loadWorks= true;
		// }catch(RuntimeException e){
		// }
		// assertTrue(loadWorks);
	}

	@Test
	public void returnReverseOrderedConnectionTreeTestFailsOrIsEmpty() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t0.2");
		myDat.add("g2\tg3\t1");
		myDat.add("g2\tg1\t0.1");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet net = new SparseNet();
		net.setParser(parser);
		net.loadNetworkData();
		HashSet<Gene> myGenes = new HashSet<Gene>();

		Gene gene1 = new Gene("id", "g1");
		myGenes.add(gene1);
		myGenes.add(new Gene("id", "g1"));
		myGenes.add(new Gene("id", "g2"));
		myGenes.add(new Gene("id", "g3"));
		net.setGenes(myGenes);
		boolean fail = false;
		try {
			net.returnReverseOrderedConnectionTree(new Gene("id", "g1"));
		} catch (RuntimeException e) {
			fail = e.getMessage()
					.equals("Exception: query gene has same symbol as some gene in provided gene list but its not the same underlying gene object.");
		}
		assertTrue(fail);
		assertTrue(net.returnReverseOrderedConnectionTree(new Gene("g4"))
				.isEmpty());
		fail = false;
		try {
			TreeMultimap<Double, Gene> myTree = net
					.returnReverseOrderedConnectionTree(gene1);
		} catch (RuntimeException e) {
			fail = e.getMessage().equals(
					"Double Entry in parsed network or connections to itself. Gene g1 involved. remove double entries");					
		}
		assertTrue(fail);
	}

	@Test
	public void returnReverseOrderedConnectionTreeTestWorks() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t0.1");
		myDat.add("g2\tg3\t1.0");
		myDat.add("g2\tg4\t0.1");
		myDat.add("g4\tg1\t0.2");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet net = new SparseNet();
		net.setParser(parser);
		net.loadNetworkData();
		HashSet<Gene> myGenes = new HashSet<Gene>();

		Gene gene1 = new Gene("id", "g1");
		myGenes.add(gene1);
		Gene gene2 = new Gene("id2", "g2");
		myGenes.add(gene2);
		Gene gene3 = new Gene("id3", "g3");
		myGenes.add(gene3);
		Gene gene4 = new Gene("id4", "g4");
		myGenes.add(gene4);
		net.setGenes(myGenes);
		System.out.println("sadfs");
		TreeMultimap<Double, Gene> out = net.returnReverseOrderedConnectionTree(gene1);
		int count = 0;
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene1);
			if (count == 1)
				assertTrue(entry.getKey() == 0.2 && entry.getValue() == gene4);
			if (count == 2)
				assertTrue(entry.getKey() == 0.1 && entry.getValue() == gene2);
			count++;
		}
		
		out = net.returnReverseOrderedConnectionTree(gene2);
		count = 0;
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene2);
			if (count == 1)
				assertTrue(entry.getKey() == 1 && entry.getValue() == gene3);
			if (count == 2)
				assertTrue(entry.getKey() == 0.1 && entry.getValue() == gene1);
			if (count == 3)
				assertTrue(entry.getKey() == 0.1 && entry.getValue() == gene4);
			count++;
		}

	}
}
