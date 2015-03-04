package ch.unil.genescore.prioritization.test;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.test.ParserStub;
import ch.unil.genescore.prioritization.PrunedConnectionTreeGetter;
import ch.unil.genescore.prioritization.SparseNet;

import com.google.common.collect.TreeMultimap;

public class PrunedConnectionTreGetterTest {

	@BeforeClass
	public static void testSetup() {

	}

	@AfterClass
	public static void testCleanup() {
	}

	// ============================================================================
	// TESTS
	
//	@Test
	public void returnReverseOrderedConnectionTreeTestFailsOrIsEmpty() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t0.2");
		myDat.add("g2\tg3\t1");
		myDat.add("g2\tg1\t0.1");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet netInternal = new SparseNet();
		netInternal.setParser(parser);
		PrunedConnectionTreeGetter net = new PrunedConnectionTreeGetter();
		net.setConnectionTreeGetter(netInternal);
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
					"Double Entry in parsed network. remove double entries");
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
		SparseNet netInternal = new SparseNet();
		netInternal.setParser(parser);
		PrunedConnectionTreeGetter net = new PrunedConnectionTreeGetter();
		net.setConnectionTreeGetter(netInternal);
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
		assertTrue(out.entries().size()==1);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene1);
		}
		
		out = net.returnReverseOrderedConnectionTree(gene2);
		count = 0;
		assertTrue(out.entries().size()==1);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene2);		
		}

	}
	
	@Test
	public void returnReverseOrderedConnectionTreeTestWorks2() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t0.1");
		myDat.add("g2\tg3\t1.0");
		myDat.add("g2\tg4\t0.1");
		myDat.add("g4\tg1\t0.2");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet netInternal = new SparseNet();
		netInternal.setParser(parser);
		PrunedConnectionTreeGetter net = new PrunedConnectionTreeGetter();
		net.setConnectionTreeGetter(netInternal);
		net.loadNetworkData();
		net.setExt(10000);
		HashSet<Gene> myGenes = new HashSet<Gene>();

		Gene gene1 = new Gene("id", "g1");
		myGenes.add(gene1);
		Gene gene2 = new Gene("id2", "g2");
		myGenes.add(gene2);
		Gene gene3 = new Gene("id3", "g3");
		myGenes.add(gene3);
		Gene gene4 = new Gene("id4", "g4");
		myGenes.add(gene4);
		
		gene1.setPosition("A", 100, 1000, true);
		gene2.setPosition("A", 10050, 20000, true);
		gene3.setPosition("B", 20000, 1000, true);
		gene4.setPosition("B", 100, 10500, true);
		
		net.setGenes(myGenes);
		System.out.println("sadfs");
		TreeMultimap<Double, Gene> out = net.returnReverseOrderedConnectionTree(gene1);
		int count = 0;
		assertTrue(out.entries().size()==2);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene1);
			if (count == 1)
				assertTrue(entry.getKey() == 0.2 && entry.getValue() == gene4);
			count++;
		}
		
		out = net.returnReverseOrderedConnectionTree(gene2);
		count = 0;
		assertTrue(out.entries().size()==2);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene2);		
			if (count == 1)
				assertTrue(entry.getKey() == 1.0 && entry.getValue() == gene3);			
			count++;
		}

	}
	
	@Test
	public void returnReverseOrderedConnectionTreeTestWorks3() {
		ArrayList<String> myDat = new ArrayList<String>();
		myDat.add("g1\tg2\t0.1");
		myDat.add("g1\tg3\t0.1");
		myDat.add("g1\tg4\t0.1");
		myDat.add("g2\tg3\t0.1");
		myDat.add("g2\tg4\t0.1");
		myDat.add("g3\tg4\t0.1");
		ParserStub.setTestData(myDat);
		ParserStub parser = new ParserStub();
		SparseNet netInternal = new SparseNet();
		netInternal.setParser(parser);
		PrunedConnectionTreeGetter net = new PrunedConnectionTreeGetter();
		net.setConnectionTreeGetter(netInternal);
		net.loadNetworkData();
		net.setExt(10000);
		HashSet<Gene> myGenes = new HashSet<Gene>();

		Gene gene1 = new Gene("id", "g1");
		myGenes.add(gene1);
		Gene gene2 = new Gene("id2", "g2");
		myGenes.add(gene2);
		Gene gene3 = new Gene("id3", "g3");
		myGenes.add(gene3);
		Gene gene4 = new Gene("id4", "g4");
		myGenes.add(gene4);
		
		gene1.setPosition("A", 100, 1000, true);
		gene2.setPosition("A", 10050, 20000, true);		
		gene3.setPosition("A", 29000, 40500, true);
		gene4.setPosition("A", 11100, 20000, true);
		
		net.setGenes(myGenes);
		System.out.println("sadfs");
		TreeMultimap<Double, Gene> out = net.returnReverseOrderedConnectionTree(gene1);
		int count = 0;
		assertTrue(out.entries().size()==2);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene1);
			if (count == 1)
				assertTrue(entry.getKey() == 0.1 && entry.getValue() == gene4);
			count++;
		}
		
		out = net.returnReverseOrderedConnectionTree(gene2);
		count = 0;
		assertTrue(out.entries().size()==1);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene2);				
			count++;
		}
		
		out = net.returnReverseOrderedConnectionTree(gene3);
		count = 0;
		assertTrue(out.entries().size()==2);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene3);				
			if (count == 1)
			assertTrue(entry.getKey() == 0.1 && entry.getValue() == gene1);				
			count++;
		}
		
		out = net.returnReverseOrderedConnectionTree(gene4);
		count = 0;
		assertTrue(out.entries().size()==2);
		for (Map.Entry<Double, Gene> entry : out.entries()) {
			if (count == 0)
				assertTrue(entry.getKey() == 1.1 && entry.getValue() == gene4);				
			if (count == 1)
			assertTrue(entry.getKey() == 0.1 && entry.getValue() == gene1);				
			count++;
		}

	}

}

