package ch.unil.genescore.vegas.test;

import java.util.ArrayList;
import java.util.Collections;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeSet;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.vegas.AllOverlappedElements;
import ch.unil.genescore.vegas.BedFileStream;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.OverlappedCollectionStream;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import ch.unil.genescore.vegas.Snp;



/**
 * Unit tests for GenewiseSnpWeightsTest
 */
public class BedFileStreamTest{

	
	static double delta;
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {	
		delta=1E-14;
	}

	@AfterClass
	public static void testCleanup() { }
	
	@Rule
    public ExpectedException thrown= ExpectedException.none();
	
	// ============================================================================
	// TESTS

	/** Test PathwayMain.run() */
	@Test
	public void  RunBedFileStream1Test() {
		//test sorted file;
		BedFileStream myBedStream = new BedFileStream("resources/test/vegas/StreamFilesTest/sortedBed.bed");
		OverlappedGenomicElement el = myBedStream.getNextAsOverlappedGenomicElement();				
		assertTrue(el.getMainElement().chr_.equals("chr1"));
		el = myBedStream.getNextAsOverlappedGenomicElement();				
		assertTrue(el.getMainElement().start_==2);
		el = myBedStream.getNextAsOverlappedGenomicElement();				
		assertTrue(el.getMainElement().start_==2);
		while(myBedStream.streamOpen()){
			el = myBedStream.getNextAsOverlappedGenomicElement();				
		}
		assertTrue(el.getMainElement().chr_.equals("chrY"));
		assertTrue(!myBedStream.streamOpen());
		thrown.expect(RuntimeException.class);
		OverlappedGenomicElement el2 = myBedStream.getNextAsOverlappedGenomicElement();
		}	
	
	@Test	
	public void  RunBedFileStream2Test() {
		//test unsorted file;
		BedFileStream myBedStream = new BedFileStream("resources/test/vegas/StreamFilesTest/unsortedBed.bed");
		OverlappedGenomicElement el = null;
		thrown.expect(RuntimeException.class);
		while(myBedStream.streamOpen()){
			el = myBedStream.getNextAsOverlappedGenomicElement();				
		}		
	}		
}
