package ch.unil.genescore.vegas.test;




import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GeneResultsSnpsOutOfBounds;


public class GeneResultsSnpsOutOfBoundsTest {

		
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();		
	}

	@AfterClass
	public static void testCleanup() { 
	}
	
	
		
	@Test		
	public void writeResultsToFileTest(){
					
		FileExportStub stub = new FileExportStub();
		LinkedHashMap<Gene, Integer> map = new LinkedHashMap<Gene, Integer>();
		map.put(new Gene("g1"), 0);
		map.put(new Gene("g2"), 2);
		
		
		//String=results.getStrings().get(0);
		GeneResultsSnpsOutOfBounds results = new GeneResultsSnpsOutOfBounds();
		results.setExporter(stub);
		results.setZeroSnpsOrAboveSnpLimit(map);
		results.writeResultsToFile("blub");
		String str=stub.getStrings().get(0);
		System.out.print(str);
		ArrayList<String> strs = stub.getStrings();
		assertTrue(stub.getStrings().get(0).equals("gene_id\tsymbol\tSNPs"));		
		assertTrue(stub.getStrings().get(1).equals("g1\tNA\t0"));
		assertTrue(stub.getStrings().get(2).equals("g2\tNA\t2"));
		assertTrue(stub.isclosed());
	}
	
	
}
