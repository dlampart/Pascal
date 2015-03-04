package ch.unil.genescore.eqtlOverlap.test;

import java.util.ArrayList;

import static org.junit.Assert.*;
import no.uib.cipr.matrix.DenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.eqtlOverlap.EqtlEvaluator;
import ch.unil.genescore.eqtlOverlap.EqtlGenePair;
import ch.unil.genescore.eqtlOverlap.EqtlResults;
import ch.unil.genescore.eqtlOverlap.GeneDataWithEqtl;
import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.test.ParserStub;

public class EqtlResultsTest {
	
	class EqtlResultsMock extends EqtlResults {
		
		@Override
		protected FileParser createParser(String filename){
			return(new ParserStub(filename));
		}
	} 
	
	// SETUP
	
		@BeforeClass
		public static void testSetup() {
		//	Settings.loadSettings();
		}

		@AfterClass
		public static void testCleanup() { }
		
		@Test
		public void test(){
			
			ArrayList<String> myDat= new ArrayList<String>();
			myDat.add("rs1\tC\tA\tNA\t77.58643142977844\t9.813427854297537E-198\tANKDD1A");
			myDat.add("rs2\tC\tA\tNA\t77.58643142977844\tee\tANKDD1A");
			myDat.add("rs2\tC\tA\tNA\t77.58643142977844\t22\tANKDD1B");
			myDat.add("rs1\tC\tA\tNA\t77.58643142977844\t2342\tANKDD1A");
			ParserStub.setTestData(myDat);
			
			EqtlResults unit = new EqtlResultsMock();
			unit.loadEqtlFile("anyName");
			assertTrue(unit.eqtlGenePairExists("ANKDD1A"));
			assertTrue(unit.eqtlGenePairExists("ANKDD1B"));
			EqtlGenePair gp =unit.getEqtlGenePair("ANKDD1A");
			assertEquals(gp.getSnpList().get("rs1").getPval(),2342,0.00001);// uses latest entry
			assertTrue(gp.getSnpList().get("rs1").getMajorAllele()=='C');
			assertTrue(gp.getSnpList().get("rs1").getMinorAllele()=='A');
			assertEquals(gp.getSnpList().get("rs1").getZscore(),77.586,0.001);
			assertTrue(!gp.getSnpList().containsKey("rs2"));
			
		}
}
