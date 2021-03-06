/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
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
