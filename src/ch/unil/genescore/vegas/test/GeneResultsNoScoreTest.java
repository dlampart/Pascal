/*******************************************************************************
 * Copyright (c) 2015 IBM Corporation and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     IBM Corporation - initial API and implementation
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
import ch.unil.genescore.vegas.GeneResultsNoScore;


public class GeneResultsNoScoreTest {

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
			GeneResultsNoScore results = new GeneResultsNoScore();
			results.add("g1");
			results.add("g2");
			
			results.setExporter(stub);
			
			results.writeResultsToFile("blub");
			String str=stub.getStrings().get(0);
			System.out.print(str);
			ArrayList<String> strs = stub.getStrings();		
			assertTrue(stub.getStrings().get(0).equals("chromosome\tstart\tend\tstrand\tgene_id\tsymbol\tScore\tStatus"));	
			assertTrue(stub.getStrings().get(1).equals("g1"));
			assertTrue(stub.getStrings().get(2).equals("g2"));
		}

}
