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

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.GeneWithItsSnps;
import ch.unil.genescore.vegas.Snp;

public class GeneWithItsSnpsTest {

	// SETUP
	
			@BeforeClass
			public static void testSetup() {
				Settings.loadSettings();
			}

			@AfterClass
			public static void testCleanup() {
			
			
				
			
			}									
			
			
			private class snpAllwaysCoding extends Snp {
				 public snpAllwaysCoding(String id) {
						super(id);						
				 }
				 public boolean isCodingForOtherGene(String thisGene){return true;}
			}
				 
				private class snpNeverCoding extends Snp {
					 public snpNeverCoding(String id) {
						super(id);
					}
					public boolean isCodingForOtherGene(String thisGene){return false;}
				}
			
			// ============================================================================
			// TESTS
					
			@Test
			public void constructorTest1(){
									

				Settings.removeCodingSnpsOfOtherGenes_=true;
				 ArrayList<Snp> mySnps = new ArrayList<Snp>();
				 mySnps.add(new snpNeverCoding("id1"));
				 mySnps.add(new snpAllwaysCoding("id2"));
				 mySnps.add(new snpNeverCoding("id3"));				 
				 GeneWithItsSnps myG = new GeneWithItsSnps(new Gene("gId"),mySnps);
				 assertTrue(myG.getGene().getId().equals("gId"));
				 
				 assertTrue(myG.getNrOfSnps()==2);
				 assertTrue(myG.getSnpList().get(0).id_.equals("id1"));
				 assertTrue(myG.getSnpList().get(1).id_.equals("id3"));				 				 
				 
			}
			
			@Test
			public void constructorTest2(){
									

				Settings.removeCodingSnpsOfOtherGenes_=false;
				 ArrayList<Snp> mySnps = new ArrayList<Snp>();
				 mySnps.add(new snpNeverCoding("id1"));
				 mySnps.add(new snpAllwaysCoding("id2"));
				 mySnps.add(new snpNeverCoding("id3"));				 
				 GeneWithItsSnps myG = new GeneWithItsSnps(new Gene("gId"),mySnps);
				 assertTrue(myG.getGene().getId().equals("gId"));
				 
				 assertTrue(myG.getNrOfSnps()==3);
				 assertTrue(myG.getSnpList().get(0).id_.equals("id1"));
				 assertTrue(myG.getSnpList().get(1).id_.equals("id2"));
				 assertTrue(myG.getSnpList().get(2).id_.equals("id3"));				 				 
				 
			}
						
				
			
}
