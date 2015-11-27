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

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Pascal;
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
