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
package ch.unil.genescore.gene.test;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeSet;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Pascal;


/**
 * Unit tests for GenomicElement.
 */
public class GenomicElementTest {
	
	/** The Pascal instance (initializes Settings) */
	private static Pascal psc = new Pascal();

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Pascal.set.resetToDefaults();
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test test statistic */
	@Test			
	public void compareToTest(){
		// test whether overlap chromosome works as expected.
		Gene g1 = new Gene("g1");				
		Gene g2 = new Gene("g2");			
		
	
		g1.setPosition("chr22", 12, 500, true);
		g2.setPosition("chr22", 20, 500, true);
		assertTrue(g1.compareTo(g2)<0);				
		g1.setPosition("chr19", 20, 500, true);
		g2.setPosition("chr2", 20, 500, true);
		assertTrue(g1.compareTo(g2)<0);
		g1.setPosition("chr1", 20, 500, true);
		g2.setPosition("chr2", 1, 1200, true);
		assertTrue(g1.compareTo(g2)<0);
		g1.setPosition("chr2", 20, 500, true);
		g2.setPosition("chr2", 20, 1200, true);
		assertTrue(g1.compareTo(g2)<0);	
		g1.setPosition("chr2", 20, 500, true);
		g2.setPosition("chr2", 20, 500, true);
		assertTrue(g1.compareTo(g2)<0);		
	}	
	/** Test test statistic */
	@Test			
	public void completelyOverlapsTest(){
		
		Gene g1 = new Gene("g1");				
		Gene g2 = new Gene("g2");			
			
		g1.setPosition("chr22", 12, 500, true);
		g2.setPosition("chr22", 20, 500, true);
		assertTrue(g1.completelyOverlapsElement(g2, 0, 0));		
			
		g1.setPosition("chr22", 12, 500, true);
		g2.setPosition("chr22", 400, 600, true);
		assertTrue(!g1.completelyOverlapsElement(g2,0,0));		
		g1.setPosition("chr22", 12, 500, true);
		g2.setPosition("chr21", 12, 50, true);
		assertTrue(!g1.completelyOverlapsElement(g2,0,0));	
		
		g1.setPosition("chr22", 12, 500, true);
		g2.setPosition("chr22", 400, 600, true);
		assertTrue(g1.completelyOverlapsElement(g2,0,100));		
		assertTrue(g2.completelyOverlapsElement(g1,-400,-50));		
		
	}
	
	/** Test test statistic */
	@Test			
	public void partiallyOverlapsTest(){
		
		Gene g1 = new Gene("g1");				
		Gene g2 = new Gene("g2");			
			
		g1.setPosition("chr22", 400, 600, true);
		g2.setPosition("chr22", 20, 500, true);
		assertTrue(g1.partiallyOverlapsElement(g2));
		g1.setPosition("chr22", 400, 600, true);
		g2.setPosition("chr22", 200, 800, true);
		assertTrue(g1.partiallyOverlapsElement(g2));		
		
		g1.setPosition("chr22", 500, 800, true);
		g2.setPosition("chr22", 500, 600, true);
		assertTrue(g1.partiallyOverlapsElement(g2));		
		
		g1.setPosition("chr22", 60, 600, true);
		g2.setPosition("chr22", 20, 50, true);
		assertTrue(!g1.partiallyOverlapsElement(g2));
		g1.setPosition("chr22", 12, 500, true);
		g2.setPosition("chr22", 400, 600, true);
		assertTrue(g1.partiallyOverlapsElement(g2));		
		g1.setPosition("chr22", 60, 500, true);
		g2.setPosition("chr21", 12, 50, true);
		assertTrue(!g1.partiallyOverlapsElement(g2));		
		
	}
	
	/** Test test statistic */
	@Test			
	public void isPastTest(){
		
		Gene g1 = new Gene("g1");				
		Gene g2 = new Gene("g2");			
		
		g1.setPosition("chr22", 20, 500, true);	
		g2.setPosition("chr22", 12, 500, true);
		
		assertTrue(!g1.isPastElement(g2,0));		
		assertTrue(!g2.isPastElement(g1,0));	
		
		g1.setPosition("chr19", 12, 500, true);
		g2.setPosition("chr2", 20, 500, true);
		assertTrue(!g1.isPastElement(g2,0));		
		assertTrue(g2.isPastElement(g1,0));
		
		g1.setPosition("chr2", 12, 400, true);
		g2.setPosition("chr2", 500, 600, true);
		assertTrue(!g2.isPastElement(g1,-110));		
		assertTrue(g2.isPastElement(g1,-90));		
		
		
	}
	
	
	/** Test behavior of compareTo(), equals(), hashCode() */
	//@Test Modify
	//TODO: discuss change to lexical order (arguments: 1 interoperability:Bedops, 2 extensibility to
	//other organisms) 
	@Test
	public void testCompareToLexical() {
	
		Gene e1 = new Gene("e1");
		Gene e2 = new Gene("e2");
		Gene e3 = new Gene("e3");
		Gene e4 = new Gene("e4");
		Gene e5 = new Gene("e5");
		Gene e6 = new Gene("e6");
		Gene e7 = new Gene("e7");
		Gene e8 = new Gene("e8");
		
		Gene e9 = new Gene("e9");
	//	Gene e9b = new Gene("e9b");
		Gene e10 = new Gene("e10");		
		Gene e10b = new Gene("e10");		
		Gene e11 = new Gene("e11");
		
		//GenomicElement e12 = new GenomicElement("e12");
		
		e1.setPosition("chr1", 1,1, true);
		e2.setPosition("chr1", 2,2,true);
		e3.setPosition("chr1", 2,5, true);
		e4.setPosition("chr12",2, 12, true);
		e5.setPosition("chr22", 2, 22, true);
		e6.setPosition("chr9", 1, 20, true);// same pos, different id	
		e7.setPosition("chr9", 1, 20, true);// same pos, different id
		e8.setPosition("chr9", 30,60, true);	
		e9.setPosition("chrA", 1, 890, true);
		e10.setPosition("chrM", 1, 30, true); // same pos, same id	
		e10b.setPosition("chrM", 1, 30, false); // same pos, same id	
		//! lexical id comparison
		e11.setPosition("chrM", 1, 30, true);// same pos, different id !!!! lexical id comparison
		
		
		// This gives indeed an exception "unknown chromosome" when adding to treeset
		//e12.setPosition("chr23", 1); 
	
		// Add in random order
		ArrayList<Gene> list = new ArrayList<Gene>();
		list.add(e5);
		list.add(e1);
		list.add(e11);
		list.add(e6);
		list.add(e9);
		list.add(e10);
		list.add(e10b);
		list.add(e4);
		list.add(e2);
		list.add(e3);
		list.add(e8);
		list.add(e7);
		
		// And shuffle to be sure :)
		Collections.shuffle(list);
				
		TreeSet<Gene> treeSet = new TreeSet<Gene>();
		treeSet.addAll(list);
	
		// Shuffle the list again
		Collections.shuffle(list);
		// Sort it
		Collections.sort(list);
		System.out.println(list.size());
		System.out.println(treeSet.size());
		assertEquals(12, list.size());
		assertEquals(11, treeSet.size());
		
		// Check order in arraylist
		for (int i=0; i<12; i++) {
			if (i < 10)
				assertEquals("e" + (i+1), list.get(i).id_);
			else
				assertEquals("e" + (i), list.get(i).id_);
			
			
		}
		
		// Check order in treeset
		int counter = 1;
		for (GenomicElement el : treeSet){
			assertEquals("e" + counter++, el.id_);
		}

	}
}