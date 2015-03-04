/*
Copyright (c) 2013 Daniel Marbach

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper available at:
http://networkinference.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
package ch.unil.genescore.vegas.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpWeightMap;

import org.apache.commons.math3.distribution.NormalDistribution;


//These tests do not cover much just late additions
public class GenomeWideScoringTest {

	
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test test statistic */
	@Test
	public void removLowMafSnpsTest() {
		SnpWeightMap myMap = new SnpWeightMap();
		String[] inputString = {"id1;3", "id2;1", "id3;4"};
		ArrayList<Snp> snpList = new ArrayList<Snp>();
		Snp snp1 = new Snp("id2");		
		
		snpList.add(snp1);
		byte[] fakeGenotype = {0,0,1,1,1,1,0,1};
		snp1.setGenotypes(fakeGenotype);
		 snp1.computeAlleleStats();
		 Snp.setGenotypeLength(8);
		 
		 
		Snp snp2 = new Snp("id5");
			
		snpList.add(snp2);
		byte[] fakeGenotype2 = {0,0,1,1,0,0,0,0};
		snp2.setGenotypes(fakeGenotype2);
		snp2.computeAlleleStats();
		
		Snp snp3 = new Snp("id3");
		
		snpList.add(snp3);
		byte[] fakeGenotype3 = {0,0,0,0,1, 1, 0,1};
		snp3.setGenotypes(fakeGenotype3);
		snp3.computeAlleleStats();
		
		Settings.useMafCutoff_=0.0;
		GenomeWideScoring GS = new GenomeWideScoring();
		GS.removeLowMafSnps(snpList);
		assertTrue(snpList.get(1)==snp2);
		System.out.println("sfd");
		Settings.useMafCutoff_=0.2;
		GS.removeLowMafSnps(snpList);
		assertTrue(snpList.get(1)==snp2);
		Settings.useMafCutoff_=0.27;
		GS.removeLowMafSnps(snpList);
		assertTrue(snpList.get(1)==snp3);
		
		
	}
}