package ch.unil.genescore.pathway.test;
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
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Collection;

import javastat.inference.nonparametric.RankSumTest;
import jsc.independentsamples.MannWhitneyTest;
import jsc.tests.H1;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.pathway.GeneScoreList;
import ch.unil.genescore.pathway.GeneSet;
import ch.unil.genescore.pathway.GeneSetLibrary;
import ch.unil.genescore.pathway.PathwayMain;


/**
 * Unit tests for GenomicElement.
 */
public class RankSumTestTest {
	
	
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

	/** Test PathwayMain.run() */
	@Test
	public void testRankSumTest(){
	double [] testdata1 = {0.8, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46};
	double [] testdata2 = {1.15, 0.88, 0.9, 0.74, 1.21};

	

	// Null constructor
	RankSumTest testclass2 = new RankSumTest();
	double pValue = testclass2.pValue("equal", testdata1, testdata2); 	
	assertEquals(pValue, 0.2544, 0.001);//copied from R
	pValue = testclass2.pValue("greater", testdata1, testdata2); 	
	assertEquals(pValue, 1-0.1272, 0.001);//copied from R
	
	
	
	MannWhitneyUTest test = new MannWhitneyUTest();
	
		double pval = test.mannWhitneyUTest(testdata1,testdata2);
		System.out.println("");
		MannWhitneyTest test3 = new	MannWhitneyTest(testdata1,testdata2);
		double pval3=test3.approxSP();		
		double pval4=test3.exactSP();
		assertEquals(pval3, 0.2446, 0.001);
		assertEquals(pval4, 0.2544, 0.001);
		
		MannWhitneyTest test4 = new	MannWhitneyTest(testdata1,testdata2, H1.GREATER_THAN);
		double pval5=test4.approxSP();		
		double pval6=test4.exactSP();
		assertEquals(pval5, 0.1227, 0.001);
		assertEquals(pval6, 0.1272, 0.001);
		
		
		System.out.println("");
		
		

	}
}
