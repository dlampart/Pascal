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
package ch.unil.genescore.pathway.test;


import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Collection;

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
public class PathwayMainTest {
	
	
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
	//@Test
	//public void testRun() {
	//	testRun(false);
	//}
	
	/** Test PathwayMain.run() when laoding pre-computed gene scores */
	@Test
	public void testRunWithLoadedScores() {
		testRun(true);
	}


	@SuppressWarnings({ "unchecked", "rawtypes" })
	/** Test PathwayMain.run() */
	private void testRun(boolean loadScores) {
		
		Settings.chromosome_ = "";
		Settings.mergeGenesDistance_ = 1;
		Settings.geneSetFile_ = "resources/test/pathway/testGeneSets.txt";
		Settings.genomeAnnotation_ = "ucsc";
		Settings.snpPvalFile_ = "resources/test/vegas/EUR.wtccc2_ms.pvals_real.txt";
		Settings.withZScore_ = false;
		if (loadScores)
			Settings.geneScoreFile_ = "resources/test/pathway/EUR.wtccc2_ms.pvals_real--testGeneSets.genescores.txt";
		else
			Settings.geneScoreFile_ = "";
		
		// Load snps, pvals, genes and gene sets
		PathwayMain pathwayMain = new PathwayMain();
		pathwayMain.run();
		GeneSetLibrary geneSetLib = pathwayMain.getGeneSetLib();
		ArrayList<GeneSet> geneSets = geneSetLib.getGeneSets();
		assertEquals(2, geneSets.size());
		assertEquals(loadScores ? 6 : 9, geneSetLib.getGenes().size());

		GeneSet set1 = geneSets.get(0);
		GeneSet set2 = geneSets.get(1);
		
		// Check that genes have been correctly loaded and sorted
		String str1 = set1.toString();
		String str2 = set2.toString();
		if (loadScores) {
			assertEquals("57449	727751	113540	146225", str1);
			assertEquals("6642	1009", str2);
		} else {
			assertEquals("729737	729759	26155	57449	727751	113540	146225", str1);
			assertEquals("729737	6642	1009", str2);
		}
		
		// Run pathway enrichment
		Settings.useSimulation_ = true;
		Settings.useChi2_ = false;
		//pathwayMain.run();

		// Check meta-genes
		GeneSet metaGenes = new GeneSet("metaGenes", (Collection)geneSetLib.getMetaGenes());
		if (loadScores)
			assertEquals("meta_113540_146225", metaGenes.toString());
		else
			assertEquals("meta_729737_729759_26155	meta_113540_146225", metaGenes.toString());

		// Check final genes in sets after removing genes without scores
		str1 = set1.toString();
		str2 = set2.toString();
		if (loadScores) {
			assertEquals(str1, "57449	727751	meta_113540_146225");
			assertEquals(str2, "6642	1009");
		} else {
			assertEquals(str1, "57449	727751	meta_113540_146225");
			assertEquals(str2, "6642	1009");
		}
		
		ArrayList<Gene> genes = Gene.sortByPosition(geneSetLib.getGenes());
		String[] ids = { "57449", "6642", "727751", "1009", "meta_113540_146225" }; // Sorted by genomic coordinates
		//String[] ids = { "6642", "meta_113540_146225", "57449", "1009", "727751" }; // Sortyed by gene score
		
		// This is copy-pasted from running this code, so only a check that nothing got broken, not that it was correct initially
		double[] geneScores = { 2.85178838E-1, 1.11902529E-1, 4.56584036E-1, 3.71473436E-1, 2.70346992E-1 }; // Sorted by genomic coordinates
		//double[] geneScores = { 1.11902529E-1, 2.70346992E-1, 2.85178838E-1, 3.71473436E-1, 4.56584036E-1 }; // Sorted by gene score
		
		int count = 0;
		for (Gene g : genes) {
			assertEquals(ids[count], g.id_);
			double[] score = g.getScore();
			assertEquals(geneScores[count], score[0], 1e-6);
			count++;
		}
		
		// Check rank-sum enrichment
		// Copy-pasted from running this code... Different from R, but I guess it's just because there are very few samples.
		// I did do another test case with many samples and result was consistent with R.
		double pval1 = set1.getRankSumPvalue();
		double pval2 = set2.getRankSumPvalue();
		assertEquals(0.5637028616507731, pval1, 1e-12);
		assertEquals(0.5637028616507731, pval2, 1e-12);
		
		// Check chi2 enrichment
		Settings.useChi2_ = true;
		geneSetLib.computeEnrichment();
		genes = new GeneScoreList(geneSetLib.getGenes(), true).getGenes();
		
		double[] geneScoresChi2 = { 1.64237442, 0.70832630, 0.27499590, 0.06418475, 0.00000000 }; // Sorted by gene score

		count = 0;
		for (Gene g : genes) {
			//assertEquals(ids[count], g.id_);
			double score = g.getChi2Stat();
			assertEquals(geneScoresChi2[count], score, 1e-6);
			count++;
		}
		
		pval1 = set1.getChi2Pvalue();
		pval2 = set2.getChi2Pvalue();
		assertEquals(0.1947126, pval1, 1e-6);
		assertEquals(1-0.5739845, pval2, 1e-6);
		assertEquals(true, set1.getDepletion());
		assertEquals(false, set2.getDepletion());
	}

}
