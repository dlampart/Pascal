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

import java.util.LinkedHashMap;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulation;


/**
 * Unit tests for Projection Pipeline.
 */

public class ProjectionVegasPipelineTest {
	
	static double delta;
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
		delta=1E-14;
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test PathwayMain.run() */
	@Test
	public void testRun() {
		testAnalyticRun(false);
		testProjectionRun(false);
	}
	

	@SuppressWarnings({ "unchecked", "rawtypes" })
	/** Test PathwayMain.run() */
	private void testAnalyticRun(boolean loadScores) {
		//major bummer!! i always used analytic vegas instead of projection approach !!!
		Settings.geneWindowUpstream_ = 110000;
		Settings.geneWindowDownstream_ = 40000;
		Settings.genomeAnnotation_="gencode";
		Settings.eigenValueFractionCut_ = 1e4;
		Settings.gencodeAnnotationFile_="resources/test/vegas/ProjectionVegasPipelineTest/gencode_only_gab4.gtf";
		Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
		Settings.withZScore_ = true;
		Settings.geneWiseSnpWeightsFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/Chr22CorrelatedStates.txt";
		//Settings.useProjectionVegas_ = true;
		Settings.useAnalyticVegas_ = true;
		Settings.useProjectionVegas_ = false;
		Settings.useMafCutoff_ = 0.0;
		Settings.fractionToBeExplained_ = 0.9999;
		Settings.chromosome_ = "chr22";	
		Settings.weightFileFormat_ =  "genewise"; 
		Settings.useImhof_ = false;
		Settings.useDavies_ = false;
		Settings.useFarebrother_ = true;
		Settings.checkOptions();
		
		GenomeWideScoring geneScorer = new GenomeWideScoring();
		
		LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();		
		ReferencePopulation myRefPop=new ReferencePopulation();
		
		geneScorer.loadGwasAndRelevantSnps(myRefPop);
		
		geneScorer.setGenes(genes.values());
		geneScorer.computeScores(false);
		
		double geneScore=genes.get("ENSG00000215568").getScore(0);
		assertEquals(geneScore,0.6168669005959458, 0.0001);
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	/** Test PathwayMain.run() */
	private void testProjectionRun(boolean loadScores) {
		
		Settings.geneWindowUpstream_ = 110000;
		Settings.geneWindowDownstream_ = 40000;
		Settings.genomeAnnotation_="gencode";
		Settings.gencodeAnnotationFile_="resources/test/vegas/ProjectionVegasPipelineTest/gencode_only_gab4.gtf";
		Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
		Settings.withZScore_ = true;
		Settings.geneWiseSnpWeightsFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/Chr22CorrelatedStates.txt";
		Settings.useProjectionVegas_ = true;
		Settings.useAnalyticVegas_ = true;
		Settings.useOnlyGwasSnps_ = false;
		Settings.chromosome_ = "chr22";	
		Settings.weightFileFormat_ =  "genewise"; 
		Settings.useImhof_ = false;
		Settings.useDavies_ = true;
		Settings.useFarebrother_ = false;
		Settings.checkOptions();
		
		GenomeWideScoring geneScorer = new GenomeWideScoring();
		
		LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();		
		ReferencePopulation myRefPop=new ReferencePopulation();
		
		geneScorer.loadGwasAndRelevantSnps(myRefPop);
		
		geneScorer.setGenes(genes.values());
		geneScorer.computeScores(false);
		
		double geneScore=genes.get("ENSG00000215568").getScore(0);
		assertEquals(geneScore,0.351490, 0.0001);
	}
}
