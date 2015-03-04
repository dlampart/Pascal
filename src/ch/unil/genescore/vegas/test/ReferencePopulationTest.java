package ch.unil.genescore.vegas.test;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.GwasSnps;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import ch.unil.genescore.vegas.ProjectionVegas;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpPositionStream;



// ============================================================================
// SETUP

@BeforeClass
public static void testSetup() {
	Settings.loadSettings();
}

@AfterClass
public static void testCleanup() { }


/**
 * Unit tests for ReferencePopulationTest
 *TODO: only partially covered
 */
public class ReferencePopulationTest {
	
		
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
		Settings.genomeAnnotation_="gencode";
		Settings.gencodeAnnotationFile_="resources/test/vegas/ProjectionVegasPipelineTest/gencode_only_gab4.gtf";
		Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
		Settings.withZScore_ = true;
		Settings.geneWiseSnpWeightsFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/Chr22CorrelatedStates.txt";
		Settings.weightFileFormat_ =  "genewise"; 
		Settings.useProjectionVegas_ = true;
		Settings.chromosome_ = "chr22";	
		
		Settings.removeCodingSnpsOfOtherGenes_=false;
		Settings.writeSnpBedFile_=false;
		Settings.writeTpedFile_=false;
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS
	
	/** Test test statistic */
	@Test
	public void initializeAllRelevantSnps(){
		//initialize gene
		LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();
		Gene myGene = genes.get("ENSG00000215568");
				
		
		
		ReferencePopulation myRefPop = new ReferencePopulation();
		myRefPop.initialize(Settings.chromosome_);
		
		//set up snps with weights
		GenewiseSnpWeights snpWeightsForEachGene = new GenewiseSnpWeights();
		snpWeightsForEachGene.loadSnpWeightsForEachGeneFromFile(null);
		Set<String> allsnpIdsLoadedFromFile = snpWeightsForEachGene.returnAllSnpIds();
		myRefPop.addToRelevantSnps(allsnpIdsLoadedFromFile);
		//set up snps with gwas-data		
		GwasSnps snps = new GwasSnps();
		snps.setHeader(true);
		snps.loadSnpPvalZval(Settings.snpPvalFile_);
		myRefPop.setGwasSnps(snps.getSnpsInList());			
		// read-position and initialize.		
		myRefPop.initializeSnps();
		myRefPop.getSnpsWithGenotypes();
		myRefPop.updateLoadedGenotypes(myGene);
		LinkedList<Snp> genotypeSnps = myRefPop.getSnpsWithGenotypes();
		HashSet<String> mySnpList = new HashSet<String>();
		for (Snp mySnp : genotypeSnps){
			mySnpList.add(mySnp.id_);
		}
		// relevant snps
		assertTrue(mySnpList.contains("rs28502153"));// should be in
		assertTrue(mySnpList.contains("rs2110439"));// should be in
		assertFalse(mySnpList.contains("rs189246845"));//contains no minor alleles
		assertFalse(mySnpList.contains("rs116250018"));//contains no minor alleles
		
		// check GWAS snps
		assertTrue(mySnpList.contains("rs5994110"));// should be in
		assertTrue(mySnpList.contains("rs2041607"));// should be in		
		System.out.println("bljbijsf");
	}
	@Test
	public void AllPosTest(){
		ReferencePopulation myRefPop = new ReferencePopulation();
		myRefPop.initialize(Settings.chromosome_);
		 SnpPositionStream posStream = myRefPop.getPositionStream(Settings.chromosome_);
		 while(posStream.streamOpen()){
			 OverlappedGenomicElement el = posStream.getNextAsOverlappedGenomicElement();
			 //System.out.println(el.getMainElement().id_);
			 //System.out.println(el.getMainElement().start_);
			 if(el.getMainElement().id_.equals(".")){
				 System.out.println(el.getMainElement().start_);
			 }
		 }
	}
}


