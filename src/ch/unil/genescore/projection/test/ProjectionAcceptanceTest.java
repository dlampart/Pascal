package ch.unil.genescore.projection.test;


import static org.junit.Assert.*;
import org.junit.Ignore;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Set;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.projection.SnpWeightCreator;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.GwasSnps;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpWeightMap;

public class ProjectionAcceptanceTest {
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

		/**test only background file*/
	
	
		/**test only background file*/
		@Test			
		public void AcceptanceTest(){
			Settings.loadSettings();
			Settings.genomeAnnotation_="gencode";
			Settings.gencodeAnnotationFile_="resources/test/projection/projectionAcceptanceTest/"
					+ "gencode_only_pick1.gtf";
			Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
			Settings.withZScore_ = true;
			Settings.geneWiseSnpWeightsFile_ = "";
			Settings.weightFileFormat_ =  "genewise"; 
			Settings.useProjectionVegas_ = true;
			Settings.chromosome_ = "chr22";				
			Settings.removeCodingSnpsOfOtherGenes_=false;
			Settings.writeSnpBedFile_=false;
			Settings.writeTpedFile_=false;
			Settings.maxSnpsPerGene_=1000;
			Settings.geneWindowDownstream_=60000;
			Settings.geneWindowUpstream_=60000;
			Settings.useAnalyticVegas_=false;
			Settings.useProjectionVegas_=true;						
			Settings.bedFilePath_ = "";						
			//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
			//		+ "genicRegionPICK1.bed";
			//Settings.bedBackgroundFilePath_ = "resources/weightPreparationPipelines/1Kweights/genicRegions.bed";
			Settings.bedBackgroundExtension_ = 60000;			
			Settings.useOnlyGwasSnps_=true;
			Settings.useMafCutoffForProjection_ = 0.05;
			Settings.conditionFraction_ = 0.01;
			Settings.fractionToBeExplained_ = 0.9999;
			Settings.useImhof_ = true;
			
			LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();
			GenewiseSnpWeights finalStruct;
			//new GenewiseSnpWeights();
			GenomeWideScoring geneScore = new GenomeWideScoring();
			geneScore.setGenes(genes.values());
			ReferencePopulation myRefPop=new ReferencePopulation();				
			SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
			finalStruct = WeightInstance.getWeights();
			geneScore.loadGwasAndRelevantSnps(myRefPop);
			geneScore.computeScores(false);
			Collection<Gene> gene = geneScore.getGenes();
			Gene myGene = genes.get("ENSG00000100151");
			assertEquals(myGene.getScore(0),0.2068, 0.001);
			System.out.println("saf");
		}
		@Test			
		public void AcceptanceTest2(){
			Settings.loadSettings();
			Settings.genomeAnnotation_="gencode";
			Settings.gencodeAnnotationFile_="resources/test/projection/projectionAcceptanceTest/"
					+ "gencode_only_AC008132.13.gtf";
			Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
			Settings.withZScore_ = true;
			Settings.geneWiseSnpWeightsFile_ = "";
			Settings.weightFileFormat_ =  "genewise"; 
			Settings.useProjectionVegas_ = true;
			Settings.chromosome_ = "chr22";				
			Settings.removeCodingSnpsOfOtherGenes_=false;
			Settings.writeSnpBedFile_=false;
			Settings.writeTpedFile_=false;
			Settings.maxSnpsPerGene_=1000;
			Settings.geneWindowDownstream_=60000;
			Settings.geneWindowUpstream_=60000;
			Settings.useAnalyticVegas_=false;
			Settings.useProjectionVegas_=true;						
			Settings.bedFilePath_ = "";	
			
			//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
			//		+ "genicRegionAC008132.13.bed";
			//Settings.bedBackgroundFilePath_ = "resources/weightPreparationPipelines/1Kweights/genicRegions.bed";
			Settings.bedBackgroundExtension_ = 60000;			
			Settings.useOnlyGwasSnps_=true;
			Settings.useMafCutoffForProjection_ = 0.05;
			Settings.conditionFraction_ = 0.01;
			Settings.fractionToBeExplained_ = 0.9999;
			Settings.useImhof_ = true;
			Settings.useDavies_ = false;
			
		
			LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();
			GenewiseSnpWeights finalStruct;
			//new GenewiseSnpWeights();
			GenomeWideScoring geneScore = new GenomeWideScoring();
			geneScore.setGenes(genes.values());
			ReferencePopulation myRefPop=new ReferencePopulation();				
			SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
			finalStruct = WeightInstance.getWeights();
			geneScore.loadGwasAndRelevantSnps(myRefPop);
			geneScore.computeScores(false);
			Collection<Gene> gene = geneScore.getGenes();
			Gene myGene = genes.get("ENSG00000182356");
			assertEquals(myGene.getScore(0),0.870373, 0.01);
			System.out.println("saf");
		}
		@Test			
		public void AcceptanceTest3(){
			//should give same value as AcceptanceTest2 // 
			// 
			Settings.loadSettings();
			Settings.genomeAnnotation_="gencode";
			Settings.gencodeAnnotationFile_="resources/test/projection/projectionAcceptanceTest/"
					+ "gencode_only_AC008132.13.gtf";
			Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
			Settings.withZScore_ = true;
			Settings.geneWiseSnpWeightsFile_ = "";
			Settings.weightFileFormat_ =  "genewise"; 
			Settings.useProjectionVegas_ = true;
			Settings.chromosome_ = "chr22";				
			Settings.removeCodingSnpsOfOtherGenes_=false;
			Settings.writeSnpBedFile_=false;
			Settings.writeTpedFile_=false;
			Settings.maxSnpsPerGene_=1000;
			Settings.geneWindowDownstream_=60000;
			Settings.geneWindowUpstream_=60000;
			Settings.useAnalyticVegas_=false;
			Settings.useProjectionVegas_=true;		
			Settings.useImhof_ = true;
			Settings.useDavies_ = false;
			Settings.useFarebrother_ = false;
			
			Settings.bedFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
					+ "genicRegionAC008132.13_extended.bed";					
			//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
			//		+ "genicRegionAC008132.13.bed";
			//Settings.bedBackgroundFilePath_ = "resources/weightPreparationPipelines/1Kweights/genicRegions.bed";
			Settings.bedBackgroundExtension_ = 60000;	
			Settings.bedWeight_=1.0;
			Settings.useOnlyGwasSnps_=true;
			Settings.useMafCutoffForProjection_ = 0.05;
			Settings.conditionFraction_ = 0.01;
			Settings.fractionToBeExplained_ = 0.9999;
			
		
			LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();
			GenewiseSnpWeights finalStruct;
			//new GenewiseSnpWeights();
			GenomeWideScoring geneScore = new GenomeWideScoring();
			geneScore.setGenes(genes.values());
			ReferencePopulation myRefPop=new ReferencePopulation();				
			SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
			finalStruct = WeightInstance.getWeights();
			geneScore.loadGwasAndRelevantSnps(myRefPop);
			geneScore.computeScores(false);
			Collection<Gene> gene = geneScore.getGenes();
			Gene myGene = genes.get("ENSG00000182356");
			assertEquals(myGene.getScore(0),0.870373, 0.01);
			System.out.println("saf");
		}
	
		@Test			
		public void AcceptanceTest4(){
			//should give same value as AcceptanceTest2 // 
			// 
			Settings.loadSettings();
			Settings.genomeAnnotation_="gencode";
			Settings.gencodeAnnotationFile_="resources/test/projection/projectionAcceptanceTest/"
					+ "gencode_only_AC008132.13.gtf";
			Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
			Settings.withZScore_ = true;
			Settings.geneWiseSnpWeightsFile_ = "";
			Settings.weightFileFormat_ =  "genewise"; 
			Settings.useProjectionVegas_ = true;
			Settings.chromosome_ = "chr22";				
			Settings.removeCodingSnpsOfOtherGenes_=false;
			Settings.writeSnpBedFile_=false;
			Settings.writeTpedFile_=false;
			Settings.maxSnpsPerGene_=1000;
			Settings.geneWindowDownstream_=60000;
			Settings.geneWindowUpstream_=60000;
			Settings.useAnalyticVegas_=false;
			Settings.useProjectionVegas_=true;		
			Settings.useImhof_ = false;
			Settings.useDavies_ = false;
			Settings.useFarebrother_ = true;
			
			Settings.bedFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
					+ "genicRegionAC008132.13_extended.bed";					
			//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
			//		+ "genicRegionAC008132.13.bed";
			//Settings.bedBackgroundFilePath_ = "resources/weightPreparationPipelines/1Kweights/genicRegions.bed";
			Settings.bedBackgroundExtension_ = 60000;	
			Settings.bedWeight_=1.0;
			Settings.useOnlyGwasSnps_=true;
			Settings.useMafCutoffForProjection_ = 0.05;
			Settings.conditionFraction_ = 0.01;
			Settings.fractionToBeExplained_ = 0.9999;
			
		
			LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();
			GenewiseSnpWeights finalStruct;
			//new GenewiseSnpWeights();
			GenomeWideScoring geneScore = new GenomeWideScoring();
			geneScore.setGenes(genes.values());
			ReferencePopulation myRefPop=new ReferencePopulation();				
			SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
			finalStruct = WeightInstance.getWeights();
			geneScore.loadGwasAndRelevantSnps(myRefPop);
			geneScore.computeScores(false);
			Collection<Gene> gene = geneScore.getGenes();
			Gene myGene = genes.get("ENSG00000182356");
			assertEquals(myGene.getScore(0),0.870373, 0.01);
			System.out.println("saf");
		}
				
		@Test			
		public void AcceptanceTest5(){
			//should give same value as AcceptanceTest2 // 
			// 
			Settings.loadSettings();
			Settings.genomeAnnotation_="gencode";
			Settings.gencodeAnnotationFile_="resources/test/projection/projectionAcceptanceTest/"
					+ "gencode_only_AC008132.13.gtf";
			Settings.snpPvalFile_ = "resources/test/vegas/ProjectionVegasPipelineTest/HDL_ONE_EuropeansChr22.tbl";
			Settings.withZScore_ = true;
			Settings.geneWiseSnpWeightsFile_ = "";
			Settings.weightFileFormat_ =  "genewise"; 
			Settings.useProjectionVegas_ = true;
			Settings.chromosome_ = "chr22";				
			Settings.removeCodingSnpsOfOtherGenes_=false;
			Settings.writeSnpBedFile_=false;
			Settings.writeTpedFile_=false;
			Settings.maxSnpsPerGene_=1000;
			Settings.geneWindowDownstream_=60000;
			Settings.geneWindowUpstream_=60000;
			Settings.useAnalyticVegas_=false;
			Settings.useProjectionVegas_=true;		
			Settings.useImhof_ = false;
			Settings.useDavies_ = false;
			Settings.useFarebrother_ = false;
			
			Settings.bedFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
					+ "genicRegionAC008132.13_extended.bed";					
			//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionAcceptanceTest/"
			//		+ "genicRegionAC008132.13.bed";
			//Settings.bedBackgroundFilePath_ = "resources/weightPreparationPipelines/1Kweights/genicRegions.bed";
			Settings.bedBackgroundExtension_ = 60000;	
			Settings.bedWeight_=1.0;
			Settings.useOnlyGwasSnps_=true;
			Settings.useMafCutoffForProjection_ = 0.05;
			Settings.conditionFraction_ = 0.01;
			Settings.fractionToBeExplained_ = 0.9999;
			
		
			LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation();
			GenewiseSnpWeights finalStruct;
			//new GenewiseSnpWeights();
			GenomeWideScoring geneScore = new GenomeWideScoring();
			geneScore.setGenes(genes.values());
			ReferencePopulation myRefPop=new ReferencePopulation();				
			SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
			finalStruct = WeightInstance.getWeights();
			geneScore.loadGwasAndRelevantSnps(myRefPop);
			geneScore.computeScores(false);
			Collection<Gene> gene = geneScore.getGenes();
			Gene myGene = genes.get("ENSG00000182356");
			assertEquals(myGene.getScore(0),0.870373, 0.01);
			System.out.println("saf");
		}
		
}

