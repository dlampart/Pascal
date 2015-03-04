package ch.unil.genescore.projection;


import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.genescore.vegas.SnpWeightMap;

public class SnpWeightCreatorTest {
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
		@Test			
		public void SnpWeightCreatorTest1(){
			
			
			Settings.genesToBeLoadedFile_ = "resources/test/projection/projectionPipelineTest/CCT8L2_ENSG";
			//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionPipelineTest/bedFileOneGene.bed";						
			Settings.bedFilePath_ = "";
			Settings.bedBackgroundExtension_=50000;
			Settings.chromosome_="chr22";
			LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation(Settings.genesToBeLoadedFile_);						
			ReferencePopulation myRefPop=new ReferencePopulation();				
			
			SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
			GenewiseSnpWeights weights = WeightInstance.getWeights();
			SnpWeightMap myMap = weights.getSnpWeightMap("ENSG00000198445");
			ArrayList<String> rs_nrs = myMap.getSnpIds();
			assertTrue(rs_nrs.contains("rs139007749"));
			assertTrue(rs_nrs.contains("rs9680652"));			
			assertTrue(!rs_nrs.contains("rs191686510"));//no minor allele occurrence
			assertTrue(!rs_nrs.contains("rs184029972"));//no minor allele occurrence
			assertTrue(rs_nrs.contains("rs4008560"));
			assertTrue(rs_nrs.contains("rs4008561"));
			assertTrue(rs_nrs.contains("rs7289232"));
			assertTrue(rs_nrs.contains("rs138591121"));
			assertTrue(!rs_nrs.contains("rs189897465"));//no minor allele occurrence									
		}
		
		// TESTS

				/**test only background and bed-file*/
				@Test			
				public void SnpWeightCreatorTest2(){
					
					
					Settings.genesToBeLoadedFile_ = "resources/test/projection/projectionPipelineTest/CCT8L2_ENSG";
					//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionPipelineTest/bedFileOneGene.bed";						
					Settings.bedFilePath_ = "resources/test/projection/projectionPipelineTest/ElementFile.bed";						
					Settings.chromosome_="chr22";
					
					LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation(Settings.genesToBeLoadedFile_);						
					ReferencePopulation myRefPop=new ReferencePopulation();				
					
					SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
					GenewiseSnpWeights weights = WeightInstance.getWeights();
					SnpWeightMap myMap = weights.getSnpWeightMap("ENSG00000198445");										
					//all overlaps found
					assertTrue(myMap.getWeight("rs9680694")==10);
					assertTrue(myMap.getWeight("rs149763016")==10);
					assertTrue(myMap.getWeight("rs5747812")==10);					
					assertTrue(myMap.getWeight("rs12172470")==10);
					assertTrue(myMap.getWeight("rs116313276")==10);
					assertTrue(myMap.getWeight("rs146206935")==10);
					assertTrue(myMap.getWeight("rs181698943")==10);
					assertTrue(myMap.getWeight("rs145889691")==10);
										
					assertTrue(!myMap.checkExistence("rs191693315")); //overlaps but no minor allele
					assertTrue(!myMap.checkExistence("rs146793829")); //overlaps but no minor allele
					assertTrue(!myMap.checkExistence("rs186564474")); //overlaps but no minor allele
					assertTrue(!myMap.checkExistence("rs186564474")); //overlaps but no minor allele					
					assertTrue(myMap.getWeight("rs139007749")==1);					
					ArrayList<String> rs_nrs = myMap.getSnpIds();
					int count=0;
					for (String rs : rs_nrs){
						if (myMap.getWeight(rs)==10){
							System.out.println(rs);
							count++;
						}
					}
					assertTrue(count==8);
					
				}
				

				/**test only background and bed-file*/
				//TODO: this takes a long time to run put into AllTest. in different suite such for long time running.
				@Test			
				public void SnpWeightCreatorTest3(){
					
					
					Settings.genesToBeLoadedFile_ = "resources/test/projection/projectionPipelineTest/CCT8L2_ENSG";
					//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionPipelineTest/bedFileOneGene.bed";						
					Settings.bedFilePath_ = "resources/test/projection/projectionPipelineTest/ElementFile.bed";						
					Settings.chromosome_="";
					
					LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation(Settings.genesToBeLoadedFile_);						
					ReferencePopulation myRefPop=new ReferencePopulation();				
					
					SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
					GenewiseSnpWeights weights = WeightInstance.getWeights();
					SnpWeightMap myMap = weights.getSnpWeightMap("ENSG00000198445");										
					//all overlaps found
					assertTrue(myMap.getWeight("rs9680694")==10);
					assertTrue(myMap.getWeight("rs149763016")==10);
					assertTrue(myMap.getWeight("rs5747812")==10);					
					assertTrue(myMap.getWeight("rs12172470")==10);
					assertTrue(myMap.getWeight("rs116313276")==10);
					assertTrue(myMap.getWeight("rs146206935")==10);
					assertTrue(myMap.getWeight("rs181698943")==10);
					assertTrue(myMap.getWeight("rs145889691")==10);
										
					assertTrue(!myMap.checkExistence("rs191693315")); //overlaps but no minor allele
					assertTrue(!myMap.checkExistence("rs146793829")); //overlaps but no minor allele
					assertTrue(!myMap.checkExistence("rs186564474")); //overlaps but no minor allele
					assertTrue(!myMap.checkExistence("rs186564474")); //overlaps but no minor allele					
					assertTrue(myMap.getWeight("rs139007749")==1);					
					ArrayList<String> rs_nrs = myMap.getSnpIds();
					int count=0;
					for (String rs : rs_nrs){
						if (myMap.getWeight(rs)==10){
							System.out.println(rs);
							count++;
						}
					}
					assertTrue(count==8);
					
				}
				/**test background and bed-file and also test filterOnBed-option*/
				@Test			
				public void SnpWeightCreatorTest4(){
					
					
					Settings.genesToBeLoadedFile_ = "resources/test/projection/projectionPipelineTest/CCT8L2_ENSG";
					//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionPipelineTest/bedFileOneGene.bed";						
					Settings.bedFilePath_ = "resources/test/projection/projectionPipelineTest/ElementFileWithId.bed";						
					Settings.chromosome_="chr22";
					Settings.filterOnBed_= true;
					
					LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation(Settings.genesToBeLoadedFile_);						
					ReferencePopulation myRefPop=new ReferencePopulation();				
					
					SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
					GenewiseSnpWeights weights = WeightInstance.getWeights();
					SnpWeightMap myMap = weights.getSnpWeightMap("ENSG00000198445");										
					//all overlaps found
					System.out.println(myMap.getWeight("rs9680694"));
					System.out.println(myMap.getWeight("rs149763016"));
					System.out.println(myMap.getWeight("rs5747812"));
					assertTrue(myMap.getWeight("rs9680694")!=10);//is not in
					assertTrue(myMap.getWeight("rs149763016")!=10);//is notin
					assertTrue(myMap.getWeight("rs5747812")==10);//is in	
					assertTrue(myMap.getWeight("rs12172470")!=10); //is in
					//
					System.out.println(myMap.getWeight("rs12172470"));
					
				
				}
				/**test background and bed-file and also test filterOnBed-option*/
				@Test			
				public void SnpWeightCreatorTest5(){
					
					
					Settings.genesToBeLoadedFile_ = "resources/test/projection/projectionPipelineTest/CCT8L2_ENSG";
					//Settings.bedBackgroundFilePath_ = "resources/test/projection/projectionPipelineTest/bedFileOneGene.bed";						
					Settings.bedFilePath_ = "resources/test/projection/projectionPipelineTest/ElementFileWithIdDouble.bed";						
					Settings.chromosome_="chr22";
					Settings.filterOnBed_= true;
					
					LinkedHashMap<String, Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation(Settings.genesToBeLoadedFile_);						
					ReferencePopulation myRefPop=new ReferencePopulation();				
					
					SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes);
					GenewiseSnpWeights weights = WeightInstance.getWeights();
					SnpWeightMap myMap = weights.getSnpWeightMap("ENSG00000198445");										
					//all overlaps found
					System.out.println(myMap.getWeight("rs9680694"));
					System.out.println(myMap.getWeight("rs149763016"));
					System.out.println(myMap.getWeight("rs5747812"));
					assertTrue(myMap.getWeight("rs9680694")!=10);//is not in
					assertTrue(myMap.getWeight("rs149763016")!=10);//is notin
					assertTrue(myMap.getWeight("rs5747812")==10);//is in	
					assertTrue(myMap.getWeight("rs12172470")!=10); //is in
					//
					System.out.println(myMap.getWeight("rs12172470"));
					
				
				}


}
