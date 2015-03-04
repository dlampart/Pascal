package ch.unil.genescore.vegas.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.GwasSnps;
import ch.unil.genescore.vegas.LdMat;
import ch.unil.genescore.vegas.LinkageDisequilibrium;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.SnpWeightMap;
import ch.unil.genescore.vegas.SnpWeightPairs;


public class LdMatTest {
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
		public void runTest(){
				
				
				SnpWeightMap myMap = new SnpWeightMap();
				String[] inputString = {"id1;3", "id2;1", "id3;4", "id4;4"};
				ArrayList<Snp> snpList = new ArrayList<Snp>();
				Snp snp1 = new Snp("id2");		
				
				snpList.add(snp1);
				byte[] fakeGenotype = {0,0,1,1,1,1,0,1};
				snp1.setGenotypes(fakeGenotype);
				 snp1.computeAlleleStats();
				 Snp.setGenotypeLength(8);
				 
				 
				Snp snp2 = new Snp("id5");
					
				snpList.add(snp2);
				byte[] fakeGenotype2 = {0,0,1,1,1,0,0,0};
				snp2.setGenotypes(fakeGenotype2);
				snp2.computeAlleleStats();
				
				Snp snp3 = new Snp("id3");
				
				snpList.add(snp3);
				byte[] fakeGenotype3 = {0,0,0,0,1, 1, 0,1};
				snp3.setGenotypes(fakeGenotype3);
				snp3.computeAlleleStats();
				
				Snp snp4 = new Snp("id4");
				
				snpList.add(snp4);
				byte[] fakeGenotype4 = {0,0,0,1,1,1, 0,0};
				snp4.setGenotypes(fakeGenotype4);
				snp4.computeAlleleStats();
				
				
				double[]  weights = {1,1,1,1};
				//}
				
				SnpWeightPairs pairs = new SnpWeightPairs(snpList,ConvenienceMethods.doubleArToArList(weights)); 
				LdMat ldMat = new LdMat( pairs,snpList);
				ldMat.processData();				
				//ldMat.setColumnWeights(weights);
				//ldMat.setRowWeights(weights);
				DenseMatrix myMat = ldMat.getWeighedCrossCorrelationMatrix();	
				assertTrue(myMat.get(1,0)==myMat.get(0,1));
				assertTrue(myMat.get(1,3)==myMat.get(3,1));
								
				double[]  weights2 = {1,2,1,1};
				SnpWeightPairs pairs2 = new SnpWeightPairs(snpList,ConvenienceMethods.doubleArToArList(weights2)); 				
				LdMat ldMat2 = new LdMat( pairs2,snpList);
				ldMat2.processData();				
				//ldMat2.setColumnWeights(weights2);
				//ldMat2.setRowWeights(weights2);
				DenseMatrix myMat2 = ldMat2.getWeighedCrossCorrelationMatrix();	
				
				assertTrue(myMat2.get(1,0)==myMat2.get(0,1)*2);
				assertTrue(myMat2.get(1,3)==myMat2.get(3,1)*2);
				//assertEquals(weightResults[1],4, 1e-10);
				
				}
			
		}

