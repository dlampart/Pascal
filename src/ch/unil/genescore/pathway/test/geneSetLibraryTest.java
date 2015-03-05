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
package ch.unil.genescore.pathway.test;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Settings;

import ch.unil.genescore.pathway.GeneSet;
import ch.unil.genescore.pathway.GeneSetLibrary;
import ch.unil.genescore.pathway.MetaGene;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;


public class geneSetLibraryTest extends GeneSetLibrary {

	
	static GeneSet  geneSet_= new GeneSet("test");
	//static GeneSetLibrary geneSetLib = new GeneSetLibrary();
	//setGenesForSimulation(ArrayList<Gene> genesForSimulation)
	//Ar
	
//	public static void main(String[] args){
//		geneSetLibraryDebug main = new geneSetLibraryDebug();
//		main.run2();
		
	//}
	
	@BeforeClass
	public static void testSetup() {
	//	Settings.loadSettings();	
	}

	@AfterClass
	public static void testCleanup() { }
	
	public class GeneSetLibraryForTest extends GeneSetLibrary{
		
		public double[] wrapHypGeomPvalue2(double[] set, double[] totSet){
			double[] out = computeHypGeomPvalue2(set, totSet);
			return out;
		}
	} 
	
	
	@Test
	public void calculateSamplingWeightsTest(){
		
		Gene g1 = new Gene("id1");g1.chr_="chr1";g1.start_=100;g1.end_=200;
		Gene g2 = new Gene("id2");g2.chr_="chr1";g2.start_=200;g2.end_=300;
		Gene g3 = new Gene("id3");g3.chr_="chr1";g3.start_=300;g3.end_=400;
		Gene g4 = new Gene("id4");g4.chr_="chr1";g4.start_=400;g4.end_=500;
		Gene g5 = new Gene("id5");g5.chr_="chr1";g5.start_=400;g5.end_=500;
		Gene g6 = new Gene("id6");g6.chr_="chr1";g6.start_=400;g6.end_=500;
		Gene g7 = new Gene("id7");g7.chr_="chr1";g7.start_=400;g7.end_=500;
		Gene g8 = new Gene("id8");g8.chr_="chr1";g8.start_=400;g8.end_=500;
		Gene g9 = new Gene("id9");g9.chr_="chr1";g9.start_=600;g9.end_=800;
		Gene g10 = new Gene("id10");g10.chr_="chr1";g10.start_=600;g10.end_=800;
		
		
		TreeSet<Gene> tg1 = new TreeSet<Gene>();
		tg1.add(g1);
		tg1.add(g2);
		//tg1.add(g4);
		
		TreeSet<Gene> tg2 = new TreeSet<Gene>();
		tg2.add(g2);
		tg2.add(g3);	
		TreeSet<Gene> tg3 = new TreeSet<Gene>();
		tg3.add(g8);
		tg3.add(g9);	
		
		MetaGene mg1 = new MetaGene(tg1);
		MetaGene mg2 = new MetaGene(tg2);
		MetaGene mg3 = new MetaGene(tg3);
		
		ArrayList<Gene> gl1 = new ArrayList<Gene>();
		gl1.add(g1);gl1.add(g3);gl1.add(g4);gl1.add(g5);
		
		
		ArrayList<Gene> gl2 = new ArrayList<Gene>();
		gl2.add(mg1);gl2.add(g3);gl2.add(g4);gl2.add(g6);
		
		ArrayList<Gene> gl3 = new ArrayList<Gene>();
		gl3.add(mg2);gl3.add(g6);gl3.add(g7);gl3.add(mg3);gl3.add(g10);
		
		
		GeneSet set1 = new GeneSet("pw1",gl1);//g1,g3,g4,g5
		GeneSet set2 = new GeneSet("pw2",gl2);//mg1:g1g2,g4,g3,g6
		GeneSet set3 = new GeneSet("pw3",gl3);//mg2:g2g3,g6,g7,mg3:g8g9:g10
		
		set1.setMetaGenes(new HashSet<MetaGene>());
		set2.setMetaGenes(new HashSet<MetaGene>());
		set3.setMetaGenes(new HashSet<MetaGene>());
		
		set2.getMetaGenes().add(mg1);
		set3.getMetaGenes().add(mg2);
		set3.getMetaGenes().add(mg3);
				
		geneSets_ = new ArrayList<GeneSet>();
		geneSets_.add(set1);
		geneSets_.add(set2);
		geneSets_.add(set3);
		genes_ = new HashSet<Gene>();
		genes_.add(g1);genes_.add(g2);genes_.add(g3);genes_.add(g4);genes_.add(g5);genes_.add(g6);genes_.add(g7);genes_.add(g8);genes_.add(g9);genes_.add(g10);
		metaGenes_ = new HashMap<String,MetaGene>();
		metaGenes_.put(mg1.id_,mg1);metaGenes_.put(mg2.id_,mg2);metaGenes_.put(mg3.id_,mg3);
		calculateSamplingWeights();
		
		//testing samplingWeightHelper_:
		assertTrue(getSamplingWeightHelper().get("id1").contains(set1));
		assertTrue(getSamplingWeightHelper().get("id1").contains(set2));
		assertTrue(!getSamplingWeightHelper().get("id1").contains(set3));
		assertTrue(!getSamplingWeightHelper().get("id2").contains(set1));
		assertTrue(getSamplingWeightHelper().get("id2").contains(set2));
		assertTrue(getSamplingWeightHelper().get("id2").contains(set3));
		assertTrue(!getSamplingWeightHelper().get("id6").contains(set1));
		assertTrue(getSamplingWeightHelper().get("id6").contains(set2));
		assertTrue(getSamplingWeightHelper().get("id6").contains(set3));
		assertTrue(getSamplingWeightHelper().get("id8").contains(set3));
		
		//testting weights of genes:		
		assertEquals(g1.getSamplingWeight(),(1.0/16+1.0/16), 0.0001);
		assertEquals(g2.getSamplingWeight(),(1.0/16+1.0/25),0.0001);
		assertEquals(g9.getSamplingWeight(),(1.0/25),0.0001);
		
		//testing weights of metagenes:
		assertEquals(mg1.getSamplingWeight(),(1.0/16+1.0/16+1.0/25), 0.0001);
		assertEquals(mg2.getSamplingWeight(),(1.0/16+1.0/16+1.0/25), 0.0001);
		assertEquals(mg3.getSamplingWeight(),(1.0/25), 0.0001);
		
		computeApproxPathwayCorrelation();
		
		assertEquals(pathwayCorMat_.get(0, 1),0.75,0.0001);
		assertEquals(pathwayCorMat_.get(1, 2),0.6708,0.001);
		assertEquals(pathwayCorMat_.get(2, 1),0.4472,0.0001);		
				
	}
	
	@Test
	public void hypgeomPvalue2Test(){
		double[] totSet = new double[50];
		for (int i=0; i<50;i++){
			totSet[i]=i;
		}
		double[] subSet = {30,40,41,42,43,44,45,46,48,49};
		Settings.hypGeomQuantiles_[0] =  0.89;
		GeneSetLibraryForTest myLib = new GeneSetLibraryForTest();
		double[] out = myLib.wrapHypGeomPvalue2(subSet,totSet);		
		double res = 0.003964583+0.0001189375;
		assertEquals(out[0],res,0.000001);	
	}
	
	@Test
	public void simulatePvalueWeightedSamplingTest(){
		
		double[] allSamplingWeights = {1,1,2,2};
		double[] allChiSqVals = {3,4,8,10};
		double[] setChiSqVals = {3,4.5};
		double p=0;
		double p1=0;
		for (int i = 0; i<10;i++){
			p1 = simulatePvalueWeightedSampling(setChiSqVals, allChiSqVals, allSamplingWeights);
			p= p+p1;
		}
		p=p/10;				
		assertEquals(p,0.933333, 0.003); //1-1/6*1/5-1/6*1/5
		System.out.println("sdf");
		p=0;
		p1=0;
		setChiSqVals[0]=17;
		setChiSqVals[1]=0;
		for (int i = 0; i<10;i++){
			p1 = simulatePvalueWeightedSampling(setChiSqVals, allChiSqVals, allSamplingWeights);
			p= p+p1;
		}
		p=p/10;				
		assertEquals(p,0.3333, 0.003);		
		System.out.println("sdf");
		
		p=0;
		p1=0;
		setChiSqVals[0]=13.5;
		setChiSqVals[1]=0;
		for (int i = 0; i<10;i++){
			p1 = simulatePvalueWeightedSampling(setChiSqVals, allChiSqVals, allSamplingWeights);
			p= p+p1;
		}
		p=p/10;				
		assertEquals(p,0.48333, 0.003);	//:: 1/6+1/6+2/30+2/24	
		System.out.println("sdf");
	}
	
}
