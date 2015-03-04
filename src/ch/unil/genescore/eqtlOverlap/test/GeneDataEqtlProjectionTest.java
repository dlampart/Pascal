package ch.unil.genescore.eqtlOverlap.test;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.*;
import no.uib.cipr.matrix.DenseMatrix;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.eqtlOverlap.EqtlEvaluator;
import ch.unil.genescore.eqtlOverlap.EqtlGenePair;
import ch.unil.genescore.eqtlOverlap.GeneDataEqtlProjection;
import ch.unil.genescore.eqtlOverlap.GeneDataWithEqtl;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.Snp;


public class GeneDataEqtlProjectionTest {


	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
	}

	@AfterClass
	public static void testCleanup() { }
	
	private byte[] fillGenotype(int genotypeLength,int freqLength){		
		byte[]geno = new byte[genotypeLength];
		int looplength;
		int count=0;
		looplength=genotypeLength/freqLength;
		for(int j=0; j<looplength; j++){
			for(int i=0; i<freqLength; i++){
				if(i<freqLength/2){
					geno[count]=1;
				}
				count++;
			}	
		}
		return(geno);
	}	
	
	@Test
	public void overlapSnps_test(){
		
	
		
		int freqLength;
	
		ArrayList<Snp> eqtlSnpList = new ArrayList<Snp>();
		Snp snp = null;
		//int i;
		byte[]geno1 = fillGenotype(128,2);
		byte[]geno2 = fillGenotype(128,4);
		byte[]geno3 = fillGenotype(128,8);
		byte[]geno4 = fillGenotype(128,16);
		
		int genotypeLength=128;
		//freqLength=2;
		
		//for(int j=0; j<genotype)
		Snp snp1 = new Snp("id1",1,1.5);snp1.setGenotypes(geno1);snp1.setGenotypeLength(128);snp1.computeAlleleStats();
		Snp snp2 = new Snp("id2",1,3);snp2.setGenotypes(geno2);	snp2.computeAlleleStats();
		Snp snp3 = new Snp("id3",1,3);snp3.setGenotypes(geno3);	snp3.computeAlleleStats();
		Snp snp4 = new Snp("id4",1,3);snp4.setGenotypes(geno1);	snp4.computeAlleleStats();
		HashSet<Snp> gSet = new HashSet<Snp>();
		gSet.add(snp1);gSet.add(snp2);gSet.add(snp3);
		HashSet<Snp> eSet = new HashSet<Snp>();eSet.add(snp4);
		HashSet<Snp> enSet = new HashSet<Snp>();
		GeneDataEqtlProjection pr = new GeneDataEqtlProjection(gSet,  eSet,enSet);
		pr.processData();
		DenseMatrix outM = pr.getCorr();
		ArrayList<Double> outZ = pr.getScores();
		System.out.println("sfds");
		
		gSet = new HashSet<Snp>();
		gSet.add(snp1);gSet.add(snp2);gSet.add(snp3);
		eSet = new HashSet<Snp>();
		enSet = new HashSet<Snp>();enSet.add(snp4);
		pr = new GeneDataEqtlProjection(gSet,eSet,enSet);
		pr.processData();
		outM = pr.getCorr();
		outZ = pr.getScores();
		System.out.println("sfds");
		
		Snp snp5 = new Snp("id5",1,3);snp5.setGenotypes(geno1);	snp5.computeAlleleStats();
		gSet = new HashSet<Snp>();
		gSet.add(snp1);gSet.add(snp2);gSet.add(snp3);
		eSet = new HashSet<Snp>();eSet.add(snp4);
		enSet = new HashSet<Snp>();enSet.add(snp5);
		pr = new GeneDataEqtlProjection(gSet,eSet,enSet);
		pr.processData();
		outM = pr.getCorr();
		outZ = pr.getScores();
		System.out.println("sfds");
		
		
	}
		
			
		
	
				
		
		
		
}
