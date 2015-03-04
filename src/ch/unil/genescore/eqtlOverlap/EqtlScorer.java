package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;

import com.sun.xml.internal.ws.policy.privateutil.PolicyUtils.Collections;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GeneData;
import ch.unil.genescore.vegas.GeneDataEqtlZscoreHolder;
import ch.unil.genescore.vegas.GeneDataFakePhenotype;
import ch.unil.genescore.vegas.GeneScoreEvaluator;
import ch.unil.genescore.vegas.GeneWithItsSnps;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulationFakeSignal;
import ch.unil.genescore.vegas.Snp;


public class EqtlScorer extends GenomeWideScoring {
	EqtlResults eqtlResults_ = null;
	protected EqtlGenePair myEqtlGenePair_ = null; 
	public EqtlScorer(){
		super();		
	}
		
		
	@Override
	protected GeneScoreEvaluator computeScore(Gene gene) {
		
		System.out.println(gene.symbol_);
		// Get the snps that are in the window around the given gene
		ArrayList<Snp> geneSnps = refpop_.findSnps(gene);
		GeneWithItsSnps geneAndSnps = new GeneWithItsSnps(gene,geneSnps);		
		
		
		
		refpop_.updateLoadedGenotypes(gene);
		
		removeLowMafSnps(geneSnps);
		if (geneSnps.size()==0){			
			Main.print("Gene has no SNPs  with maf larger than <mafcutoff>\n");
			return null;
		}
		
		if(!eqtlResults_.eqtlGenePairExists(gene.symbol_)){
			Main.print("Gene has no eqtl associated\n");
			return(null);
			
		}
		
		//EqtlGenePair myEqtlGenePair = eqtlResults_.getEqtlGenePair(gene.symbol_);
		eqtlResults_.fillEqtlGenePair(gene.symbol_,myEqtlGenePair_);
		myEqtlGenePair_.overlapSnps(geneSnps);
		double[] eqtlZscores = myEqtlGenePair_.getOverlappedEqtlZscores();
		ArrayList<Snp> overlappedSnpList = myEqtlGenePair_.getOverlappedSnpList();
		if(overlappedSnpList.size()==0){
			Main.print("gwas snps and eqtl-snps don't overlap (within the chosen gene-window). \n");
			return(null);
		}
		GeneDataEqtlZscoreHolder currentGeneData = null;
		if (Settings.useFakePhenotype_){			
			ReferencePopulationFakeSignal fakeSignalGenerator = new ReferencePopulationFakeSignal();
			fakeSignalGenerator.runFakeSignal();
			double[] fakeSignal = fakeSignalGenerator.getSignal();
			currentGeneData = new GeneDataFakePhenotypeWithEqtl(geneSnps, fakeSignal, eqtlZscores);	
	}else
				currentGeneData = new GeneDataWithEqtl(overlappedSnpList, eqtlZscores);
						
		currentGeneData.processData();				
		GeneScoreEvaluator myEvaluator = new EqtlEvaluator(currentGeneData);
		myEvaluator.computeScore();
		gene.setScore(myEvaluator.getScore());
		gene.calcChi2StatFromScore();
		return myEvaluator;				
	}	
	

	public void setEqtlResults(EqtlResults eqtlResults){
		eqtlResults_=eqtlResults;
	}
	public void setEqtlGenePair(EqtlGenePair pair){
		myEqtlGenePair_=pair;
	}		
	
}

