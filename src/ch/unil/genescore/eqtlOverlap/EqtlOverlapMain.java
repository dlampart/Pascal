package ch.unil.genescore.eqtlOverlap;

import java.util.LinkedHashMap;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulation;


/**
 * Main class handling eqtl-Overlap analysis 
 * (loads data, computes overlap, writes output).
 */
public class EqtlOverlapMain  {
	
	String name_ =null;
	 LinkedHashMap<String, Gene> genes_ =null;
	 EqtlScorer geneScorer_ = null;
		// ============================================================================
		// PUBLIC METHODS
	
		public void fillEqtlScorer(){
			
			name_ = extractName(Settings.snpPvalFile_, Settings.eqtlFile_);
			
			Main.println("LOADING INPUT FILES");
			Main.println("-------------------\n");		

			// Load gene sets and genes, optionally with scores
			genes_ = loadGenes();
			if(Settings.runEqtlProjection_){
				geneScorer_ = new  EqtlScorerProjection();
				
			}
			else{geneScorer_ = new EqtlScorer();}
			 
			if (Settings.onlyTopOverlappedEqtl_){
				geneScorer_.setEqtlGenePair(new  TopSnpEqtlGenePair());
			}			
			else {			
				geneScorer_.setEqtlGenePair(new EqtlGenePairDirectionIgnored());
		//	if (!Settings.withZScore_){
		//			throw new RuntimeException("Eqtl-analysis needs z-scores in gwas-file to work");
		//		}									
			}
			geneScorer_.setGenes(genes_.values());	
			ReferencePopulation myRefPop=new ReferencePopulation();
			myRefPop.loadGwasAndRelevantSnps();		
			geneScorer_.setReferencePopulation(myRefPop);
			//SnpWeightCreator WeightInstance = new SnpWeightCreator(myRefPop,genes_);
			//finalStruct = WeightInstance.getWeights();		
			
		}
	public void fillEqtlResults(){
		EqtlResults myEqtlResults = new EqtlResults();
		myEqtlResults.loadEqtlFile(Settings.eqtlFile_);
		geneScorer_.setEqtlResults(myEqtlResults);
		
	}	
	
	public void run(){
		fillEqtlScorer();
		fillEqtlResults();
		geneScorer_.setAdditionalOutputFileSuffix(".eqtl");
		geneScorer_.computeScores();
	}

	
	
 private LinkedHashMap<String, Gene> loadGenes() {		
		// The genes from the genome annotation
		LinkedHashMap<String, Gene> genes;		
			// Load genome annotation (genes without scores)
			GeneAnnotation annot = GeneAnnotation.createAnnotationInstance();
			genes = annot.loadAnnotation();
			// Remove excluded genes (this updates genes_)
			if (!Settings.excludedGenesFile_.equals("")) {
				int numRemoved = annot.removeGenes(Settings.excludedGenesFile_);
				Main.println("- " + numRemoved + " genes removed because they are in the excluded genes file");
			}
		return genes;
	}
	
	/** Extract the name for this run from the snp and eqtl annotation files */
	private String extractName(String geneScoreFile, String eqtlFile) {
		
		String gwasName = Utils.extractBasicFilename(geneScoreFile, false);
		String functName = Utils.extractBasicFilename(eqtlFile, false);
		
		functName = functName.replace("_nodeProperties", "");
		functName = functName.replace("_undir", "");
		
		return gwasName + "--" + functName;
	}
	public EqtlScorer getGeneScorer(){
		return(geneScorer_);
	}
}
