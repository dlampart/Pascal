package ch.unil.genescore.prioritization;

import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.io.FilenameUtils;

import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.pathway.GeneScoreList;
import ch.unil.genescore.pathway.PathwayMain;
import ch.unil.genescore.vegas.DistributionMethods;
import ch.unil.genescore.vegas.GeneScoreEvaluator;


public class PrioritizationMain extends PathwayMain {

	
	PrunedConnectionTreeGetter net_ = null;	
	TreeMultimap<Double,Gene> processedGenes_ = null;
	String networkName_ = null;
	/** Constructor */
	public PrioritizationMain() {
		
		// Prefix used for output files
		//name_ = extractName(Settings.snpPvalFile_, Settings.geneSetFile_);
		name_ = extractName(Settings.snpPvalFile_, Settings.geneSetFile_);
		
		Main.println("LOADING INPUT FILES");
		Main.println("-------------------\n");		

		// Load gene sets and genes, optionally with scores
		genes_ = loadGenes();
	}
	//public PrioritizationMain(){
		
		//super();
	//}

	/** Run enrichment analysis */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void run() {

		
			Main.println("COMPUTING GENE SCORES");
			Main.println("---------------------\n");		
			
			
			// Sort genes by position 	
	//		if (!Settings.loadScoresFromFiles_){
				
		//		ArrayList<Gene> genes = new ArrayList<Gene>(genes_.values());
			//	Collections.sort(genes);
			// Compute scores for individual genes (even those that may later be merged, useful for output)
			//	geneScorer_.setGenes(genes);			
		//		geneScorer_.computeScores();
		//	}
			//else {
				Main.println("Loading genescores form files:");
				Main.println(Settings.geneScoreFile_);
				Main.println(Settings.metaGeneScoreFile_);
				genes_=loadScoresFromFile(Settings.geneScoreFile_);
				
		//	}
			setupNetwork(genes_.values());
			rankGeneScores(genes_.values());
			getNetworkName();
			runOverGenes();			
			enrichNetwork();
	}
	
	private void getNetworkName(){
		networkName_= FilenameUtils.getBaseName(Settings.netPath_);
	}
	
	private void enrichNetwork(){		
		processedGenes_ = net_.pruneTreeByDistance(processedGenes_);	
		rankGeneScores(processedGenes_.values());
		int n = processedGenes_.size();
		PrioritizationGeneResultsScore printer = new PrioritizationGeneResultsScore(); 
		printer.setExporter("overall");
		
		TreeMultimap<Double,Gene> tmpTree = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
		int count=0;
		for(Map.Entry<Double,Gene> entry : processedGenes_.entries()){
			count++;
			tmpTree.put(DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(count/(double)(n)), entry.getValue());
			if((count % 50)==0){
				AnalyticPrioritizationNoRemoval evaluator = new AnalyticPrioritizationNoRemoval(tmpTree);	
				evaluator.computeScore();
				 printer.writeLine(evaluator);				 
			}			
		}						
	}
	
	private String getOutputString(String localInfix){
		String inputData = null;
		if (Settings.loadScoresFromFiles_)
			inputData = FilenameUtils.getBaseName(Settings.geneScoreFile_);
		else
			inputData = Settings.gwasName_;
		String outStr = Settings.outputDirectory_ + "/" + inputData + "." + Settings.outputSuffix_ + "." + networkName_ + localInfix +".prioritization" + Settings.chromFileExtension_  + ".txt";
		return outStr;
	}	
	
	private void runOverGenes(){
		
		processedGenes_ = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
		PrioritizationGeneResultsScore printer = new PrioritizationGeneResultsScore(); 
		printer.setExporter("genescores");
		int count=0;
		for (Gene currentGene : genes_.values()){
			count++;
			System.out.println(count);
			TreeMultimap<Double,Gene> currentTree = net_.returnReverseOrderedConnectionTree(currentGene);	
			if (currentTree.size()>1){
				
				AnalyticPrioritization evaluator = new AnalyticPrioritization(currentTree);
				evaluator.computeScore();
				processedGenes_.put(1-evaluator.getScore()[0],currentGene);
				printer.writeLine(evaluator, currentGene);				
			}
		}
		
	}
	
	private void setupNetwork(Collection<Gene> genes){
		net_ = new PrunedConnectionTreeGetter();		
		net_.setExt((int) (Settings.mergeGenesDistance_*1E6));
		net_.setConnectionTreeGetter(Settings.getConnectionGetter());
		net_.setGenes(genes);
		net_.loadNetworkData();
		
	}
	
	/** Transform gene scores to ranks */
	private void rankGeneScores(Collection<Gene> genes) {
		
		// Transform scores to empirical probabilities and map to 1-df chi2 distribution
		//ChiSquaredDistribution chi2 = new ChiSquaredDistribution(1);
		double prevScore = -1;
		int count = 1;

		// Sort genes by score
		GeneScoreList rankedGenes = new GeneScoreList(genes, true);
		
		for (Gene g : rankedGenes.getGenes()) {
			// Check that genes are ordered
			double score[] = g.getScore();
			if (score[0] < prevScore)
				throw new RuntimeException("Genes are not ordered by score");
			prevScore = score[0];
			
			double empiricalProb = count / ((double) genes.size()+1);
			g.setNormalizedScore(empiricalProb);
			count++;
		}
	}
	
	

}

