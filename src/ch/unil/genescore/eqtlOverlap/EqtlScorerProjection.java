package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.Genome;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.AnalyticVegas;
import ch.unil.genescore.vegas.GeneData;
import ch.unil.genescore.vegas.GeneDataFakePhenotype;
import ch.unil.genescore.vegas.GeneScoreEvaluator;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.GeneDataInterface;

public class EqtlScorerProjection extends EqtlScorer {
	
	NeighboursCollection geneNeighbours_ = null;
	public EqtlScorerProjection(){
		super();
		geneNeighbours_	= new NeighboursCollection();
	}

	
	
	public HashSet<Snp> getEqtlInGenes(Collection<String> symbolList){
		HashSet<Snp> snpList = new HashSet<Snp>();
		EqtlGenePair myPair = null;
		for(String symbol : symbolList){
			if(eqtlResults_.eqtlGenePairExists(symbol)){
				//myPair =eqtlResults_.fillEqtlGenePair(symbol);
				eqtlResults_.fillEqtlGenePair(symbol, myEqtlGenePair_);
				myEqtlGenePair_.overlapSnps(refpop_.getSnpsWithGenotypes());
				snpList.addAll(myEqtlGenePair_.getOverlappedSnpList());
			}
		}	
		return(snpList);		
	}
	
	private HashSet<Snp> GetEqtlSnpsInNeighbourhoodGenes(Gene gene){
		ArrayList<String> symbolList = gene.getSymbolList();
		HashSet<String> neighbourSymbols  = geneNeighbours_.getAllNeighbourName(symbolList);
		HashSet<Snp> snpList = getEqtlInGenes(neighbourSymbols);		
		return snpList;
	}
	
	@Override
	protected GeneScoreEvaluator computeScore(Gene gene) {
	
		System.out.println(gene.symbol_);
		// Get the snps that are in the window around the given gene
		HashSet<Snp> geneSnps = refpop_.findSnpsHash(gene);
		refpop_.updateLoadedGenotypes(gene);
	
		removeLowMafSnps(geneSnps);
	
	
		
	
	Set<Snp> eqtlSnps = getEqtlInGenes(gene.getSymbolList());
	removeLowMafSnps(eqtlSnps);
	Set<Snp> eqtlNeighbourSnps = GetEqtlSnpsInNeighbourhoodGenes(gene);
	removeLowMafSnps(eqtlNeighbourSnps);
	
	geneSnps.addAll(eqtlSnps);
	geneSnps.removeAll(eqtlNeighbourSnps);
	
	if (geneSnps.size()==0){			
		Main.print("Gene has no SNPs  with maf larger than <mafcutoff>\n");
		return null;
	}
	if(!eqtlResults_.eqtlGenePairExists(gene.symbol_))
		Main.print("Gene has no eqtl associated\n");
	//	return(null);
		
	
	
		GeneDataInterface currentGeneData = null;
		currentGeneData = new GeneDataEqtlProjection(geneSnps,eqtlSnps,eqtlNeighbourSnps);
		currentGeneData.processData();
		DenseMatrix corr = currentGeneData.getCorr();
		ArrayList<Double> snpScores =  currentGeneData.getScores();
		GeneScoreEvaluator evaluator = null;
											
			//evaluator = new AnalyticVegas(snpScores, corr);
		evaluator = new AnalyticVegas(snpScores, corr);
			long t0 = System.currentTimeMillis();
			boolean success = evaluator.computeScore();
			System.out.println(success);
			long t1 = System.currentTimeMillis();
			
			// Set gene score
			gene.setScore(evaluator.getScore());
			gene.calcChi2StatFromScore();

			// Print info to console
			if (Settings.verbose_) {
				Main.print("\t" + Utils.padRight(Utils.chronometer(t1 - t0), 22));
				if (evaluator != null)
					Main.print(evaluator.getConsoleOutput());
				Main.println();
			}

			// If the score could not be computed (exception in analytic vegas)
			if (!success) {
				String str = gene.toString();
				str += "\t" + evaluator.getNoScoreOutput();
				noScores_.add(str);
				// Print error if not in verbose mode (otherwise it is already written above)
				if (!Settings.verbose_)
					Main.println(str);
			}

			return evaluator;
	}
	
	@Override
	public  void setGenes(Collection<Gene> genes) { genes_ = genes; 
		geneNeighbours_.findNeighbours(genes_);
	}
	
	
}
