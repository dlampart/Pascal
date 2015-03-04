package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.AnalyticVegas;
import ch.unil.genescore.vegas.GeneData;
import ch.unil.genescore.vegas.GeneDataFakePhenotype;
import ch.unil.genescore.vegas.GeneScoreEvaluator;
import ch.unil.genescore.vegas.GeneWithItsSnps;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulationFakeSignal;
import ch.unil.genescore.vegas.Snp;

public class EqtlScorerTopEqtlSnp extends EqtlScorer {
		
	public EqtlScorerTopEqtlSnp(){
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
		
		EqtlGenePair myEqtlGenePair = eqtlResults_.getEqtlGenePair(gene.symbol_);
		//myEqtlGenePair.overlapSnpsWithoutChecking(geneSnps);		
		myEqtlGenePair.getSnpWithTopEqtlValue(geneSnps);
		ArrayList<Snp> overlappedSnpList = myEqtlGenePair.getOverlappedSnpList();
		if(overlappedSnpList.size()==0){
			Main.print("gwas snps and eqtl-snps don't overlap (within the chosen gene-window). \n");
			return(null);
		}
		if(overlappedSnpList.size()>1){
			throw new RuntimeException(" within top-eqtl procedure: more than one SNP per gene");		
		}
		GeneData currentGeneData = null;
		if (Settings.useFakePhenotype_){			
			ReferencePopulationFakeSignal fakeSignalGenerator = new ReferencePopulationFakeSignal();
			fakeSignalGenerator.runFakeSignal();
			double[] fakeSignal = fakeSignalGenerator.getSignal();
			currentGeneData = new GeneDataFakePhenotype(overlappedSnpList, fakeSignal);
		}
		else{
			currentGeneData = new GeneData(overlappedSnpList);				
			}
			currentGeneData.processData();
			DenseMatrix corr = currentGeneData.getCorr();
			ArrayList<Double> snpScores =  currentGeneData.getScores();
			GeneScoreEvaluator evaluator = null;
												
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
}
