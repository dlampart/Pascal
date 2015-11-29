/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.vegas;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;

import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.main.Pascal;


/**
 * Summarize p-values of multiple snps at a locus
 */
public class GenomeWideScoring {


	/** Appended to output files right before the file extension */
	private String additionalOutputFileSuffix_ = "";
	
	/** The loaded genes (all genes if genesToBeLoaded_ is empty) */
	protected Collection<Gene> genes_ = null;	
	
	/** The reference population instance, used to load snp positions and genotypes */
	protected ReferencePopulation refpop_ = null;
	
	/** Computes the gene scores */
	protected GeneScoreEvaluator evaluator;

	/** Writes gene scoring results */
	private GeneScoreWriter geneScoreWriter;
	
	/** Flag set false when computeScore() has been run first time (to display info only on first run) */
	private boolean firstRun_ = true;		
	protected double[] fakeSignal_ = null;
	
	/** Directory where correlation matrices are saved */
	private File corMatDir;
	/** Directory where gene SNPs are saved */
	private File geneSnpsDir;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GenomeWideScoring() {
		
		// Initialize gene score evaluator based on the selected method in settings
		if (Pascal.set.useAnalyticVegas_)
			evaluator = new AnalyticVegas();
		else if (Pascal.set.useMaxVegas_)
			evaluator = new MaxSimulAndAnalyticVegas();
		else if (Pascal.set.useMaxEffVegas_)
			evaluator = new MaxEffVegas();
        else
			throw new RuntimeException("No gene scoring method selected");
		
		additionalOutputFileSuffix_ = ConvenienceMethods.addDotAfter(Pascal.set.outputSuffix_) + evaluator.getTypeString();
		if(Pascal.set.useFakeSignal_){
			ReferencePopulationFakeSignal fakeSignalGenerator = new ReferencePopulationFakeSignal();
			fakeSignalGenerator.runFakeSignal();
			fakeSignal_ = fakeSignalGenerator.getSignal();
		}
		
		if (Pascal.set.writeCorFiles_) {
			corMatDir = new File(Pascal.set.outputDirectory_, "correlationMatrices");
			corMatDir.mkdirs();
		}
		if (Pascal.set.writeGenewiseSnpFiles_) {
			geneSnpsDir = new File(Pascal.set.outputDirectory_, "geneSNPs");
			geneSnpsDir.mkdirs();
		}
	}


	// ----------------------------------------------------------------------------

	

	private class ChromosomeSwitchChecker {
		
		String prevChrom_ = "not initialized";				
		public void check(GenomicElement el, ReferencePopulation refpop){
			if (!el.chr_.equals(prevChrom_)){
				
				if (Pascal.set.verbose_) {
					Pascal.println();
					Pascal.println("Chromosome " + el.chr_);
					Pascal.println("----------------");
				} else {
					Pascal.println("- " + el.chr_ + " ...");
				}
				// Initialize refpop
				refpop.initialize(el.chr_);
				if (Pascal.set.verbose_)
					Pascal.println();

				
			}
			prevChrom_ = el.chr_;						
		}
	}

	
	/** 
	 * Compute scores for the given list of genes, which must be ordered by chromosome and starting position.
	 * Gene scores are saved in the gene instances (Gene.score_).
	 * Results are written to file if the flag is set.
	 */
	public void computeScores() {

		if (genes_ == null || genes_.size() == 0) {
			Pascal.warning("GenomeWideScoring.computeScores(): No genes specified");
			return;
		}
		
		//GeneResultsSnpsOutOfBounds_ =  new GeneResultsSnpsOutOfBounds();
		
		if (firstRun_) {
			// Print info on method and parameters used to console
			printConsoleMethodInfo();
			// Load MTJ, display info
			// TODO @David I removed this option, if you want we can only show if Pascal.set.verbose_
			//if (Pascal.set.writeDetailedErrorOutput_){
				initializeMtj();
			//}
		}
		// Open output file
		geneScoreWriter = new GeneScoreWriter(evaluator, additionalOutputFileSuffix_);
		
		// Print header for results displayed on console
		if (Pascal.set.verbose_)
			printConsoleHeader();
		else
			Pascal.println("Computing gene scores for chromosome:");

		// TODO @David refactoring (this class is for nothing)
		ChromosomeSwitchChecker switchChecker = new ChromosomeSwitchChecker();

		// Compute score for each gene
		for (Gene gene : genes_) {
			// Check if we passed on to the next chromosome
			switchChecker.check(gene, refpop_);
			// Compute the score
			if (computeScore(gene))
				geneScoreWriter.writeScore(evaluator, gene);
			
// TODO @David refactoring
//			GeneScoreEvaluator evaluator = null;
//			evaluator = computeScore(gene);
//			// Evaluator is only null if the gene should be skipped (zero snps or exceeding max number of snps)
//			if (evaluator == null)
//				continue;			
//			scores_.writeLine(evaluator, gene);
		}
		geneScoreWriter.close();
		firstRun_ = false;
	}

	
	// ----------------------------------------------------------------------------
	private GeneData setupGeneData(GeneWithItsSnps gene){
		GeneData currentGeneData = null;
		if (Pascal.set.useFakePhenotype_){				
				currentGeneData = new GeneDataFakePhenotype(gene.getSnpList(), fakeSignal_);
			}
			else{
				currentGeneData = new GeneData(gene.getSnpList());			
			}
		return(currentGeneData);
		
	}	
		
	
	/** Compute score for the given gene */
	protected boolean computeScore(Gene gene) {				
		
		// Get the snps that are in the window around the given gene			
		GeneWithItsSnps geneAndSnps = new GeneWithItsSnps(gene, refpop_.findSnps(gene));
					
		if (Pascal.set.verbose_) 
			 geneAndSnps.printGeneNameAndNrOfSnps();			
		
		refpop_.updateLoadedGenotypes(gene);				
		removeLowMafSnps(geneAndSnps.getSnpList());				
		if (!geneAndSnps.checkNrOfSnps()) {
			String line = gene.id_ + "\t" + gene.getSymbolOrNA() + "\t" + geneAndSnps.getNrOfSnps();
			geneScoreWriter.writeSnpsOutOfBounds(line);
			return false;
		}
		
		GeneData currentGeneData = setupGeneData(geneAndSnps);			
		currentGeneData.processData();	
			
		if (geneSnpsDir != null) {		
			File file = new File(geneSnpsDir, "snpVals_gene" +"_"+ gene.symbol_ + "_" + gene.id_ + ".txt");
			currentGeneData.writeGeneSnpsToFile(file);
        }		
		
		 if (corMatDir != null) {
			 File file= new File(corMatDir, "corMat_snpVals_gene" + "_" + gene.symbol_ +"_"+ gene.getId() + ".txt");
			 currentGeneData.writeCovMatToFile(file);
		 }		 
		calculateScore(currentGeneData, gene);
		return true;
	}
	
	public void calculateScore(GeneData currentGeneData, Gene gene){
		
		//GeneScoreEvaluator evaluator = Pascal.set.getGeneScoreEvaluator();
        evaluator.setDataFromGeneData(currentGeneData);
		// Compute gene score
		long t0 = System.currentTimeMillis();
		boolean success = evaluator.computeScore();
		long t1 = System.currentTimeMillis();
		
		// Set gene score
		gene.setScore(evaluator.getScore());
		gene.calcChi2StatFromScore();

		// Print info to console
		if (Pascal.set.verbose_) {
			Pascal.print("\t" + ConvenienceMethods.padRight(Pascal.utils.chronometer(t1 - t0), 22));
			Pascal.print(evaluator.getConsoleOutput());
			Pascal.println();
		}

		// TODO @David But evaluator was not set to null here, so this function would always return not null, correct?
		// What's the point of splitting this function off and why is it public?
		// If the score could not be computed (exception in analytic vegas)
		if (!success) {
			String str = gene.toString();
			str += "\t" + evaluator.getNoScoreOutput();
			//noScores_.add(str);
			geneScoreWriter.writeNoScore(str);
			if (Pascal.set.verbose_)
				Pascal.println(str);
		}

		//return evaluator;
	}
	
	// ----------------------------------------------------------------------------

	/** Print info to console about relevant parameters and which gene scoring method is used */
	private void printConsoleMethodInfo() {
		
		// Print info
		if (Pascal.set.useAnalyticVegas_) {
			Pascal.println("- Gene scoring method: analytic VEGAS");
			if (Pascal.set.snpWeightingDelta_.get(0) != 0.0 || Pascal.set.snpWeightingDelta_.size() > 1)
				Pascal.println("- SNP weighting with delta: " + ConvenienceMethods.array2string(Pascal.set.snpWeightingDelta_, ","));

		} else if (Pascal.set.useSimulationVegas_) {
			Pascal.println("- Gene scoring method: original VEGAS using MC simulation (not recommended)");
//			if (Pascal.set.testStatisticNumSnps_ > 0)
//				Pascal.println("- Test statistic: using up to " + Pascal.set.testStatisticNumSnps_ + " SNPs within gene windows");
//			else
//				Pascal.println("- Test statistic: using all SNPs within gene windows");
		
		} 
		
        else if (Pascal.set.useMaxVegas_) {
            Pascal.println("TBD print max vegas params");
        } 
        else if (Pascal.set.useMaxEffVegas_) {
            Pascal.println("TBD print maxeff vegas params");
        } 
 
        else {
			throw new RuntimeException("No gene scoring method selected");
		}
	}
	
	
	// ----------------------------------------------------------------------------

	/** Print a header for the gene score computation to the console */
	private void printConsoleHeader() {

		if (Pascal.set.useAnalyticVegas_)
			Pascal.println("Format: gene_id, symbol, #snps, runtime, [warnings]");
		else if (Pascal.set.useSimulationVegas_)
			Pascal.println("Format: gene_id, symbol, #snps, #simulations, runtime, [warnings]");		            
        else if (Pascal.set.useMaxVegas_)
            Pascal.println("TBD");
        else if (Pascal.set.useMaxEffVegas_) {
            Pascal.println("TBD");
        } 
   
		else
			throw new RuntimeException("No gene scoring method selected");

		Pascal.println("------------------------------------------------------------------------------");
	}


	// ----------------------------------------------------------------------------

	/** Load MTJ and display info */
	private void initializeMtj() {
		// This is just some bogus operation to load MTJ jni so that potential warning msgs are displayed here
		// TODO is there a way to check if the native libs could be loaded and display a more understandable msg for the user?		
		ConsoleHandler handler = new ConsoleHandler();				
			
			handler.setLevel(Level.FINEST);		
			Logger.getLogger("com.github.fommil.jni.JniLoader").addHandler(handler);
			Logger log =Logger.getLogger("com.github.fommil.jni.JniLoader");
		//		log.addHandler(handler);
			log.setLevel(Level.FINEST);
		
		
		Pascal.println("- Attempting to load native matrix library");
		try {
			SymmDenseEVD.factorize(new UpperSymmDenseMatrix(1));
		} catch (Exception e) { }
		Pascal.println();
		//handler.setLevel(Level.SEVERE);
		//Do not display further info from MTJ		
		//Logger.getLogger("com.github.fommil.jni.JniLoader").setLevel(Level.WARNING);		
		//Logger.getLogger("com.github.fommil.jni.JniLoader").setLevel(Level.FINEST);
	}
	
	public void removeLowMafSnps(Collection<Snp> snpList){
		
		Iterator<Snp> it = snpList.iterator();
		while (it.hasNext()){
			Snp snp = it.next();
			if (snp.getMaf() < Pascal.set.useMafCutoff_ || snp.getAlleleSd()==0)
				it.remove();
		}
		
	}
	// ============================================================================
	// GETTERS AND SETTERS
	public void addAdditionalOutputFileSuffix(String ext) { additionalOutputFileSuffix_ = additionalOutputFileSuffix_ + ext; }
	
	//public Genome getSnps() { return snps_; }
	public Collection<Gene> getGenes() { return genes_; }
	//public  void setGenes(List<Gene> genes) { genes_ = genes; }
	public  void setGenes(Collection<Gene> genes) { genes_ = genes;	}
	public  void setReferencePopulation(ReferencePopulation refpop) { refpop_=refpop;}


	public GeneScoreEvaluator getEvaluator() {
		return evaluator;
	}

}
