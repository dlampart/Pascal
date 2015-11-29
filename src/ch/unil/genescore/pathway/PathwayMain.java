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
package ch.unil.genescore.pathway;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.main.Pascal;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.gpsutils.FileParser;


/**
 * Main class handling pathway / gene set enrichment analysis 
 * (loads data, computes enrichment, writes output).
 */
public class PathwayMain {

	/** Reference to Pascal */
	private Pascal pascal;
	
	/** The name, used for output files (default, extracted from gwas and gene set filenames) */
	protected String name_ = null; 
	/** The library of gene sets */
	private GeneSetLibrary geneSetLib_ = null;
	
	/** The gene score evaluator */
	protected GenomeWideScoring geneScorer_ = null;
	/** All used genes */	
	protected LinkedHashMap<String, Gene> genes_ = null;
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public PathwayMain(Pascal pascal) {
		
		this.pascal = pascal;
		
		// Prefix used for output files
		name_ = extractName(Pascal.set.snpPvalFile_, Pascal.set.geneSetFile_);
		
		Pascal.println("LOADING INPUT FILES");
		Pascal.println("-------------------\n");		

		// Load gene sets and genes, optionally with scores
		genes_ = loadGenes();
		
		
		// Load snps and pvals, initialize reference population		
		geneScorer_ = new GenomeWideScoring();			
		ReferencePopulation myRefPop=new ReferencePopulation();
		myRefPop.loadGwasAndRelevantSnps();	
		
		geneScorer_.setGenes(genes_.values());	
		geneScorer_.setReferencePopulation(myRefPop);
		

	}

	
	// ----------------------------------------------------------------------------

	/** Run enrichment analysis */
	public void run() {

		
			Pascal.println("COMPUTING GENE SCORES");
			Pascal.println("---------------------\n");		

			// Sort genes by position 	
			if (!Pascal.set.loadScoresFromFiles_){
				
				ArrayList<Gene> genes = new ArrayList<Gene>(genes_.values());
				Collections.sort(genes);
			// Compute scores for individual genes (even those that may later be merged, useful for output)
				geneScorer_.setGenes(genes);			
				geneScorer_.computeScores();
			}
			else {
				Pascal.println("Loading genescores form files:");
				Pascal.println(Pascal.set.geneScoreFile_.getPath());
				Pascal.println(Pascal.set.metaGeneScoreFile_.getPath());
				genes_=loadScoresFromFile(Pascal.set.geneScoreFile_);
				
			}
			geneSetLib_ = new GeneSetLibrary(Pascal.set.geneSetFile_, genes_);
			
			// Print info on gene overlap
			Pascal.println();
			Pascal.println("Loaded:");
			Pascal.println("- " + geneSetLib_.getGeneSets().size() + " gene sets");
			Pascal.println("- " + geneSetLib_.getGenes().size() + " genes with pathway / gene set annotation");			
			Pascal.println("Ignored:");
			Pascal.println("- " + geneSetLib_.getNotFoundInGenomeAnnot().size() + " genes not found in loaded genome annotation (" + GeneAnnotation.getAnnotationName() + ")");			
			Pascal.println();	
		
		if (!(Pascal.set.mergeGenesDistance_ < 0)) {
			Pascal.println("CREATING META-GENES");
			Pascal.println("-------------------\n");	
			Pascal.println("Merging genes of the same pathway that are <" + Pascal.set.mergeGenesDistance_ + "mb apart:");
		
			// Merge neighboring genes in meta-genes
			geneSetLib_.createMetaGenes();
		
			// Output info
			int[] minMaxTot = geneSetLib_.countMergedGenes();
			int numMeta = geneSetLib_.getMetaGenes().size();
			Pascal.println("- " + numMeta + " meta-genes created");
			Pascal.println("- " + minMaxTot[2]/(double)numMeta + " genes per meta-gene on average (min=" + minMaxTot[0] + ", max=" + minMaxTot[1] + ")");
			Pascal.println();
			Pascal.println("Computing gene scores for meta-genes ...");

			// Sort meta-genes by position
			ArrayList<Gene> metagenes = new ArrayList<Gene>(geneSetLib_.getMetaGenes());
			if (!Pascal.set.loadScoresFromFiles_){
				Collections.sort(metagenes);
				// Compute scores
				geneScorer_.setGenes(metagenes);
				geneScorer_.addAdditionalOutputFileSuffix(".fusion");
				geneScorer_.computeScores();
				}
			else if(Pascal.set.mergeGenesDistance_ >= 0){
				LinkedHashMap<String, Gene> metaGenesFromFile;
				metaGenesFromFile = loadScoresFromFile(Pascal.set.metaGeneScoreFile_);
				for (Gene metaGene : metagenes){
					String currentId=metaGene.id_;
					if(metaGenesFromFile.containsKey(currentId)){
						metaGene.copyScores(metaGenesFromFile.get(currentId));
					}
					
				}
				
			}
		}		
		Pascal.println("COMPUTING GENE SET ENRICHMENT");
		Pascal.println("-----------------------------\n");		
		Pascal.println("Computing enrichment p-values for " + geneSetLib_.getGeneSets().size() + " gene sets ...");
		
		// Remove genes without scores
		geneSetLib_.removeMetaGenesWithoutScore();
		
		// Compute enrichment
		
		boolean doApproxAnalysis = false;
		if(doApproxAnalysis){
			geneSetLib_.computeApproxPathwayCorrelationLowerBound();
			geneSetLib_.writePathwayCorrelationMat(new File(Pascal.set.outputDirectory_, "pathwayCorLowerBoundMat.txt"));
			geneSetLib_.computeApproxPathwayCorrelation();
			geneSetLib_.writePathwayCorrelationMat(new File(Pascal.set.outputDirectory_, "pathwayCorMat.txt"));
		}

		long t0 = System.currentTimeMillis();
		geneSetLib_.computeEnrichment();
		long t1 = System.currentTimeMillis();
				
		
		Pascal.println("- Runtime: " + Pascal.utils.chronometer(t1-t0));
		Pascal.println("- " + geneSetLib_.countNumSignificantSets() + " gene sets with p-value < " + Pascal.set.writeSignificanceThreshold_);
		Pascal.println();
		
		// Write result
		File file = new File(Pascal.set.outputDirectory_, name_ + ".txt");
		geneSetLib_.writeOutput(file);		
	}

	
	// ============================================================================
	// PRIVATE METHODS
	public void removeGenesWithoutScore() {
	
		Iterator<String> iter = genes_.keySet().iterator();
		while (iter.hasNext()){
			String currentId  = iter.next();
			Gene currentGene = genes_.get(currentId);
		if(currentGene.getChi2Stat() == -1){
			iter.remove();
		}
		}
		
	}

	
	public HashMap<String, Integer>  getHeaderHash(String[] fieldnames){
		
		HashMap<String, Integer> headerIndex = new HashMap<String, Integer>();
		for (int i=0; i<fieldnames.length; i++){
			headerIndex.put(fieldnames[i],i);
		}
		return headerIndex;
	}
	
	public LinkedHashMap<String, Gene> loadScoresFromFile(File file){
		LinkedHashMap<String, Gene> Genes = new LinkedHashMap<String, Gene>();		
		FileParser parser = new FileParser(Pascal.log, file);
		String[] header = parser.readLine();
		if (header == null){
		parser.close();
		return Genes;
		}
		HashMap<String, Integer> hash = getHeaderHash(header);
		
		while(true){			
			String[] nextLine = parser.readLine();			
			if (nextLine == null)			
				break;		
			String geneId = nextLine[hash.get("gene_id")];
			double score = Double.valueOf(nextLine[hash.get("pvalue")]);
			String chr =  nextLine[hash.get("chromosome")];
			String geneSymbol =  nextLine[hash.get("gene_symbol")];
			int start = Integer.valueOf(nextLine[hash.get("start")]);
			int end = Integer.valueOf(nextLine[hash.get("end")]);
			boolean posStrand = nextLine[hash.get("strand")].equals("+");
											
			Gene currentGene = new Gene(geneId,geneSymbol);			
			currentGene.setPosition(chr, start, end, posStrand);			
			currentGene.setScore(score);
			currentGene.calcChi2StatFromScore();
			Genes.put(geneId,currentGene);			
		}
		parser.close();
		return Genes;
	 
		
		
	}
	
	protected LinkedHashMap<String, Gene> loadGenes() {
		
		// The genes from the genome annotation
		LinkedHashMap<String, Gene> genes;		
			// Load genome annotation (genes without scores)
			GeneAnnotation annot = GeneAnnotation.createAnnotationInstance();
			genes = annot.loadAnnotation();
			// Remove excluded genes (this updates genes_)
			if (!Pascal.set.excludedGenesFile_.equals("")) {
				int numRemoved = annot.removeGenes(Pascal.set.excludedGenesFile_);
				Pascal.println("- " + numRemoved + " genes removed because they are in the excluded genes file");
			}

		
	
		return genes;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Extract the name for this run from the snp and pathway annotation files */
	protected String extractName(File geneScoreFile, File geneSetFile) {
		
		String gwasName = Pascal.utils.extractBasicFilename(geneScoreFile.getName(), false);
		String functName = Pascal.utils.extractBasicFilename(geneSetFile.getName(), false);
		
		functName = functName.replace("_nodeProperties", "");
		functName = functName.replace("_undir", "");
		String condDot="";
		if (!Pascal.set.outputSuffix_.equals("")){
			condDot=".";
		}
		return gwasName + ".PathwaySet--" + functName + "--" + Pascal.set.outputSuffix_ + condDot + pascal.getGenomeWideScoring().getEvaluator().getTypeString() + Pascal.set.chromFileExtension_;
	}
	
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public GeneSetLibrary getGeneSetLib() { return geneSetLib_; }
	
}
