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
package ch.unil.genescore.pathway;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulation;


/**
 * Main class handling pathway / gene set enrichment analysis 
 * (loads data, computes enrichment, writes output).
 */
public class PathwayMain {

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
	public PathwayMain() {
		
		// Prefix used for output files
		name_ = extractName(Settings.snpPvalFile_, Settings.geneSetFile_);
		
		Main.println("LOADING INPUT FILES");
		Main.println("-------------------\n");		

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
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void run() {

		
			Main.println("COMPUTING GENE SCORES");
			Main.println("---------------------\n");		

			// Sort genes by position 	
			if (!Settings.loadScoresFromFiles_){
				
				ArrayList<Gene> genes = new ArrayList<Gene>(genes_.values());
				Collections.sort(genes);
			// Compute scores for individual genes (even those that may later be merged, useful for output)
				geneScorer_.setGenes(genes);			
				geneScorer_.computeScores();
			}
			else {
				Main.println("Loading genescores form files:");
				Main.println(Settings.geneScoreFile_);
				Main.println(Settings.metaGeneScoreFile_);
				genes_=loadScoresFromFile(Settings.geneScoreFile_);
				
			}
			geneSetLib_ = new GeneSetLibrary(Settings.geneSetFile_, genes_);
			
			// Print info on gene overlap
			Main.println();
			Main.println("Loaded:");
			Main.println("- " + geneSetLib_.getGeneSets().size() + " gene sets");
			Main.println("- " + geneSetLib_.getGenes().size() + " genes with pathway / gene set annotation");			
			Main.println("Ignored:");
			Main.println("- " + geneSetLib_.getNotFoundInGenomeAnnot().size() + " genes not found in loaded genome annotation (" + GeneAnnotation.getAnnotationName() + ")");			
			Main.println();	
		
		if (!(Settings.mergeGenesDistance_ < 0)) {
			Main.println("CREATING META-GENES");
			Main.println("-------------------\n");	
			Main.println("Merging genes of the same pathway that are <" + Settings.mergeGenesDistance_ + "mb apart:");
		
			// Merge neighboring genes in meta-genes
			geneSetLib_.createMetaGenes();
		
			// Output info
			int[] minMaxTot = geneSetLib_.countMergedGenes();
			int numMeta = geneSetLib_.getMetaGenes().size();
			Main.println("- " + numMeta + " meta-genes created");
			Main.println("- " + minMaxTot[2]/(double)numMeta + " genes per meta-gene on average (min=" + minMaxTot[0] + ", max=" + minMaxTot[1] + ")");
			Main.println();
			Main.println("Computing gene scores for meta-genes ...");

			// Sort meta-genes by position
			ArrayList<Gene> metagenes = new ArrayList<Gene>(geneSetLib_.getMetaGenes());
			if (!Settings.loadScoresFromFiles_){
				Collections.sort(metagenes);
				// Compute scores
				geneScorer_.setGenes(metagenes);
				geneScorer_.setAdditionalOutputFileSuffix(".meta");
				geneScorer_.computeScores();
				}
			else if(Settings.mergeGenesDistance_ >= 0){
				LinkedHashMap<String, Gene> metaGenesFromFile;
				metaGenesFromFile = loadScoresFromFile(Settings.metaGeneScoreFile_);
				for (Gene metaGene : metagenes){
					String currentId=metaGene.id_;
					if(metaGenesFromFile.containsKey(currentId)){
						metaGene.copyScores(metaGenesFromFile.get(currentId));
					}
					
				}
				
			}
		}		
		Main.println("COMPUTING GENE SET ENRICHMENT");
		Main.println("-----------------------------\n");		
		Main.println("Computing enrichment p-values for " + geneSetLib_.getGeneSets().size() + " gene sets ...");
		
		// Remove genes without scores
		geneSetLib_.removeMetaGenesWithoutScore();
		
		// Compute enrichment
		
		boolean doApproxAnalysis = false;
		if(doApproxAnalysis){
		geneSetLib_.computeApproxPathwayCorrelationLowerBound();
		geneSetLib_.writePathwayCorrelationMat("", "pathwayCorLowerBoundMat.txt");
		geneSetLib_.computeApproxPathwayCorrelation();
		geneSetLib_.writePathwayCorrelationMat("", "pathwayCorMat.txt");
		}
		
		long t0 = System.currentTimeMillis();
		geneSetLib_.computeEnrichment();
		long t1 = System.currentTimeMillis();
				
		
		Main.println("- Runtime: " + Utils.chronometer(t1-t0));
		Main.println("- " + geneSetLib_.countNumSignificantSets() + " gene sets with p-value < " + Settings.writeSignificanceThreshold_);
		Main.println();
		
		// Write result
		String filename = Settings.outputDirectory_ + "/" + name_ + ".txt";
		geneSetLib_.writeOutput(filename);		
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
	
	public LinkedHashMap<String, Gene> loadScoresFromFile(String filename){
		LinkedHashMap<String, Gene> Genes = new LinkedHashMap<String, Gene>();		
		FileParser parser = new FileParser(filename);
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
	
	protected LinkedHashMap loadGenes() {
		
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
	
	
	// ----------------------------------------------------------------------------

	/** Extract the name for this run from the snp and pathway annotation files */
	protected String extractName(String geneScoreFile, String geneSetFile) {
		
		String gwasName = Utils.extractBasicFilename(geneScoreFile, false);
		String functName = Utils.extractBasicFilename(geneSetFile, false);
		
		functName = functName.replace("_nodeProperties", "");
		functName = functName.replace("_undir", "");
		
		return gwasName + "--" + functName + "--" + Settings.outputSuffix_;
	}
	
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public GeneSetLibrary getGeneSetLib() { return geneSetLib_; }
	
}
