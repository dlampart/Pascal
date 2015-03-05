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
package ch.unil.genescore.gene;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map.Entry;

import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;


/**
 * Gene annotation, provides functionality to load IDs and coordinates of genes
 */
abstract public class GeneAnnotation {

	/** The file with the genome annotation */
	protected String annotationFile_ = null;
	
	/** User specified set of genes to be loaded (leave empty for all genes, boolean indicates if gene was found in the annotation) */
	protected HashMap<String, Boolean> genesToBeLoaded_ = null;
	/** Load only genes from this chromosome (leave empty for all chromosomes) */
	protected String chromosomeToBeLoaded_ = null;
	/** Ignore X and Y chromosomes */
	protected boolean ignoreAllosomes_ = false;
	/** Ignore mitochondiral genes */
	protected boolean ignoreChrM_ = false;
	/** Flag indicates whether only protein coding genes should be loaded */
	protected boolean loadOnlyProteinCoding_ = false;
	
	/** The genes that were loaded from the annotation file */
	protected LinkedHashMap<String, Gene> genes_ = null;
	/** Set of genes that were excluded */
	protected HashSet<String> excludedGenes_ = null;
	
	
	// ============================================================================
	// STATIC METHODS

	/** Create an instance of the appropriate subclass as specified in Settings.genomeAnnotation_ */
	public static GeneAnnotation createAnnotationInstance() {
		
		// TODO allow annotation to be loaded from a bed file, create appropriate subclass
		if (Settings.genomeAnnotation_.equals("gencode"))
			return new GeneAnnotationGencode();
		else if (Settings.genomeAnnotation_.equals("ucsc"))
			return new GeneAnnotationUcsc();
		else if (Settings.genomeAnnotation_.equals("bed"))
			return new GeneAnnotationBed();
			throw new RuntimeException("Settings.genomeAnnotation_ must be 'gencode','ucsc' or 'bed'.");
	}
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneAnnotation(String annotationFile) {
		
		annotationFile_ = annotationFile; 
		chromosomeToBeLoaded_ = Settings.chromosome_;
		ignoreAllosomes_ = Settings.ignoreAllosomes_;
		ignoreChrM_ = Settings.ignoreChrM_;
		loadOnlyProteinCoding_ = Settings.loadOnlyProteinCodingGenes_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Load the specified genes (genesToBeLoaded_) from the annotation file */
	public LinkedList<Gene> loadAnnotation(String genesToBeLoadedFile) {
		
		// Load the set of genes to be considered
		loadGenesToBeLoaded(genesToBeLoadedFile);
		// Load the specified genes from the annotation
		loadAnnotation();
		
		// If not all specified genes were found, print warning
		if (genesToBeLoaded_ != null && genesToBeLoaded_.size() != genes_.size()) {
			int numNotFound = genesToBeLoaded_.size() - genes_.size();
			
			String genesNotFound = "";
			for (Entry<String,Boolean> entry : genesToBeLoaded_.entrySet())
				if (!entry.getValue())
					genesNotFound += entry.getKey() + " ";

			Main.println("   - " + numNotFound + " genes were not found in the annotation: " + genesNotFound);
		}
		//Don't need the hashmap put into list
		LinkedList<Gene> genes = new LinkedList<Gene>();
		for(Gene gene : genes_.values()){
			genes.add(gene);
		}
		Collections.sort(genes);
		return genes;
	}

	// ----------------------------------------------------------------------------

	/** Load the specified genes (genesToBeLoaded_) from the annotation file */
	abstract public LinkedHashMap<String, Gene> loadAnnotation();

	
	// ----------------------------------------------------------------------------

	/** Load the specified set of genes (genesToBeLoaded_) */
	public void loadGenesToBeLoaded(String filename) {
				
		// Return if no gene file was specified (all genes from the annotation will be loaded)
		if (filename == null ||
				filename.equals(" ") ||
				filename.equals(""))
			return;

		genesToBeLoaded_ = new HashMap<String, Boolean>();
		FileParser parser = new FileParser(filename);
		
		while (true) {
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			if (nextLine.length != 1)
				parser.error("Expected one column (gene ID)");
			
			genesToBeLoaded_.put(nextLine[0].toUpperCase(), false);
		}
		parser.close();
	}


	// ----------------------------------------------------------------------------

	/** Write the genes to a file with gene id, symbol and position */
	public void writeGeneList(String filename) {
		
		FileExport writer = new FileExport(filename);
		String prevChr = null;
		int prevStart = -1;
		
		String header = "geneId\tgeneSymbol\tchromosome\tstart\tend\tstrand";
		writer.println(header);
		
		for (Gene gene : genes_.values()) {
			if (prevStart == -1 || prevChr != gene.chr_) {
				prevStart = gene.start_;
				prevChr = gene.chr_;
			}
			
			if (gene.start_ < prevStart)
				Main.error("Genes are not ordered by genomic position");
			
			String nextLine = gene.id_ + "\t" + gene.symbol_ + "\t" + 
					gene.chr_ + "\t" + gene.start_ + "\t" + gene.end_ + "\t" + 
					(gene.posStrand_ ? "+" : "-");
			
			writer.println(nextLine);
		}
		writer.close();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Return true if this chromosome should be skipped */
	public boolean skipChromosome(String chr) {
		
		boolean notSpecifiedChrom = chromosomeToBeLoaded_ != null && !chromosomeToBeLoaded_.equals("") && !chr.equals(chromosomeToBeLoaded_);
		boolean ignoreChrom = (ignoreAllosomes_ && (chr.equals("chrX") || chr.equals("chrY"))) ||
				(ignoreChrM_ && chr.equals("chrM"));
		boolean isValid = Chromosome.isValidId(chr);
		
		return !isValid || notSpecifiedChrom || ignoreChrom;
	}
	

	// ----------------------------------------------------------------------------

	/** Remove the genes in the given file from the loaded set of genes */
	public int removeGenes(String excludedGenesFile) {
		
		excludedGenes_ = loadGeneList(excludedGenesFile);
		int numRemoved = 0;
		
		for (String gene : excludedGenes_) {
			Gene g = genes_.remove(gene);
			if (g != null)
				numRemoved++;
		}
		return numRemoved;
	}

	
	// ----------------------------------------------------------------------------
	
	/** Load a list of genes (one line per gene, id in first column) */
	public static HashSet<String> loadGeneList(String excludedGenesFile) {
		
		HashSet<String> excludedGenes = new HashSet<String>();
		
		if (excludedGenesFile == null || excludedGenesFile.length() == 0)
			return excludedGenes;
		
		FileParser parser = new FileParser(excludedGenesFile);

		// Parse header
		String[] header = parser.readLine();
		if (!header[0].equals("gene_id"))
			Main.error("Expected header line with first field 'gene_id' (tab-separated)");
		
		// First line
		while (true) {
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			
			// Gene id
			String id = nextLine[0];
			id = GeneAnnotationGencode.removeEnsemblVersion(id);
			excludedGenes.add(id);
		}
		parser.close();
		
		return excludedGenes;
	}

	
	/** Get name of genome annotation for console output ('UCSC known gene', 'GENCODE genes') */
	public static String getAnnotationName() {
		
		if (Settings.genomeAnnotation_.equalsIgnoreCase("ucsc"))
			return "UCSC known genes";
		else if (Settings.genomeAnnotation_.equalsIgnoreCase("gencode"))
			return "GENCODE genes";
		else
			throw new RuntimeException("Unknown genomeAnnotation:" + Settings.genomeAnnotation_);
	}
	
	
	// ============================================================================
	// PROTECTED METHODS

	/** Check if the given gene ID or symbol are listed in genesToBeLoaded_, if yes return true and set the corresponding entry true */
	protected boolean checkGenesToBeLoaded(String geneId, String geneName) {
	
		if (genesToBeLoaded_ == null)
			return true;
		
		// Check if this gene is part of the specified gene set
		String specifiedGene = null;
		if (genesToBeLoaded_.containsKey(geneId))
			specifiedGene = geneId;
		else if (genesToBeLoaded_.containsKey(geneName))
			specifiedGene = geneId;

		// Flag this gene as found
		if (specifiedGene != null) {
			genesToBeLoaded_.put(specifiedGene, true);
			return true;

		} else {
			return false;
		}
	}
	
}
