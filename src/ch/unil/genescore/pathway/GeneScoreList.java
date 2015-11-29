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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileParser;


/**
 * 
 */
public class GeneScoreList {

	/** Genes ranked by score */
	private ArrayList<Gene> rankedGenes_ = null;
	/** Flag defines whether genes are sorted in ascending or descending order */
	private boolean ascending_ = true;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor loading genes from file */
	public GeneScoreList(File geneScoreFile, File excludedGenesFile) {
			
		// Load genes and scores
		loadGeneScores(geneScoreFile, excludedGenesFile);
	}

	
	/** Constructor creating list from given gene set */
	public GeneScoreList(Collection<Gene> genes, boolean ascending) {
		
		rankedGenes_ = new ArrayList<Gene>(genes);
		ascending_ = ascending;
		sortGeneList(0);
	}

	
	// ----------------------------------------------------------------------------

	/** Remove genes that are not in the given gene set, return list of removed genes */
	public ArrayList<String> intersect(Set<String> geneSet) {
		
		// The genes not in the given set
		ArrayList<String> removedGenes = new ArrayList<String>();
		
		Iterator<Gene> iter = rankedGenes_.iterator();
		while (iter.hasNext()) {
			Gene gene = iter.next();
			
			if (!geneSet.contains(gene.id_)) {
				removedGenes.add(gene.id_);
				iter.remove();
			}
		}
		return removedGenes;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Sort the ranked gene list according to the i'th gene score */
	public void sortGeneList(int i) {
		
		// Comparator to sort genes by score
		final class GeneScoreComparator implements Comparator<Gene> {
			// The index of the gene score used for the comparision
			private int scoreIndex_ = -1;
			// Ascending or descending
			private int ascendingSign = ascending_ ? 1 : -1;
			// Constructor
			GeneScoreComparator(int scoreIndex) {
				scoreIndex_ = scoreIndex;
			}
			// Compare
			@Override
			public int compare(Gene g1, Gene g2) {
				return ascendingSign * Double.compare(g1.getScore(scoreIndex_), g2.getScore(scoreIndex_));
			}
		}
		Collections.sort(rankedGenes_, new GeneScoreComparator(i));
	}

	
	// ============================================================================
	// PRIVATE METHODS
		
	/** 
	 * Load the gene scores, initialize genes_. 
	 * Removes genes in excludedGenes_. 
	 * Translates the gene score ids to idTypeGeneSets if necessary. (TODO, see below)
	 */
	private void loadGeneScores(File geneScoreFile, File excludedGenesFile) {
		
		// Load genes to be excluded
		HashSet<String> excludedGenes = GeneAnnotation.loadGeneList(excludedGenesFile);

//		boolean translateToEntrez = Settings.idTypeGeneSets_.equalsIgnoreCase("entrez") && 
//									Settings.idTypeGeneScores_.equalsIgnoreCase("ensembl");
//		if (translateToEntrez)
//			// TODO Can lead to several genes with same ID, screws up LinkedHashMap approach in PathwayMain
//			throw new NotImplementedException();
		
		rankedGenes_ = new ArrayList<Gene>();
		FileParser parser = new FileParser(Pascal.log, geneScoreFile);
		
		// Parse header, check format: chromosome	start	end	strand	gene_id	gene_symbol	pvalue
		String[] header = parser.readLine();
		int numCol = header.length;
		if (numCol < 7)
			parser.error("Expected at least 7 columns");
		if (!header[0].equals("chromosome"))
			parser.error("Expected header line with column 1 'chromosome'");
		if (!header[1].equals("start"))
			parser.error("Expected header line with column 2 'start'");
		if (!header[2].equals("end"))
			parser.error("Expected header line with column 3 'end'");
		if (!header[3].equals("strand"))
			parser.error("Expected header line with column 4 'strand'");
		if (!header[4].equals("gene_id"))
			parser.error("Expected header line with column 5 'gene_id'");
		if (!header[5].equals("gene_symbol"))
			parser.error("Expected header line with column 6 'gene_symbol'");
		if (!header[6].equals("pvalue"))
			parser.error("Expected header line with column 7 'pvalue'");
		
		int numExcluded = 0;
		int numNoScore = 0;
		int numNoAnnot = 0;
		
		while (true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;

			// Check number of lines
			if (nextLine.length != header.length)
				parser.error("Line does not have the same number of columns as the header");
			
			// Gene id
			String geneId = nextLine[4];
			
			// Check if it should be excluded
			if (excludedGenes.contains(geneId)) {
				numExcluded++;
				continue;
			}

			// Parse score. TODO add support for multiple gene scores
			double[] score = new double[1];
			try {
				score[0] = Double.valueOf(nextLine[6]);
			} catch (Exception e) {
				numNoScore++;
				continue;
			}
			
			// Parse position
			String chr = nextLine[0];
			int start = Integer.parseInt(nextLine[1]);
			int end = Integer.parseInt(nextLine[2]);
			boolean posStrand = GenomicElement.isPosStrand(nextLine[3]);

			// Translate from ensembl to entrez
			Collection<String> idSet;
//			if (translateToEntrez) {
//				idSet = GeneIdMapping.getInstance().ensembl2entrez(geneId);
//			} else {
				idSet = new ArrayList<String>(1);
				idSet.add(geneId);
//			}
			
			for (String id : idSet) {
				Gene gene = new Gene(id, nextLine[5]);
				gene.setPosition(chr, start, end, posStrand);
				gene.setScore(score);
				rankedGenes_.add(gene);
			}
		}
		sortGeneList(0);
		
		if (numExcluded > 0)
			Pascal.println("- " + numExcluded + " genes excluded");
		if (numNoScore > 0)
			Pascal.println("- " + numNoScore + " genes with NA score");
		if (numNoAnnot > 0)
			Pascal.println("- " + numNoAnnot + " genes not found in the genome annotation");
	}
	
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public ArrayList<Gene> getGenes() { return rankedGenes_; }
	public Gene getGene(int i) { return rankedGenes_.get(i); }
	public int getNumGenes() { return rankedGenes_.size(); }
	
}
