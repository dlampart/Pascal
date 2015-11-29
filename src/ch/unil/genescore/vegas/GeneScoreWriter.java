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

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileExport;

/**
 * Write gene scoring results / errors to files
 */
public class GeneScoreWriter {
	
	/** Suffix appended to filename right before the file extension */
	private String outputFileSuffix;
	
	/** File for genes with scores */
	private FileExport scores;
	/** File for genes without scores */
	private FileExport noScores;
	/** File for genes with zero or too many SNPs */
	private FileExport snpsOutOfBounds;
	
	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor */
	public GeneScoreWriter(GeneScoreEvaluator evaluator, String outputFileSuffix) {
		
		this.outputFileSuffix = outputFileSuffix;
		
		// Create the file for the gene scores (other files will only be created if needed)
		File file = new File(Pascal.set.outputDirectory_, 
				Pascal.set.gwasName_ + ConvenienceMethods.addDotBefore(outputFileSuffix) 
				+ ".genescores" + Pascal.set.chromFileExtension_  + ".txt");
		scores = new FileExport(Pascal.log, file);
		scores.println("chromosome\tstart\tend\tstrand\tgene_id\tgene_symbol" + evaluator.getResultsAsStringHeader());
	}
	
	
	// ----------------------------------------------------------------------------

	/** Close files, print info */
	public void writeScore(GeneScoreEvaluator evaluator, Gene gene) {

		scores.println(gene.toString() + evaluator.getResultsAsString());
		scores.flush();
	}

	
	// ----------------------------------------------------------------------------

	/** Close files, print info */
	public void writeNoScore(String str) {

		// Open file and write header
		if (noScores == null) {
			File file = new File(Pascal.set.outputDirectory_, 
					Pascal.set.gwasName_ + ConvenienceMethods.addDotBefore(outputFileSuffix) 
					+ ".scoreComputeError" + Pascal.set.chromFileExtension_ + ".txt");		
			noScores = new FileExport(Pascal.log, file);
			noScores.println("chromosome\tstart\tend\tstrand\tgene_id\tsymbol\tScore\tStatus");		
		}
		// Write line
		noScores.println(str);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Close files, print info */
	public void writeSnpsOutOfBounds(String str) {

		// Open file and write header
		if (snpsOutOfBounds == null) {
			File file = new File(Pascal.set.outputDirectory_,
					Pascal.set.gwasName_ + ConvenienceMethods.addDotBefore(outputFileSuffix) 
					+ ".numSnpError" + Pascal.set.chromFileExtension_ + ".txt");
			snpsOutOfBounds = new FileExport(Pascal.log, file);
			snpsOutOfBounds.println("gene_id\tsymbol\tSNPs");
		}
		// Write line
		noScores.println(str);
	}

	
	// ----------------------------------------------------------------------------

	/** Close files, print info */
	public void close() {
		
		// Score
		scores.close();
		// No score
		if (noScores != null) {
			Pascal.warning("Gene score computation did not converge at specified precision for some genes");
			noScores.close();
		}
		// Out of bounds
		if (snpsOutOfBounds != null) {
			Pascal.warning("Writing genes without SNPs or exceeding the maximum number of SNPs");
			snpsOutOfBounds.close();
		}
	}
	
	
	// ============================================================================
	// PRIVATE METHODS

	
}
