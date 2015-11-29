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
package ch.unil.genescore.gene;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;



import java.util.HashSet;

import ch.unil.genescore.main.Pascal;
import ch.unil.genescore.vegas.Snp;
import ch.unil.genescore.vegas.DistributionMethods;


/**
 * A gene and it's score(s) / p-value(s)
 */
public class Gene extends GenomicElement {

	/** The gene name / symbol (use GenomicElement.id_ for the ensembl/entrez ID) */
	public String symbol_ = null;
	
	/** 
	 * The association scores / p-values.
	 * Keep this as a vector, we could potentially do pathway analysis over multiple gene scores, 
	 * e.g. obtained using different snp weighting strategies.
	 */
	private double[] score_ = null;
	
	/** score that has to been rank-normalized  */
	private double normalizedScore_ = 1;
	
	/** Score mapped to chi2 stat */
	private double chi2Stat_ = -1;
	
	/** sampling weight; only used for pathway methods.*/
	private double samplingWeight_ = 0;

	
	// ============================================================================
	// STATIC METHODS

	/** Create a list of genes sorted by position */
	static public ArrayList<Gene> sortByPosition(Collection<Gene> genes) {
	
		ArrayList<Gene> genesSorted = new ArrayList<Gene>(genes);
		Collections.sort(genesSorted);
		return genesSorted;
	}
	

	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Gene(String id) {

		super(id);
	}

	
	/** Constructor */
	public Gene(String id, String symbol) {
		
		super(id);
		symbol_ = symbol;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the TSS (start if on + strand, end if on - strand) */
	public int getTss() {
		
		if (posStrand_)
			return start_;
		else
			return end_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the the symbol or "NA" if the symbol is NULL */
	public String getSymbolOrNA() {
		
		return (symbol_ == null) ? "NA" : symbol_;
	}

	
	// ----------------------------------------------------------------------------

	/** String representation of the gene */
	public String toString() {
		
		return chr_ + "\t" + start_ + "\t" + end_ + "\t" + getPosStrandStr() + "\t" + id_ + "\t" + getSymbolOrNA();
	}
	
	
	// ----------------------------------------------------------------------------
	/** Get the snps that are in the window around the given gene */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public ArrayList<Snp> findSnps(Genome snps) {
		
		// Define the window around the gene (depends on orientation)
		int start;
		int end;
		if (posStrand_) {
			start = start_ - Pascal.set.geneWindowUpstream_;
			end = end_ + Pascal.set.geneWindowDownstream_;
		} else {
			start = start_ - Pascal.set.geneWindowDownstream_;
			end = end_ + Pascal.set.geneWindowUpstream_;
		}
		
		// Get the set of snps assigned to this gene as an array
		ArrayList<Snp> geneSnps = (ArrayList) snps.getElementsIn(chr_, start, end);
		return geneSnps;
	}
	
	//@SuppressWarnings({ "rawtypes", "unchecked" })
	public HashSet<Snp> findSnpsHash(Genome snps) {
		
		HashSet<Snp> hashset = new HashSet<Snp>();
		ArrayList<Snp> geneSnps = findSnps(snps);
		for(Snp snp : geneSnps)
			hashset.add(snp);
		return(hashset);
	}

	
	public void copyScores(Gene inputGene){
		setScore(inputGene.getScore());
		setChi2Stat(inputGene.getChi2Stat());
		
	}

	public void calcChi2StatFromScore() { 	
		chi2Stat_ = DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(score_[0]);	
	}
	// ============================================================================
	// GETTERS AND SETTERS
		
	public void setScore(double[] x) { score_ = x; }
	public void setScore(double x) { 
		double[] xar = new double[1];
		xar[0]=x;
		score_ = xar; }
	public double[] getScore() { return score_; }
	public double getScore(int i) { return score_[i]; }
	public double getNormalizedScore() { return normalizedScore_; }
	public void setNormalizedScore(double normalizedScore) { normalizedScore_=normalizedScore;}

	
	public void setChi2Stat(double x) {
		chi2Stat_ = x; }
	public double getChi2Stat() { return chi2Stat_; }
	public void setsamplingWeight(double x) {
		samplingWeight_ = x; }
	public double getSamplingWeight(){return samplingWeight_;}
	public void updateSamplingWeight(double x){samplingWeight_=samplingWeight_ + x;} 

	public ArrayList<String> getSymbolList(){
		ArrayList<String> myList = new ArrayList<String>();
		myList.add(symbol_);
		return(myList);
	}
}
