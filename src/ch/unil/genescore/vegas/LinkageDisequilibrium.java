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

import java.util.ArrayList;
import java.util.Collection;

import no.uib.cipr.matrix.DenseMatrix;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import ch.unil.genescore.main.Main;


/**
 * Linkage disequilibrium (LD) computation.
 */
public class LinkageDisequilibrium {

	
	
	// ============================================================================
		// PUBLIC METHODS
			
	
	// ============================================================================
		// STATIC METHODS
		
	/** Compute correlation (LD) r values for the given snps (use snpCorrelation_ to avoid recomputing previous values) */
	static public RealMatrix computeCorrelationMatrix(ArrayList<Snp> geneSnps) {

		// The correlation matrix for the given set of snps
		RealMatrix corr = MatrixUtils.createRealMatrix(geneSnps.size(), geneSnps.size());

		// For each snp pair
		for (int i=0; i<geneSnps.size(); i++) {
			for (int j=i+1; j<geneSnps.size(); j++) {
				Snp snp_i = geneSnps.get(i);
				Snp snp_j = geneSnps.get(j);

				double r = computeCorrelation(snp_i, snp_j);
				corr.setEntry(i, j, r);
				corr.setEntry(j, i, r);
			}
		}	
		
		// Set diagonal to 1
		for (int i=0; i<geneSnps.size(); i++)
			corr.setEntry(i, i, 1);

		return corr;
	}
	/** Compute rectangular 'cross'-correlation (LD) r values for the two given snp lists( */
	static public RealMatrix computeCrossCorrelationMatrix(ArrayList<Snp> geneSnps1,ArrayList<Snp> geneSnps2) {

		// The correlation matrix for the given set of snps
		RealMatrix corr = MatrixUtils.createRealMatrix(geneSnps1.size(), geneSnps2.size());
		// For each snp pair
		for (int i=0; i<geneSnps1.size(); i++) {
			for (int j=0; j<geneSnps2.size(); j++) {
				Snp snp_i = geneSnps1.get(i);
				Snp snp_j = geneSnps2.get(j);
				double r = computeCorrelation(snp_i, snp_j);
				corr.setEntry(i, j, r);				
			}
		}					
		return corr;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Compute the LD / correlation (r) between the two snps 
	 * correlation matrix is calculated by dividing by n instead n-1
	 * because this leads to a positive semi-definite matrix for sure.
	 * */
	static public double computeCorrelation(Snp snp1, Snp snp2) {
		
		byte[] g1 = snp1.getGenotypes();
		byte[] g2 = snp2.getGenotypes();
		if (g1 == null || g2 == null)
			throw new RuntimeException("Genotypes have not been loaded");
		
		int n = g1.length;
		if (g2.length != n || Snp.getGenotypeLength() != n)
			Main.error("Genotype vectors must be same length");
		
		double preR = 0;
		// if missing data present algorithm should behave as it would with mean-imputation
		for (int i=0; i<n; i++)
			if (g1[i] != -9 && g2[i] != -9)
			preR += g1[i] * g2[i];
			else if (g1[i] != -9 && g2[i] == -9)
				preR += g1[i]*snp2.getAlleleMean();
			else if (g1[i] == -9 && g2[i] != -9)
				preR += snp1.getAlleleMean()*g2[i];
			else // (g1[i] == -9 && g2[i] == -9)
				preR += snp1.getAlleleMean()*snp2.getAlleleMean();
		double r = preR/g1.length;					
		r= (r - snp1.getAlleleMean()*snp2.getAlleleMean())/(snp1.getAlleleSd()*snp2.getAlleleSd());
		return r;
	}
	
	/** loads genotype form List into a DenseMatrix object. it also normalizes the data  and mean imputes*/
	static public DenseMatrix loadGenotypeIntoDenseMatrix(ArrayList<Snp> geneSnps){
		int numberOfSnps = geneSnps.size();
		int lengthOfGenotype = Snp.getGenotypeLength();
		double[][] dataArray = new double[lengthOfGenotype][numberOfSnps];
		double genotype;
		for (int j=0; j < numberOfSnps; ++j){
			Snp currentSnp = geneSnps.get(j);
			byte[] currentGenotype = currentSnp.getGenotypes();
			if (currentGenotype == null)
				throw new RuntimeException();
			for (int i=0; i < lengthOfGenotype; ++i){
				
				if (currentGenotype[i] == -9)
					genotype = currentSnp.getAlleleMean();
				else
					genotype = currentGenotype[i];
				dataArray[i][j] = (genotype-currentSnp.getAlleleMean())/currentSnp.getAlleleSd();
			}
		}
		return new DenseMatrix(dataArray);
	}
	static public DenseMatrix computeCorrelationMatrixMTJ(ArrayList<Snp> geneSnps){
		DenseMatrix normalizedDataMat = loadGenotypeIntoDenseMatrix(geneSnps);
		DenseMatrix  correlationMat = new DenseMatrix(geneSnps.size(),geneSnps.size());
		double invNumberOfRows = 1.0/((double) normalizedDataMat.numRows());		
		normalizedDataMat.transAmult(invNumberOfRows, normalizedDataMat, correlationMat);		
		return correlationMat;
	} 
	static public DenseMatrix computeCrossCorrelationMatrixMTJ(ArrayList<Snp> geneSnpsDesired, ArrayList<Snp> geneSnpsKnown) {
		DenseMatrix normalizedGenotypeMatKnown = loadGenotypeIntoDenseMatrix(geneSnpsKnown);
		DenseMatrix normalizedGenotypeMatDesired = loadGenotypeIntoDenseMatrix(geneSnpsDesired);
		DenseMatrix  crossCorrelationMat = new DenseMatrix(geneSnpsDesired.size(),geneSnpsKnown.size());
		double invNumberOfRows = 1.0/((double) normalizedGenotypeMatKnown.numRows());
		normalizedGenotypeMatDesired.transAmult(invNumberOfRows, normalizedGenotypeMatKnown, crossCorrelationMat);
		return crossCorrelationMat;
	
	}		
	
	static public DenseMatrix getZscores(Collection<Snp> geneSnps){
		DenseMatrix ZscoreMat = new DenseMatrix(geneSnps.size(),1);
		int count=0;
		for (Snp snp : geneSnps){			
			ZscoreMat.set(count, 0, snp.getZscore());
			count++;
		}
		return ZscoreMat;
		
	}
	
	
	// ============================================================================
	// PRIVATE METHODS
		

	// ============================================================================
	// GETTERS AND SETTERS

	
}
