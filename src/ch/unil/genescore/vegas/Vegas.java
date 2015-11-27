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
import java.io.IOException;
import java.util.ArrayList;
import java.util.ListIterator;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Pascal;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;

/**
 * Summarize p-values of multiple snps at a locus
 */
public abstract class Vegas extends GeneScoreEvaluator {

	/** The set of snps to be considered (correspond to a gene or locus of interest) */
	//protected ArrayList<Snp> geneSnps_ = null;
	protected ArrayList<Double> snpScores_ = null;
	
	/** The number of gene snps */
	protected int numSnps_ = -1;
	/** The gene score */
	protected double geneScore_ = -1;
	
	/** The covariance matrix (ld matrix) */
	protected UpperSymmDenseMatrix covariance_ = null;
	/** By how much the covariance matrix was changed to make it positive definite (sum of squares, 0 if the original was positive definite) */
	protected double positiveDefiniteCorrection_ = 0;
	/** The tolerance that was used to make the covariance matrix positive definite */
	protected double tolerance_ = 0;
	
	/** The test statistic for the actual data */
	protected double testStatisticReal_ = -1;
	
	protected double[] weights_; 

	
	// ============================================================================
	// ABSTRACT METHODS

	/** Compute test statistic for k'th p-value of the snps */
	protected abstract double computeTestStatisticRealSubclass();

		
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Vegas() { }

	
	/** Constructor */
	public Vegas(ArrayList<Double> snpScores, DenseMatrix ld) {
		
		setSnpScores(snpScores);
		setCovariance(ld);
	}

	/** Constructor */
	public Vegas(ArrayList<Double> snpScores, UpperSymmDenseMatrix ld) {
		
		setSnpScores(snpScores);
		setCovariance(ld);
	}


	// ----------------------------------------------------------------------------

	/** Write a detailed gene report to file */
/*	public void writeGeneReport(String filename) {
		
		Main.println();
		FileExport writer = new FileExport(filename);
		
		// The header
		String header = "";
		for (Snp snp : geneSnps_)
			header += "\t" + snp.id_;
		writer.println(header);
		
		// Positions
		String nextLine = "Position";
		for (Snp snp : geneSnps_)
			nextLine += "\t" + snp.start_;
		writer.println(nextLine);
		
		// P-values
		nextLine = "Pvalue";
		for (Snp snp : geneSnps_)
			nextLine += "\t" + Utils.toStringScientific10(snp.getPval());
		writer.println(nextLine);
		
		// Correlation matrix
		for (int i=0; i<numSnps_; i++) {
			nextLine = geneSnps_.get(i).id_;
			for (int j=0; j<numSnps_; j++)
				nextLine += "\t" + Utils.toStringScientific10(covariance_.get(i, j));
			writer.println(nextLine);
		}
		writer.close();
	}
*/	
		
	// ----------------------------------------------------------------------------

	/** Assign the snps */
	public void setSnpScores(ArrayList<Double> snpScores) {
		
		snpScores_ = snpScores;
		numSnps_ = snpScores_.size();
	}

		
	// ============================================================================
	// PROTECTED METHODS
	
	/** Compute test statistic for all snp p-values (real and phenotype permutations if available) */
	public void computeTestStatisticReal() {
		
		// Transform snp p-values to chi-squared test statistics		
		//for (int i=0; i<numSnps_; i++)
			//snpScores_[i].computeChiSquaredStatistics();
		
		testStatisticReal_ = computeTestStatisticRealSubclass();
	}
	
	
	// ----------------------------------------------------------------------------

//	/**
//	 * Method implemented following corpcor R package function make.positive.definite(). 
//	 * Computes the nearest positive definite of a real symmetric matrix, using the algorithm 
//	 * of NJ Higham (1988, Linear Algebra Appl. 103:103-118).
//	 * 
//	 * tol:
//	 * Tolerance for singular values and for absolute eigenvalues - only those with values larger 
//	 * than tol are considered non-zero (default: tol = max(dim(m))*max(D)*.Machine$double.eps).
//	 * Note that the definition tol= max(dim(m))*max(D)*.Machine$double.eps is exactly compatible 
//	 * with the conventions used in "Octave" or "Matlab".
//	 *
//	 *	make.positive.definite(m, tol) {
//     *		d = dim(m)[1]
//     *		if (dim(m)[2] != d) 
//     *			stop("Input matrix is not square!")	
//     *		es = eigen(m)
//     *		esv = es$values
//     *		if (missing(tol)) 
//     *			tol = d * max(abs(esv)) * .Machine$double.eps
//     *		delta = 2 * tol
//     *		tau = pmax(0, delta - esv) # pmax: max of corresponding elements in the two vectors
//     *		dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
//     *		return(m + dm)
//     *	}
//	 */
//	protected void makeCovariancePositiveDefinite() {
//		
//		int d = numSnps_;
//		assert d == covariance_.getColumnDimension() && d == covariance_.getRowDimension();
//		
//		// Apache commons eigen decomposition		
//		//EigenDecomposition eigen = new EigenDecomposition(matrix);
//		//double[] eigenvalues = eigen.getRealEigenvalues();
//
//		// Colt eigen decomposition
//		DenseDoubleEigenvalueDecomposition eigen = new DenseDoubleEigenvalueDecomposition(Utils.apache2colt(covariance_));
//		DoubleMatrix1D eigenValues = eigen.getRealEigenvalues();
//		RealMatrix V = Utils.colt2apache(eigen.getV());
//		RealMatrix VT = V.transpose();
//		
//		assert eigenValues.size() == d;
//		
//		// Get the leading eigenvalue (max absolute value)
//		double leadingEigenValue = 0;
//		for (int i=0; i<d; i++) {
//			double x = Math.abs(eigenValues.get(i));
//			if (x > leadingEigenValue)
//				leadingEigenValue = x;
//		}
//		
//		// Tolerance
//		tolerance_ = d * leadingEigenValue * Precision.EPSILON;
//		double delta = 2 * tolerance_;
//		DiagonalMatrix tau = new DiagonalMatrix(d);
//		for (int i=0; i<eigenValues.size(); i++)
//			tau.setEntry(i, i, Math.max(0, delta - eigenValues.get(i)));
//
//		// dm = V * tau * V_t
//		RealMatrix dm = tau.preMultiply(V).multiply(VT);
//		assert dm.getColumnDimension() == d && dm.getRowDimension() == d;
//		
//		// Compute the sum of squares (Frobenius norm) of the matrix that will be added to the covariance
//		positiveDefiniteCorrection_ = 0;
//		for (int i=0; i<d; i++)
//			for (int j=0; j<d; j++)
//				positiveDefiniteCorrection_ += Math.pow(dm.getEntry(i, j), 2);
//		
//		// Add dm to the covariance matrix
//		covariance_ = covariance_.add(dm);
//	}

	
	// ----------------------------------------------------------------------------

	/** Compute the mean LD (correlation) of the SNP pairs (mean of upper triangular part of covariance_) */
	protected double meanLd() {

		double sum = 0;
		int count = 0;
		
		for (int i=0; i<numSnps_; i++) {
			for (int j=i+1; j<numSnps_; j++) {
				count++;
				sum += covariance_.get(i, j);
			}
		}
		return (count == 0) ? -1 : sum/count;
	}
    protected void writeMatrixMTJ(AbstractMatrix mat, String outName) {
        
        String filename = Settings.outputDirectory_ + "/" + "run_" + outName + ".ld";
        System.out.println(filename); 
        File file = new File(filename);
        try {
            file.createNewFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
        FileExport writer = new FileExport(filename);
        
        for (int i = 0; i < mat.numRows(); i++){
            String line = "";
            for (int j = 0; j < mat.numColumns(); j++){
                if (j == 0){
                    line += mat.get(i,j);
                }
                else {
                    line += "\t" + mat.get(i,j);
                }
            }
            writer.println(line);
        }
        writer.close();
        Pascal.println("");
    }
    
    //public double[] getZscores(){
    //    double[] zScoresList = new double[snpScores_.length];
    //    Iterator<Double> itr = snpScores_.length;
     //   int count=0;
     //   while(itr.hasNext()){
     //       zScoresList[count]=itr.next().getZscore();
      //      count++;
     //   }
      //  return zScoresList;
  //  }

    public void setDataFromGeneData(GeneData Dat){
    	
    	setCovariance(Dat.getCorr());
    	setWeights(Dat.getWeights());
    	setSnpScores(Dat.getScores());
    }
	// ============================================================================
	// GETTERS AND SETTERS
    
	
    public void setCovariance(DenseMatrix covariance) { covariance_ = new UpperSymmDenseMatrix(covariance); }
    public void setCovariance(UpperSymmDenseMatrix covariance) { covariance_ = covariance; }

    public double[] getWeights() { return weights_; } 
	public void setWeights(double[] t) { weights_ = t; }
    
	public double getTestStatisticReal() { return testStatisticReal_; } 
	public void setTestStatisticReal(double t) { testStatisticReal_ = t; }
	//public double getScore() { return geneScore_; }
		
	public double getPositiveDefiniteCorrection() { return positiveDefiniteCorrection_; }
	//public double[][] getCovariance(){ 
	//	double[][] covs = covariance_.
	//	return covs;
	//}
	
	public double[] getScore() {
		double[] result = { geneScore_ };
		return result;
	}

}
