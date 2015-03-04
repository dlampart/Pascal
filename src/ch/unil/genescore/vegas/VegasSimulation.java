/*
Copyright (c) 2013 Daniel Marbach

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper available at:
http://networkinference.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
package ch.unil.genescore.vegas;

import java.util.ArrayList;
import java.util.PriorityQueue;

import no.uib.cipr.matrix.DenseMatrix;

import org.apache.commons.math3.random.CorrelatedRandomVectorGenerator;

import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;


/**
 * Summarize p-values of multiple snps at a locus
 */
public class VegasSimulation extends Vegas {
	
	/** The number of random samples that were used to compute the empirical p-values */
	protected int numSamples_ = 0;
	/** Random vector generator for multivariate normal distribution */
	protected CorrelatedRandomVectorGenerator multivariateNormalRng_ = null;
	
	/** The number of snps to be used for the test statistic (<= numSnps, depends on Settings.testStatisticNumSnps) */
	protected int numSnpsTestStat_ = -1;
	/** Compute the gene p-value */
	protected EmpiricalPvalue genePval_ = null;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public VegasSimulation(ArrayList<Snp> geneSnps, DenseMatrix ld) {
		
		// overrides the initialization of the super class
		setGeneSnps(geneSnps);
		setCovariance(ld);
	}


	// ----------------------------------------------------------------------------

	/** Compute the score for this snp set */
	public boolean computeScore() {
				
		// If there is only one snp, assign it's p-values to the gene and return
		if (numSnps_ == 1) {
			geneScore_ = geneSnps_.get(0).getPval();
			Main.print("\t0"); // zero simulations were done
			return true;
		}
		
		// Compute the test statistic for the real p-values and the phenotype permutations if available
		computeTestStatisticReal();
		// Compute empirical p-values
		computeEmpiricalPvalues();
		
		// Get result
		geneScore_ = genePval_.getPval();
		return true;
	}

	
	// ----------------------------------------------------------------------------

	/** Perform Monte Carlo simulation to compute empirical p-values */
	public void computeEmpiricalPvalues() {
		
		// Initialize the empirical p-value instances
		genePval_ = new EmpiricalPvalue(testStatisticReal_);

		// For each stage (increasing number of samples)
		for (int i=0; i<Settings.adaptiveNumSamples_.size(); i++) {
			numSamples_ = Settings.adaptiveNumSamples_.get(i);
			Main.print((i == 0) ? "\t" : ",");
			Main.print(Utils.toStringScientific(numSamples_));
			
			// Initialize empirical pval
			genePval_.start();

			// Sample from the null model
			for (int s=0; s<numSamples_; s++)
				sampleTestStatistic();

			// Check stopping criteria
			if (genePval_.stop(numSamples_))
				break;
		}
	}

	
	// ----------------------------------------------------------------------------

	/** Get a string representation of the results */
	public String getResultsAsString() {
			
		String line = "\t" + Utils.toStringScientific10(geneScore_) + "\t" + numSnps_;
		// The number of samples in scientific notation
		line += "\t" + Utils.toStringScientific(numSamples_);
		
		if (Settings.writeDetailedOutput_)
			line += "\t" + Utils.toStringScientific10(meanLd());
		
//		// Comma-separated list of snps used
//		line += "\t" + geneSnps_.get(0).id_;
//		for (int i=1; i<numSnps_; i++)
//			line += "," + geneSnps_.get(i).id_;
		
		return line;
	}

	
	// ----------------------------------------------------------------------------

	/** Get a header line corresponding to getResultsAsString() */
	public String getResultsAsStringHeader() {

		String header = "\tpvalue\tnumSnps\tsimulations";
		if (Settings.writeDetailedOutput_)
			header += "\tavgSnpCorrelation";
		return header;
	}
	
		
	// ----------------------------------------------------------------------------

	/** Set the covariance matrix and initialize the rng */
	public void setCovariance(DenseMatrix covariance) {
		
		throw new RuntimeException("TBD");
//		// Set the covariance matrix
//		super.setCovariance(covariance);
//		//TODO: This will crash!! I just added a RealMatrix-cast to trick the compiler 
//		// Multivariate normal vector generator
//		try {
//			// 1e-12: diagonal elements threshold under which columns are considered to be dependent on previous ones and are discarded
//			// This allows some positive semidefinite matrices to pass that wouldn't when using 0 or Precision.EPSILON (1e-16)
//			multivariateNormalRng_ = new CorrelatedRandomVectorGenerator((RealMatrix) covariance_, 1e-12, new GaussianRandomGenerator(Settings.wellRng_));
//			
//		} catch (NonPositiveDefiniteMatrixException exception) {
//			makeCovariancePositiveDefinite();
//			// Now the matrix should be positive definite and we can use 0 as diagonal element threshold...
//			multivariateNormalRng_ = new CorrelatedRandomVectorGenerator((RealMatrix) covariance_, tolerance_, new GaussianRandomGenerator(Settings.wellRng_));			
//		}
	}
	

	// ----------------------------------------------------------------------------

	/** Assign the snps */
	public void setGeneSnps(ArrayList<Snp> geneSnps) {
		
		super.setGeneSnps(geneSnps);
		
		if (Settings.testStatisticNumSnps_ > 0 && Settings.testStatisticNumSnps_ < numSnps_)
			numSnpsTestStat_ = Settings.testStatisticNumSnps_;
		else
			numSnpsTestStat_ = numSnps_;
	}

		
	// ============================================================================
	// PRIVATE METHODS
	
	/** Sample the test statistic under the null model */
	private void sampleTestStatistic() {
				
		// Generate a new random vector
		double[] x = multivariateNormalRng_.nextVector();
		// Take square
		for (int k=0; k<x.length; k++)
			x[k] = Math.pow(x[k], 2);
		
		// Compute test statistic / sum
		double testStat = sumTopN(x);
		// Add it to the pvalue instances
		genePval_.addSample(testStat);
	}

	
	// ----------------------------------------------------------------------------

	/** Compute test statistic for k'th p-value of the snps */
	protected double computeTestStatisticRealSubclass() {
		if(true)
			throw new RuntimeException("make sure we sorted out chisquare zscores confusion before comming here");
		
		double[] chi2 = new double[numSnps_];
		for (int i=0; i<numSnps_; i++)
			chi2[i] = geneSnps_.get(i).getChi2Stat();
		
		return sumTopN(chi2);
	}


	// ----------------------------------------------------------------------------

	/** Sum the top numSnpsTestStat_ values */
	private double sumTopN(double[] chi2) {
		
		PriorityQueue<Double> topN = null;
		if (numSnps_ > numSnpsTestStat_)
			topN = getTopN(chi2);
		
		double testStat = 0;
		if (topN != null) {
			for (int i=0; i<numSnpsTestStat_; i++)
				testStat += topN.poll();
			assert topN.size() == 0;
			
		} else {
			for (int i=0; i<chi2.length; i++)
				testStat += chi2[i];
		}
		
		return testStat;
	}

	
	// ----------------------------------------------------------------------------

	/** Compute test statistic for k'th p-value of the snps */
	private PriorityQueue<Double> getTopN(double[] v) {
		
		PriorityQueue<Double> topN = new PriorityQueue<Double>(numSnpsTestStat_);
		int i = 0;
		for (; i<numSnpsTestStat_; i++)
			topN.add(v[i]);
		
		for (; i<v.length; i++) {
			if (v[i] > topN.peek()) {
				topN.poll();
				topN.add(v[i]);
			}
		}
		assert topN.size() == numSnpsTestStat_;
		return topN;
	}
	
	
	// ============================================================================
	// GETTERS AND SETTERS
			
	
}
