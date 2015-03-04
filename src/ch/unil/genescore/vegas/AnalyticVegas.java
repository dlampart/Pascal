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
import java.util.Collections;
import java.util.Iterator;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;


/**
 * Analytic version of VEGAS
 */
public class AnalyticVegas extends Vegas {

	/**
	 * Convergence Status enum: 
	 * CONV: convergence
	 * ALMOST_CONV:low error but not not much smaller than value (almost convergence)
	 * NOT_CONV:non-convergence
	 * NOT_RUN: not run (because equivalent to only one snp)
	*/
	protected enum Status {CONVERGED, FAIL, DAVIES_SUCCESS, DAVIES_FAIL_FAREBROTHER_SUCCESS,
			DAVIES_FAIL_FAREBROTHER_FAIL, DAVIES_LOWPRECISION_FAREBROTHER_FAIL,
			DAVIES_LOWPRECISION_FAREBROTHER_SUCCESS, NOT_RUN};
	
	/** Convergence status */
	protected Status status_ = null;
	
	private int auxStatus_ = 0;	
	private int mainStatus_ = 0;	
	private double auxScore_ = 0;	
	private double mainScore_ = 0;	
	
	protected double[] lambda_;
	
	///** The weights for the SNPs (relevant only if snpWeightingDelta != 0) */
	//weights_ = null;
;
	/** Used for the weighted version */	
	protected UpperSymmDenseMatrix matToDecomposeMTJ_ = null;

		
	// ============================================================================
	// PUBLIC METHODS
	
    /** Initialization */
    static {
    	if (Settings.useImhof_){
    		System.loadLibrary("jni-gsl");
    	}
    }

    
	/** Constructor */
    // TODO delete? The initialization below is necessary...
    public AnalyticVegas(){
    	
    }
 
   
    public AnalyticVegas(ArrayList<Double> snpScores, DenseMatrix ld) {		
		super(snpScores, ld);	
    }
    public AnalyticVegas(ArrayList<Double> snpScores, UpperSymmDenseMatrix ld) {		
		super(snpScores, ld);	
    }
/*
		matToDecomposeMTJ_ = new UpperSymmDenseMatrix(snpScores.size());
		double[] emptyWeights = new double[snpScores.size()];
		for (int i = 0 ; i < emptyWeights.length ; i++){
			emptyWeights[i]=1;
		}
		initialize(emptyWeights);
	}
*/    
	public AnalyticVegas(ArrayList<Double> snpScores, DenseMatrix ld, double[] snpWeights) {
		
		super(snpScores, ld);			
		setWeights(snpWeights);
	}
	
	// ----------------------------------------------------------------------------

	/** Compute genotypic variances for all snps; put into snpVar_ */
//	public void initialize(double[] weights) {
		
	//setSnpWeights(weights);
//		constructMatToDecomposeMTJ();		
	//}

	/** Compute genotypic variances for all snps; put into snpVar_ */
	//TODO:  take snpWeighting outside
	/*public void initialize(double delta) {
		
		computeWeightsFromVariances(delta);		
		constructMatToDecomposeMTJ();		
	}*/

	
	// ----------------------------------------------------------------------------
	/**  check if snpWeights have been set, if not set to default value */
	public void checkWeights(){
		
		if (weights_==null){
			weights_ = new double[snpScores_.size()];
			for (int i = 0 ; i < weights_.length ; i++){
				weights_[i]=1;
			}
		}
	}
	/** Compute the score for this snp set */
	public boolean computeScore() {
		checkWeights();
		constructMatToDecomposeMTJ();	
		// Initialize
		reinitialize();				
		// Compute the test statistic for the real p-values and the phenotype permutations if available
		
		computeTestStatisticReal();
		// Compute VEGAS gene p-values analytically		
		computeAnalyticVegasPvalues();	
		
		return status_ != Status.DAVIES_FAIL_FAREBROTHER_FAIL;
	}

		
	// ----------------------------------------------------------------------------

	/** Get string representation of result that will be written to the output file */
	  public String getResultsAsString() {
		  
		  String line = "\t" + numSnps_ + "\t" + Utils.toStringScientific10(geneScore_) + "\t" + getStatusString();		 
	        return line;
	    }

	    /** Get a header line corresponding to getResultsAsString() */
	    public String getResultsAsStringHeader() {

	        //String header = "\tpvalue\tnumSnps";
	    	String header = "\tnumSnps\tpvalue\tStatus";
	  //      if (Settings.writeDetailedOutput_)
	   //         header += "\tavgSnpCorrelation\tstatus";
	        return header;
	    }
	
	// ----------------------------------------------------------------------------

	/** Get output to be printed on console after computing score for a gene */
	public String getConsoleOutput() { 
		
		if (Settings.writeDetailedOutput_) {
			return getStatusString();
		
		}
		else return "";
	}
	

	
	// ----------------------------------------------------------------------------

	/** Get output to be printed for genes where no score could be computed (error) */
	public String getNoScoreOutput() {
		
		String str = Utils.toStringScientific10(geneScore_) + "\t" + getStatusString();		
		return str; 
	}
	
	
	/** Get a header line corresponding to getNoScoreOutput() */
	public String getNoScoreOutputHeader() { 
	
		return "Score\tStatus";
	}


	// ============================================================================
	// PRIVATE METHODS

	/** Reinitialize results before the next run */
	protected void reinitialize() {

		geneScore_ = -1;
		testStatisticReal_ = -1;
		
		status_ = Status.NOT_RUN;
		mainScore_ = Double.NaN;
		mainStatus_ = -1;
		auxScore_ = Double.NaN;
		auxStatus_ = -1;
		}
	
	// ----------------------------------------------------------------------------

	/** Compute VEGAS gene p-values analytically */
	protected void computeAnalyticVegasPvalues() {
				
        // Initialize at NA in case we fail on the way
        geneScore_ = Double.NaN;

		// Get the positive eigenvalues of correlation matrix accounting for 99.99% of the variance
		lambda_ = computeLambda();		
		// If there is only one lambda (i.e., few perfectly correlated snps),
		// special treatment.
		// before take the mean p-value
		// TODO: check if code below makes any difference (Code before interfered with encapsulation)
		if (lambda_.length < 1)
			throw new RuntimeException("zero eigen values: should be captured beforehand");
        	
		if (lambda_.length == 1){
			//ChiSquaredDistribution chiSquareDist1 = new ChiSquaredDistribution(1);			
			double chiStat = testStatisticReal_/lambda_[0];
			double result = DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(chiStat);
			//double result = 1- chiSquareDist1.cumulativeProbability(chiStat);
			
				
			geneScore_ = result;
			status_ = Status.NOT_RUN;
			return;
		} 
		WeightedChisquareAlgorithm mainAlgo = null;
		if (useOnlyOneAlgorithm()){
			if (Settings.useDavies_){
			mainAlgo = new Davies(lambda_);			
			}
			else if (Settings.useFarebrother_){
			mainAlgo = new Farebrother(lambda_);			
			}
			else if (Settings.useImhof_){
				mainAlgo = new Imhof(lambda_);
			}
		mainScore_ = mainAlgo.probQsupx(testStatisticReal_);
		mainStatus_ = mainAlgo.getIfault();
		}
		else {
			mainAlgo = new Davies(lambda_);
			mainScore_ = mainAlgo.probQsupx(testStatisticReal_);
			mainStatus_ = mainAlgo.getIfault();
			if (mainScore_ < 1e-12 | mainStatus_ != 0){
				WeightedChisquareAlgorithm auxAlgo = new Farebrother(lambda_);
				auxScore_ = auxAlgo.probQsupx(testStatisticReal_);
				auxStatus_ = auxAlgo.getIfault();
								
			}
		}		
		setStatusAndScore();		
	}
	private boolean useOnlyOneAlgorithm(){
		return (Settings.useDavies_ || Settings.useFarebrother_ || Settings.useImhof_);
	}
		private void setStatusAndScore(){
			if (useOnlyOneAlgorithm()){
				if (mainStatus_!= 0){
					status_=Status.FAIL;
					geneScore_ = Double.NaN;
				}
				else{
					status_=Status.CONVERGED;
					geneScore_=mainScore_;
				}
				
			}
			else{
				if (mainStatus_ == 0 && mainScore_ > 1e-12){
					status_=Status.DAVIES_SUCCESS;
					geneScore_=mainScore_;
				}
				if (mainStatus_ == 0 && mainScore_ <  1e-12 && auxStatus_ !=0 ){
					status_=Status.DAVIES_LOWPRECISION_FAREBROTHER_FAIL;
					geneScore_=1e-12;
				}
				if (mainStatus_ != 0 && auxStatus_ !=0 ){
					status_=Status.DAVIES_FAIL_FAREBROTHER_FAIL;
					geneScore_=Double.NaN;
				}
				if (mainStatus_ != 0 && auxStatus_ == 0 ){
					status_=Status.DAVIES_FAIL_FAREBROTHER_SUCCESS;
					geneScore_ = auxScore_;
				}
				if (mainStatus_ == 0 && auxStatus_ == 0 && mainScore_ < 1e-12){
					status_=Status.DAVIES_LOWPRECISION_FAREBROTHER_SUCCESS;
					geneScore_ = auxScore_;
				}
				if (mainStatus_ != 0 && auxStatus_ != 0){
					status_=Status.DAVIES_FAIL_FAREBROTHER_FAIL;
					geneScore_ = 1;
				}
			}
		}	
	
	
	
	// ----------------------------------------------------------------------------
		//TODO:  take snpWeighting outside
	/** Compute genotypic variances for all snps; put into snpVar_ */
	/*private void computeWeightsFromVariances(double delta){

		int counter=0;
		Iterator<Double> iter = snpScores_.iterator();
		snpWeights_ = new double[snpScores_.size()];
		while (iter.hasNext()) {
			double maf = iter.next().getMaf();
			assert(maf <= 1.0 && maf >= 0.0);
		
			//snpWeights_[counter] = Math.max(tmp*(1-tmp), 0.001);
			snpWeights_[counter] = Math.pow(maf*(1.0-maf), delta);

		
		
			counter++;
		}
	}
	*/

	// ----------------------------------------------------------------------------
		public void constructMatToDecomposeMTJNew() {
			matToDecomposeMTJ_ = new UpperSymmDenseMatrix(covariance_);
		}
//TODO in release version: remove this and replace with much simpler above.		
	public void constructMatToDecomposeMTJ() {
		
		// constructWeightMat
		DoubleFactory1D factory1d = DoubleFactory1D.dense;
		DoubleFactory2D factory2d = DoubleFactory2D.dense;
		DoubleMatrix2D weightMat = factory2d.diagonal(factory1d.make(weights_));	
		DenseMatrix weightMatMTJ = new DenseMatrix(weightMat.toArray());
		//TODO: change back: just for debugging.
		DenseMatrix tmp = new DenseMatrix(covariance_);
		
		UpperSymmDenseMatrix mtjMat = new UpperSymmDenseMatrix(MTJConvenienceMethods.regularizeMat(tmp,0.0));
		SymmDenseEVD evd = null;
		try {
			evd = SymmDenseEVD.factorize(mtjMat);
		} catch (NotConvergedException e) {
			Main.println("Error: Eigendecomposition failed to converge");
		}
		DenseMatrix eigenVects = evd.getEigenvectors();	
		double[] eigenVals = evd.getEigenvalues();
		
		DenseMatrix eigenValsSqrtMat = new DenseMatrix(eigenVals.length, eigenVals.length);		
		for (int i = 0;i < eigenVals.length; ++i){			
			eigenValsSqrtMat.set(i,i, Math.sqrt(Math.max(eigenVals[i],0)));		
		}
		
		DenseMatrix sqrtMat = new DenseMatrix(eigenVals.length, eigenVals.length);		
		DenseMatrix interim = new DenseMatrix(eigenVals.length, eigenVals.length);
		eigenVects.mult(eigenValsSqrtMat,sqrtMat);
				
		interim=MTJConvenienceMethods.doGammaTLambdaGammaMTJ(sqrtMat,weightMatMTJ);
		matToDecomposeMTJ_ = new UpperSymmDenseMatrix(interim);
		
	}

	
	// ----------------------------------------------------------------------------
//TODO: Hack to use zscores to calculate chi-square stat if availabe 
	//TODO:: problem here::do we get snpScores as zscores or pvals or chisquare
	/** Compute test statistics with weights */
	protected double computeTestStatisticRealSubclass() {
		double myTest=0;
		if (Settings.withZScore_){
			//double[] zscores = new double[numSnps_];
			//for (int i=0; i<numSnps_; i++)
				//zscores[i] = geneSnps_.get(i).getZscore();						
			for (int i=0; i<numSnps_; i++)
				myTest += snpScores_.get(i)*snpScores_.get(i)*weights_[i];	
		}
		else{
			double[] chi2 = new double[numSnps_];
			for (int i=0; i<numSnps_; i++)
				chi2[i] = snpScores_.get(i);
			
			for (int i=0; i<numSnps_; i++)
				myTest += chi2[i]*weights_[i];	
		}
		return myTest;				
	}




	// ----------------------------------------------------------------------------

	/** 
	 * Get positive eigenvalues of correlation matrix accounting for 99.99% of the variance.
	 * Computes lambdas from matToDecompose_ instead of covariance_ 
	 */
	public double[] computeLambda() {

		// Apache eigendecomposition
		// Amazingly, this doesn't work, fails to converge for many genes while R and Colt never have problems
		//EigenDecomposition eigen = new EigenDecomposition(covariance_);
		//double[] eigenValues = eigen.getRealEigenvalues();

		// Colt eigendecomposition
		//DenseDoubleEigenvalueDecomposition eigen = new DenseDoubleEigenvalueDecomposition(matToDecompose_);

		// MTJ eigendecomposition
		SymmDenseEVD evd = null;
		try {
			evd = SymmDenseEVD.factorize(matToDecomposeMTJ_);
		} catch (NotConvergedException e) {
			Main.println("Error: Eigendecomposition failed to converge");
		}
		double[] ev = evd.getEigenvalues();

		// Sort so that leading values are first
		ArrayList<Double> eigenValues = new ArrayList<Double>(ev.length);
		// I'm adding them in reverse order because that's how it seems they're ordered by EigenDecomposition
		// TODO sorted correctly?
		for (int i=0; i<ev.length; i++)
			eigenValues.add(ev[i]);
		    //eigenValues.add(ev.get(N-1-i));				
		Collections.sort(eigenValues, Collections.reverseOrder());

		// Sum of positive eigenvalues
		double sum = 0;
		for (Double x : eigenValues)
			if (x > 0)
				sum += x;

		double cutoff = sum / (Settings.eigenValueFractionCut_);
		// The index where the eigenvalues are below the cutoff		
		int k = 0;
		for (Double x : eigenValues) {
			sum -= x;
			k++;
			if (sum < cutoff)
				break;			
		}

		// Construct lambda, the i first eigenvalues
		assert k >= 0;
		double[] lambda = new double[k];

		for (int i=0; i<k; i++)
			lambda[i] = eigenValues.get(i);

		return lambda;
	}

	
	// ----------------------------------------------------------------------------

	/** Reinitialize results before the next run */
	protected String getStatusString() {
		
		return status_.name();
	}
		
	
	// ============================================================================
	// SETTERS AND GETTERS
//	public void setSnpWeights(double[] snpWeights){weights_= snpWeights;}
//	public double[] getSnpWeights(){return weights_;}
	//public int getErrorCode() { return imhofStatus_; }
	//public double[] getLambda(){return lambda_;}
	public void setCovariance(DenseMatrix covariance) { covariance_ = new UpperSymmDenseMatrix(covariance); }
	public UpperSymmDenseMatrix getMatToDecompose(){return matToDecomposeMTJ_;}
	//public DenseMatrix getMatrixToDecompose(){return matrixToDecompose_;}

	
}
