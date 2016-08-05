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

import com.sun.jna.Native;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.AnalyticVegas.Status;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
/**
 * used to calculate max-Statistic
 * 
 * */
public class MaxVegas extends AnalyticVegas {
	
	
	protected enum Status {CONVERGED, FAIL, NOT_RUN, TOO_MANY_SNPS, NOT_POSDEF, EFF_NR_OF_TEST_ADJ};
	Status status_ = null;
	
	static MvnPack INSTANCE = (MvnPack) Native.loadLibrary("mvtpack", MvnPack.class);
	//
	//MvnPackDirectMapping INSTANCEDirect = new MvnPackDirectMapping();
	protected double[] upper_ = null;
	protected int[] pruningIndices_ = null;
	protected UpperSymmDenseMatrix correl_ = null;
	protected double pruningCutoff_;	
	//private double[] weights_;
	protected double[] variances_;
	
	/**
	 * used to calculate max-Statistic
	 * 
	 * inputs: 
	 * snpScores: weighted snpScores;
	 * ld: weighted covariance matrix.
	 * weigths: weights used
	 * 
	 * uses variables correl_ and upper_ to process matrices and score vector
	 * 
	 * analysis steps:
	 * 
	 * this is all way too complex because I made a reasoning mistake: actually, you could just take in the weighted values and matrices. take the max of weigthted z-scores. replace these
	 * z-scores with max-zscore. scale cov-mat to diag one and the max-zscore such that cov-mat has diag=I. pruning would then have to be done outside of maxVegas though.
	 * 
	 * 
	 * 1) deweight() : gets correl_ and upper_ into the state as if no weighting had happened. (stupid but has just grown that way) 
	 * 2) stretchToVarOne(): scales only correl_ and npt upper_ such that correl_ has diag 1. NOTE: badly imputed snp-scores don't get upweighted in the process.
	 * need to filter badly imputed snps to avoid instabilities.
	 * 3) prune(): prunes highly correlated snps: Implementation problem: want to keep highest weight and its actual snp-value (strategy:
	 * if more than one snp have highest value just pick first.)
	 * 4) calculateWeightedMaximum(): stretch scores according to weights:
	 *  (x_1_weight1, x_2_weight2,..) 
	 *  4.1)stretch according to variances 
	 * take maximum:
	 *  max(x_1_weight1, x_2_weight2,..)
	 * reverse stretch on same maximum  for each dimension separately according to
	 * weights. To bring every dimension again down to variance 1.
	 * 
	 * NOTES:::
	 * 
	 * deweight() UNNECCESSARY
	 * because we scale to one. weights_ can be chosen arbitrarily. E.G. the weighting can include downweigts on badly imputed snps. but this has to be calculated
	 * beforehand and be given to this function
	 * */
	
		public MaxVegas(){
			super();
			pruningCutoff_=Settings.maxPruningCutoff_;		
			String myJna = System.getProperty("jna.library.path");
		}
		public MaxVegas(ArrayList<Double> snpScores, DenseMatrix ld, double[] weights) {			
		super(snpScores, ld);			
		pruningCutoff_=Settings.maxPruningCutoff_;		
		weights_ = weights;
		System.setProperty("jna.debug_load", "true");						
				
	}
		
		public MaxVegas(ArrayList<Double> snpScores, UpperSymmDenseMatrix ld, double[] weights) {			
			super(snpScores, ld);			
			pruningCutoff_=Settings.maxPruningCutoff_;		
			weights_ = weights;
			System.setProperty("jna.library.path", "lib/fortranlibs/");
			
			
			
			
					
		}	
		
	public MaxVegas(ArrayList<Double> snpScores, DenseMatrix ld, double pruningCutoff) {			
		super(snpScores, ld);
		pruningCutoff_=pruningCutoff;		
		System.setProperty("jna.debug_load", "true");	
		//INSTANCE = (MvnPack) Native.loadLibrary("/Users/dlampart/Documents/workspace/genescore/fortranlibs/libmvtpack.dlyb", MvnPack.class);
		//INSTANCE = (MvnPack) Native.loadLibrary("fortranlibs/libmvtpack.dlyb", MvnPack.class);
		INSTANCE = (MvnPack) Native.loadLibrary("fortranlibs/mvtpack", MvnPack.class);
				
	}
	
	/**finds pruning indices i.e. which indices to keep in pruning step; additionally modifies upper_ such that upper_ values are replaced with maximums 
	 * among snps that are correlated such that it leads to pruning */
	private void findPruningIndices(){
		
		int n=correl_.numColumns();		
		boolean[] indices = new boolean[n];		
		for (int i=0;i<n ; i++){
			//check if i has already been pruned:
			if (indices[i]){
				continue;
			}
			for (int j=(i+1);j<n; j++){
				if(Math.abs(correl_.get(i,j)) > pruningCutoff_){
					indices[j]=true;
					if(weights_[j]>weights_[i]){
						weights_[i]=weights_[j];
						upper_[i]=upper_[j];
					}
				}
			}			
		}
		int counter=0;
		for (boolean index : indices){
			if (!index)
				counter++;
		}
		pruningIndices_ = new int[counter];		
		int counter2=0;
		
		for (int i=0; i< n; i++){			
			if (!indices[i]){				
				pruningIndices_[counter2]=i;
				counter2++;
			}
		}
	}
	private void pruneOnIndices(){
		correl_=MTJConvenienceMethods.getSubDenseMatrixSymmMTJ( correl_, pruningIndices_);
		upper_=ConvenienceMethods.subIndexing(upper_,pruningIndices_);
		weights_=ConvenienceMethods.subIndexing(weights_,pruningIndices_);
	}
	
	//TODO: implement a pruning function.
	public void prune(){
		findPruningIndices();
		pruneOnIndices();
	}	
		
	
	

	private double[] deweightSnpScores(){
		
		double[] upper  = new double[numSnps_];
		if (!Settings.withZScore_){
			for (int i=0; i< numSnps_; i++)
				upper[i]=Math.sqrt(snpScores_.get(i)/weights_[i]);
		}
		else{	//throw new RuntimeException("with chisquare stat not implemented atm");}				
			for (int i=0; i< numSnps_; i++)
			upper[i]=Math.abs(snpScores_.get(i)/weights_[i]);
		}
		return upper;
	}
	public void deweight(){
		correl_ = MTJConvenienceMethods.deScale(covariance_, weights_);
		upper_ = deweightSnpScores();
		
	}
	public void stretchToVarOne(){
		
		 variances_=MTJConvenienceMethods.getDiagDar(correl_);
		
			for (int i=0; i<numSnps_; i++)
			
				
			upper_[i]=upper_[i]/(Math.sqrt(variances_[i]));				
			correl_= MTJConvenienceMethods.scaleToDiagOne(correl_);
	}
	
	public void calculateWeightedMaximum(){
		double maxStat=0;
		for (int i=0; i<upper_.length; i++)
			upper_[i]=(upper_[i]*weights_[i])*Math.sqrt(variances_[i]);
		for (int i=0; i<upper_.length; i++)
			maxStat=Math.max(maxStat, upper_[i]);
		for (int i=0; i<upper_.length; i++)
			upper_[i]=(maxStat/(weights_[i]))/Math.sqrt(variances_[i]);

	}
	
	@Override
	public boolean computeScore() {
		
		deweight();
		stretchToVarOne();
		prune();
		
		calculateWeightedMaximum();			
		if (upper_.length==0){
			geneScore_=1;
			status_=Status.TOO_MANY_SNPS;
			return 	status_ != Status.FAIL;
		}
		runGenescoreCalculation();
		return status_ != Status.FAIL;
	}
	public void runGenescoreCalculation(){
		if(Settings.regFactorMax_==0){
			calculateGenescore();
			if (status_==Status.NOT_POSDEF){
				correl_= MTJConvenienceMethods.regularizeMat(correl_, 0.0001);
				correl_= MTJConvenienceMethods.scaleToDiagOne(correl_);
				calculateGenescore();
			}
		}else{
			if(Settings.regFactorMax_<0){
				throw new RuntimeException("reg factor can't be smaller than 0.");
			}
			correl_= MTJConvenienceMethods.regularizeMat(correl_,Settings.regFactorMax_);
			correl_= MTJConvenienceMethods.scaleToDiagOne(correl_);
			calculateGenescore();
		}
	}
	
	protected void calculateGenescore(){	
		int n = upper_.length;
		int[] infin= new int[n];
		double[] delta= new double[n];
		double[] lower= new double[n];
		for (int i=0; i< n; i++){
			infin[i]=2;
			delta[i]=0;
			lower[i]=(-1)*upper_[i];
		}

		int maxpts = 1000;
		double abseps=1e-15;
		double releps=1e-3;
		DoubleByReference abseps_ref = new DoubleByReference(abseps);
		DoubleByReference releps_ref = new DoubleByReference(releps);		
		IntByReference maxpts_ref = new IntByReference(maxpts);		
		DoubleByReference error = new DoubleByReference(0);
		DoubleByReference value = new DoubleByReference(0);
		IntByReference inform = new IntByReference(0);

		double[] correlflat = new double[n*(n-1)/2];
		for( int i = 0; i < n; i++ ) {
			for( int j = 0; j < i; j++ ) {
				correlflat[(j+1) + (i-1)*i/2 - 1] = correl_.get(i, j);
			}
		}					

		if (upper_.length<1000){
			long t0 = System.currentTimeMillis();
			MvnPackDirectMapping.mvtdst_(new IntByReference(n), new IntByReference(0), lower, upper_, infin, correlflat, delta, maxpts_ref, abseps_ref, releps_ref, error, value, inform);
			long t1 = System.currentTimeMillis();
			long timeDiff=t1-t0;
			geneScore_=(1-value.getValue());

			double err = error.getValue();		
			status_ = evaluateMvtdstStatus(inform.getValue());
			///do checks for errors.
			if(Settings.externalConvergenceCheck_){
				for(int i = 0; i< 8; i++){
					if(err*3 < geneScore_ && status_==Status.CONVERGED){
						break;
					}

					t0 = System.currentTimeMillis();
					MvnPackDirectMapping.mvtdst_(new IntByReference(n), new IntByReference(0), lower, upper_, infin, correlflat, delta, maxpts_ref, abseps_ref, releps_ref, error, value, inform);
					t1 = System.currentTimeMillis();
					timeDiff=t1-t0;

					geneScore_=(1-value.getValue());		
					err = error.getValue();
					status_ = evaluateMvtdstStatus(inform.getValue());
				}
			}
			if(geneScore_==0){			
				geneScore_=computeEffectivePval();	
				status_=Status.EFF_NR_OF_TEST_ADJ;
			}
		}
		else{
			geneScore_=1;
			geneScore_=computeEffectivePval();	
			status_=Status.EFF_NR_OF_TEST_ADJ;
			//status_=Status.TOO_MANY_SNPS;			
		}
	}

	protected double  computeEffectivePval(){
		int effTest = computeEffectiveNrOfTests();
		double maxStat = 0;
	
		for (double score : snpScores_){
			if (Settings.withZScore_){
				score=Math.pow(score,2);
			}
			if (maxStat< score)
				maxStat=score;
		}
		double pval=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(maxStat);
		double pvalEff=Math.min(pval*effTest,1);
		return pvalEff;
	}
	
	private int computeEffectiveNrOfTests(){
		
		constructMatToDecomposeMTJNew();
		double[] lambdas=computeLambda();
		double nrOfSnps=snpScores_.size();
		double totSums=0;		
		int i=0;
		while(nrOfSnps*0.995 > totSums){
			totSums=totSums+lambdas[i];		
			i++;
		}		
		return(i);
		
		
	}
		

	Status evaluateMvtdstStatus(int status){
	Status myStatus = null;
	if(status==3){
		//throw new RuntimeException("mvtdst: correlation mat not positive definite");
		myStatus = Status.NOT_POSDEF;
	}
	if(status==2){
		throw new RuntimeException("mvtdst: too many (or not enough?)variables");
	}
	if(status==0){
		myStatus = Status.CONVERGED;
	}
	if(status==1){
		//TODO: change back
		myStatus = Status.CONVERGED;
		//myStatus = Status.FAIL;
	}	
	return myStatus;
		
}
	public boolean checkTooManySnps(){
		if (status_==Status.TOO_MANY_SNPS){
			return true;
		}
		else {return false;
		}		
	}
	
	@Override
	public String getResultsAsString() {
		String line = "\t" + numSnps_ + "\t" + Utils.toStringScientific10(geneScore_) + "\t" + status_.toString();		 
		return line;
	}
	

	@Override
	public String getResultsAsStringHeader() {

     	String header = "\tnumSnps\tpvalue\tStatus";
        return header;
	
	}
	public void setPruningCutoff(double pruningCutoff){
		pruningCutoff_=pruningCutoff;
	}
		
}
