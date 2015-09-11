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

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.MaxVegas.Status;

import com.sun.jna.Native;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

public class MaxVegasWithoutPruning extends AnalyticVegas {
	
	protected enum Status {CONVERGED, FAIL, NOT_RUN, TOO_MANY_SNPS, NOT_POSDEF};
	Status status_ = null;
//	static MvnPack INSTANCE = (MvnPack) Native.loadLibrary("mvtpack", MvnPack.class);;
	MvnPackDirectMapping INSTANCEDirect = new MvnPackDirectMapping();
	protected double[] upper_ = null;
	
	protected UpperSymmDenseMatrix correl_ = null;
	
	public MaxVegasWithoutPruning(ArrayList<Double> snpScores, DenseMatrix ld) {			
		super(snpScores, ld);			
		System.setProperty("jna.debug_load", "true");
		
		//INSTANCE = (MvnPack) Native.loadLibrary("/Users/dlampart/Documents/workspace/genescore/fortranlibs/libmvtpack.dlyb", MvnPack.class);
		//INSTANCE = (MvnPack) Native.loadLibrary("fortranlibs/libmvtpack.dlyb", MvnPack.class);
		
		correl_= new UpperSymmDenseMatrix(ld);
		String myJna = System.getProperty("jna.library.path");
		System.out.println(myJna);
	    //System.setProperty("jna.library.path", <path to your library>);
		//INSTANCE = (MvnPack) Native.loadLibrary("mvtpack", MvnPack.class);
		System.out.println(myJna);
		
				
	}
	@Override
	public boolean computeScore() {
		
		
		
		calculateWeightedMaximum();
		stretchToVarOne();
		if (upper_.length==0){geneScore_=1;status_=Status.TOO_MANY_SNPS;return 	status_ != Status.FAIL;}
		runGenescoreCalculation();
		return status_ != Status.FAIL;
	}

	public void calculateWeightedMaximum(){
		double maxStat=0;
		for (int i=0; i<snpScores_.size(); i++){
			maxStat=Math.max(Math.abs(snpScores_.get(i)),maxStat);
		}
		upper_ = new double[snpScores_.size()];
		for (int i=0; i<snpScores_.size(); i++){
			upper_[i] = maxStat;
		}				
	}
	public void stretchToVarOne(){
		
		 double[] variances=MTJConvenienceMethods.getDiagDar(correl_);
		
			for (int i=0; i<numSnps_; i++)							
				upper_[i]=upper_[i]/(Math.sqrt(variances[i]));				
			correl_= MTJConvenienceMethods.scaleToDiagOne(correl_);
	}
	public void runGenescoreCalculation(){
		
		calculateGenescore();
		if (status_==Status.NOT_POSDEF){
			correl_= MTJConvenienceMethods.regularizeMat(correl_, 0.0001);
			correl_= MTJConvenienceMethods.scaleToDiagOne(correl_);
			calculateGenescore();
		}
		
	}
	private void calculateGenescore(){	
		int n = upper_.length;
		int[] infin= new int[n];
		double[] delta= new double[n];
		double[] lower= new double[n];
		for (int i=0; i< n; i++){
			infin[i]=2;
			delta[i]=0;
			lower[i]=(-1)*upper_[i];
		}
		
		int maxpts = 100000;
		double abseps=1e-3;
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
			System.out.println("start");
	//		long t0 = System.currentTimeMillis();
	//		INSTANCE.mvtdst_(new IntByReference(n), new IntByReference(0), lower, upper_, infin, correlflat, delta, maxpts_ref, abseps_ref, releps_ref, error, value, inform);
	//		long t1 = System.currentTimeMillis();
	//		long timeDiff=t1-t0;
	//		System.out.println(timeDiff);
	//		System.out.println("sec time");
		
		long t0 = System.currentTimeMillis();
		System.out.println("number of variables"+ n);
		MvnPackDirectMapping.mvtdst_(new IntByReference(n), new IntByReference(0), lower, upper_, infin, correlflat, delta, maxpts_ref, abseps_ref, releps_ref, error, value, inform);
		long t1 = System.currentTimeMillis();
		long timeDiff=t1-t0;
		System.out.println(timeDiff);
		geneScore_=(1-value.getValue());
		double err = error.getValue();
		
		status_ = evaluateMvtdstStatus(inform.getValue());
		}
		else{
			geneScore_=1;
			status_=Status.TOO_MANY_SNPS;			
		}
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
}
