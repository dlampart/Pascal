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

import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.LowerTriangDenseMatrix;
import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import no.uib.cipr.matrix.UpperTriangDenseMatrix;
import no.uib.cipr.matrix.Vector.Norm;
import ch.unil.genescore.main.Settings;

import java.util.*;

public class VegasSimulationNew extends AnalyticVegas {

	static Random r = new Random(Settings.randomSeed_);
	double maxStat_;
	boolean tooManySimuls_ = false;
	double simulRuns_ = 0;
	double maxAllowedSimulRuns_ = 0;
	public VegasSimulationNew(ArrayList<Double> snpScores, DenseMatrix ld, double maxAllowedSimulRuns) {			
		super(snpScores, ld);			
		maxAllowedSimulRuns_=maxAllowedSimulRuns;
	
	}
	
	public VegasSimulationNew(ArrayList<Double> snpScores, UpperSymmDenseMatrix ld, double maxAllowedSimulRuns) {			
		super(snpScores, ld);			
		maxAllowedSimulRuns_=maxAllowedSimulRuns;
	
	}
	@Override
	protected double computeTestStatisticRealSubclass() {
		maxStat_=0;
		if (Settings.withZScore_){
			for (int i=0; i<snpScores_.size(); i++){				
				maxStat_=Math.max(Math.abs(snpScores_.get(i)),maxStat_);
			}
		}
		else{
			for (int i=0; i<snpScores_.size(); i++){
				maxStat_=Math.max(Math.sqrt(snpScores_.get(i)),maxStat_);
			}
		
		}
		
		return 0;
	}
	
	protected double computeTestStatisticSimul(DenseVector w ) {
		double maxStat_=0;
		if (Settings.withZScore_){
			for (int i=0; i<snpScores_.size(); i++){
				maxStat_=Math.max(Math.abs(snpScores_.get(i)),maxStat_);
			}
		}
		else{
			for (int i=0; i<snpScores_.size(); i++){
				maxStat_=Math.max(Math.sqrt(snpScores_.get(i)),maxStat_);
			}
		
		}
		return 0;
	}

	@Override
	public boolean computeScore() {
		computeTestStatisticRealSubclass();
		matToDecomposeMTJ_ = MTJConvenienceMethods.regularizeMat(covariance_, 0.0001);
		
		DenseCholesky chol = DenseCholesky.factorize(matToDecomposeMTJ_);
		UpperTriangDenseMatrix U = chol.getU();
		
		DenseVector normVect = new DenseVector(matToDecomposeMTJ_.numColumns());
		
		DenseVector corVect = new DenseVector(matToDecomposeMTJ_.numColumns());
		double count=0;
		simulRuns_=100;
		while (count<100){			
			count=0;
			simulRuns_=simulRuns_*10;
			if (simulRuns_ > maxAllowedSimulRuns_){
				tooManySimuls_=true;
				count=0;
				break;
			}
		//	System.out.println("nr of Simuls " + simulRuns_);
			for (int i=0; i<simulRuns_; i++){			
				simulateNormals(normVect);				
				U.transMult(normVect, corVect);
				double maxSimul = corVect.norm(Norm.Infinity);		
				if (maxSimul> maxStat_){
					count++;
				} 
			}
		}
		geneScore_ = count/simulRuns_;
		return true;
	}
	
	private void simulateNormals(DenseVector normVect) {
		for (int i=0; i<normVect.size();i++){
			normVect.set(i,r.nextGaussian());			
		}
		
	}
	/** Reinitialize results before the next run */
	protected String getStatusString() {
		if (tooManySimuls_){
			return "SimulFail";
		}
		return "Simul";
	}
	public boolean checkSimulFail() {
		return tooManySimuls_;			
	}
	
}
