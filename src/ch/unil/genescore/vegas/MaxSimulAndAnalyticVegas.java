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
import ch.unil.genescore.vegas.VegasSimulationNew;;

public class MaxSimulAndAnalyticVegas extends AnalyticVegas {

	AnalyticVegas myAnalyticVegas_ = null;
	
	
	
	public MaxSimulAndAnalyticVegas(){

	}
	
	public MaxSimulAndAnalyticVegas(ArrayList<Double> snpScores, DenseMatrix ld) {
		super(snpScores, ld);
		
	}
	
	public MaxSimulAndAnalyticVegas(ArrayList<Double> snpScores, DenseMatrix ld, double[] weights) {			
		super(snpScores, ld);					
		weights_ = weights;						
	}
	
	@Override
	public boolean computeScore(){
	//	if (weights_==null){
	//		setDefaultWeights();
		//}
		VegasSimulationNew simulVegas = new VegasSimulationNew(snpScores_,covariance_,(1E5));
		simulVegas.computeScore();
		if (simulVegas.checkSimulFail()){
			simulVegas= null;		
			MaxVegas maxVegas= new MaxVegas(snpScores_,covariance_, weights_);
			maxVegas.computeScore();
			if(maxVegas.checkTooManySnps()){
				
				simulVegas = new VegasSimulationNew(snpScores_,covariance_,(1E7));
				simulVegas.computeScore();
				//myAnalyticVegas_ = simulVegas;
			}
			else {
				myAnalyticVegas_ = maxVegas;
			}
		}
		else{
			myAnalyticVegas_ = simulVegas;
		}
		copyResultsOver();
		//MaxVegasWithoutPruning analyticVegas= new MaxVegasWithoutPruning(snpScores_,covariance_, weights_);
		return true;
	}
	private void copyResultsOver(){
		/**copies results from lower-level myAnalyticVegas_ to this*/
		geneScore_=myAnalyticVegas_.geneScore_;	
		status_=myAnalyticVegas_.status_;	
		tolerance_=myAnalyticVegas_.tolerance_;
		positiveDefiniteCorrection_=myAnalyticVegas_.positiveDefiniteCorrection_;
		numSnps_=myAnalyticVegas_.numSnps_;
		
	}
	
	private void setDefaultWeights(){
		double[] weights = new double[snpScores_.size()];
		for (int i=0; i<snpScores_.size(); i++){
			weights[i]=1;					
		}		
		weights_ = weights;			
	}

	@Override
	public String getResultsAsString() {
		 
		return myAnalyticVegas_.getResultsAsString();
	}
	@Override
	public String getResultsAsStringHeader() {
		return myAnalyticVegas_.getResultsAsStringHeader();		
	}
	
	@Override
	public String getTypeString() {
		return "max";		
	}
}
	
