/*******************************************************************************
 * Copyright (c) 2015 IBM Corporation and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     IBM Corporation - initial API and implementation
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
}
	
