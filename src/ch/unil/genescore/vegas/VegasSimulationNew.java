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

	static Random r = new Random();
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
