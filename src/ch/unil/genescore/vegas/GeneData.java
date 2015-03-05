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

import ch.unil.genescore.main.Settings;
import no.uib.cipr.matrix.DenseMatrix;

public class GeneData implements GeneDataInterface {
	protected DenseMatrix cov_; 
	protected ArrayList<Double> snpScores_;
	protected ArrayList<Snp> snpList_;
	
	public GeneData(ArrayList<Snp> snpList){
		snpScores_ = new ArrayList<Double>();
		snpList_=snpList;		
	}
	
	
	protected void setSnpScores(){
		snpScores_ = new ArrayList<Double>();
		if (Settings.withZScore_){
			for (Snp snp : snpList_){
				snpScores_.add(snp.getZscore());				
			}
		}
		else{
			for (Snp snp : snpList_){
				
				snpScores_.add(snp.getChi2Stat());			
			}
		}		
	}
	
	public DenseMatrix getCorr(){return cov_;}
	public ArrayList<Double> getScores(){return snpScores_;}
		
	 public void writeCovMatToFile(String fileName){
	    	
		 //ArrayList<String> valStrings = new ArrayList<String>(); 
	    	 
	    //for (int i=0 ; i < snpList_.size() ; i++){
	    //	String currentString = "";
	   // 	currentString += String.valueOf(snpScores_);	    		
	  //  	valStrings.add(currentString);
	//    }
	    String header="scores";
	  //  WritingMethods.writeSnpPosWithValToFile(snpList_, valStrings,header, fileName,Settings.writeGenewiseSnpFiles_);
	    //String corfileName = "corMat_" + fileName;
    	WritingMethods.writeMTJ(cov_, fileName,Settings.writeGenewiseSnpFiles_);
	 }
	 
	 
	 public void writeGeneSnpsToFile(String fileName){
	    	
		 ArrayList<String> valStrings = new ArrayList<String>(); 
	    	 
	    for (int i=0 ; i < snpScores_.size() ; i++){
	    	String currentString = "";
	    	currentString = String.valueOf(snpScores_.get(i));	    		
	    	valStrings.add(currentString);
	    }
	    String header="scores";
	    WritingMethods.writeSnpPosWithValToFile(snpList_, valStrings,header, fileName,Settings.writeGenewiseSnpFiles_);
	 }
	 
	 public void processData(){
		 cov_ = LinkageDisequilibrium.computeCorrelationMatrixMTJ(snpList_);
		setSnpScores();
	 }
	 public double[] returnWeights(boolean dummy){
		 int n= snpScores_.size();
		 double[] weights = new double[n];
		 for (int i=0; i<n; i++)
			 	weights[i]=1;
		 return weights;
	 }

	@Override
	public double[] getWeights() {
		double[] d = new double[snpScores_.size()];
		for(int i=0; i<snpScores_.size();i++){
			d[i]=1;
		}
		return d;
	}	 
}
