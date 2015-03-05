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
import java.util.HashMap;

public class SnpWeightPairs {

	//Snp snp1;
	private ArrayList<Snp> snps_=null;
	private ArrayList<Double> weights_ = null;
	
	public SnpWeightPairs(){
		
	}
	
	
	public SnpWeightPairs(ArrayList<Snp> snps,ArrayList<Double> weights){
		
		snps_ = snps;
		weights_ = weights;
	}
	public void fuseSnpWeightPairs(SnpWeightPairs secondSnpWeightPairs){
		HashMap<String, Snp> snpHash = new HashMap<String, Snp>();
		HashMap<String, Double> weightHash = new HashMap<String, Double>();
		String currentId=null;
		for (int i=0;i<snps_.size(); i++){
			currentId=snps_.get(i).id_;
			if (snpHash.containsKey(currentId)){
				throw new RuntimeException("error: snpWeightPairs should only contain unique entries");
			}
			snpHash.put(currentId, snps_.get(i));
			weightHash.put(currentId, weights_.get(i));			
		}		
		for (int i=0;i<secondSnpWeightPairs.getSnps().size();i++){
			currentId=secondSnpWeightPairs.getSnps().get(i).id_;
			double otherWeight = secondSnpWeightPairs.getWeights().get(i);
			if (snpHash.containsKey(currentId)){
				double currentWeight = weightHash.get(currentId);				
				double maxWeight = Math.max(currentWeight, otherWeight);
				weightHash.put(currentId, maxWeight);
			}
			else {				
				snpHash.put(currentId, secondSnpWeightPairs.getSnps().get(i));
				weightHash.put(currentId, otherWeight);
			}
			snpHash.put(currentId, snps_.get(i));
			weightHash.put(currentId, weights_.get(i));	
		}		
		ArrayList<Snp> snps = new ArrayList<Snp>();
		ArrayList<Double> weights  = new ArrayList<Double>();
		for (String id : snpHash.keySet()){
			snps.add(snpHash.get(id));
			weights.add(weightHash.get(id));			
		}
		snps_ = snps;
		weights_ = weights;
	}
	
	public ArrayList<Snp> getSnps(){return snps_;}
	public ArrayList<Double> getWeights(){return weights_;}
	public double[] getWeightsNotList(){
		int n=weights_.size();
		double[] outAr  = new double[n];
		for (int i=0 ; i < n ; i++)
			outAr[i]=weights_.get(i);		
		return outAr;
	}	
}
