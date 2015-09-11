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
