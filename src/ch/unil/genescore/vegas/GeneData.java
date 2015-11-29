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

import java.io.File;
import java.util.ArrayList;

import ch.unil.genescore.main.Pascal;
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
		if (Pascal.set.withZScore_){
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
		
	 public void writeCovMatToFile(File file){
	    	
    	WritingMethods.writeMTJ(cov_, file);
	 }
	 
	 
	 public void writeGeneSnpsToFile(File file){
	    	
		 ArrayList<String> valStrings = new ArrayList<String>(); 
	    	 
		 for (int i=0 ; i < snpScores_.size() ; i++){
			 String currentString = "";
			 currentString = String.valueOf(snpScores_.get(i));	    		
			 valStrings.add(currentString);
		 }
		 String header="scores";
		 WritingMethods.writeSnpPosWithValToFile(snpList_, valStrings, header, file);
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
