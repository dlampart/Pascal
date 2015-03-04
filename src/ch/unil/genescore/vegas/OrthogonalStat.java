package ch.unil.genescore.vegas;

import java.util.ArrayList;
import java.util.Arrays;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;

// TODO include in release? delete or just hide?
public class OrthogonalStat extends GeneScoreEvaluator {
	
	int nrComponents_ = 1;
	private ArrayList<Double> snpScores_ = null;
	private UpperSymmDenseMatrix covariance_ = null;
	private DenseVector corStats_ = null;
	private DenseVector indepStats_ = null;
	private double[] eigVals_ = null; 
	protected double[] indepScaledStats_ = null;	
	private static NormalDistribution normDist_ = new NormalDistribution();
	protected double geneScore_ = 1;
	
	public OrthogonalStat(){};
	public OrthogonalStat(ArrayList<Double> snpScores, DenseMatrix ld, int nrComponents) {
	
			setSnpScores(snpScores);
			setCorStats();
			setCovariance(ld);
			orthogonalizeStat();
			scaleIndepStat();	
	}

	public OrthogonalStat(ArrayList<Double> snpScores, DenseMatrix ld){
		 this(snpScores, ld, 1);
	}
			
		/** Assign the snps */
	public void setSnpScores(ArrayList<Double> snpScores) {
		snpScores_ = snpScores;
			
	}
	public void setCovariance(DenseMatrix covariance) {
		covariance_ = new UpperSymmDenseMatrix(covariance);	
			// Set the covariance matrix			
	}
	private void setCorStats(){
		corStats_ = new DenseVector(snpScores_.size());
		int count=0;
		for (double zscore : snpScores_){			
			corStats_.set(count, zscore);
			count++;
		}
	}
	public void setNrComponents(int nrComponents){		
		nrComponents_ = nrComponents;
	}
	public double getGeneScore(){
		return geneScore_;
	}
	
	private void orthogonalizeStat() {
		UpperSymmDenseMatrix mtjMat = new UpperSymmDenseMatrix(new DenseMatrix(covariance_));		
		SymmDenseEVD evd = null;
		try {
			evd = SymmDenseEVD.factorize(mtjMat);
		} catch (NotConvergedException e) {
			Main.println("Error: Eigendecomposition failed to converge");
		}
		DenseMatrix eigenVects = evd.getEigenvectors();	
		eigVals_ = evd.getEigenvalues();
		eigenVects.transpose();
		indepStats_ = new DenseVector(eigVals_.length);
		eigenVects.mult(corStats_, indepStats_);	
		
}
	private void scaleIndepStat() {
		
		double cur = -1;
		int nrOfRelVals=0;
		double totSum=0;
		for (double eigVal : eigVals_){
			totSum=totSum+eigVal;
		}
		double relSum=0;
		double leftoverVar = 1-Settings.varExplained_;
		nrOfRelVals = eigVals_.length;		
		for (double eigVal : eigVals_){
			assert(cur <= eigVal);
			cur = eigVal;
			relSum=relSum+eigVal;			
			if (relSum/totSum < leftoverVar){
				nrOfRelVals--;				
			}
		}		
		//for (double eigVal : eigVals_){
		//	assert(cur <= eigVal);
		//	cur = eigVal;
		//	if (eigVal >= cutoff_)
		//		nrOfRelVals++;
		//}
		indepScaledStats_ = new double[nrOfRelVals];
		for (int i= 0; i <  indepScaledStats_.length; ++i){
			int revI = eigVals_.length - i - 1;
			indepScaledStats_[i] = Math.abs(indepStats_.get(revI) / Math.sqrt(eigVals_[revI]));
		}
		
		Arrays.sort(indepScaledStats_);
		 System.out.println(geneScore_);
	}
	public boolean computeScore(){
		
		//disables nrComponents_ use instead varExplained
		//TODO: take out nrComponents_
		nrComponents_=indepScaledStats_.length;
		//
		int index = indepScaledStats_.length-nrComponents_;
		if (index < 0){						
			geneScore_ = 1.0;
			return true;
		}
		//trick to get to lower tail;
		 double singleTailProb = normDist_.cumulativeProbability((indepScaledStats_[index])*(-1)*2);
		 BinomialDistribution bin = new BinomialDistribution(indepScaledStats_.length, singleTailProb);
		 geneScore_ = 1.0 - bin.cumulativeProbability(nrComponents_-1);		
		 // TODO return true if success, false if fail
		return true;
	}
	
	public String getResultsAsString() {
		
		// Tab-separated list of gene scores
		String line = "";					
				line += "\t" + Utils.toStringScientific10(geneScore_);						
		return line;
	}
	
	public String getResultsAsStringHeader() {
		
		// TODO
		return "\tpval";
	}
	
	/** Get output to be printed on console after computing score for a gene */
	// TODO
	public String getConsoleOutput() {
		
		return "";
	}

	
	/** Not used for this subclass */
	protected double computeTestStatisticReal(int k) {
		throw new RuntimeException("Should never be called");
	}
	
	
	/** Get gene score(s) */
	public double[] getScore() {
		
		double[] result = { geneScore_ };
		return result;
	}
	@Override
	public void setDataFromGeneData(GeneData Dat) {
		setCovariance(Dat.getCorr());
		setSnpScores(Dat.getScores());    	
		setCorStats();		
		orthogonalizeStat();
		scaleIndepStat();
    	
		
	}
	


}
