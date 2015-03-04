package ch.unil.genescore.vegas;

import java.util.ArrayList;
import java.util.Arrays;

import no.uib.cipr.matrix.DenseMatrix;
/** class to handle our Correlation Matrices grouped with the corresponding snps and weights*/
public class LdMat {

	
	DenseMatrix rowGenotypes_ = null;
	DenseMatrix columnGenotypes_ = null;
	private ArrayList<Snp> rowSnpList_ = null;
	private ArrayList<Snp> columnSnpList_ = null;
	private double[] rowWeights_ = null;
	private double[] columnWeights_ = null;
	
	SnpWeightPairs rowSnpWeightPairs_ = null;
	//public LdMat (ArrayList<Snp> rowSnpList){
		
	//	rowSnpList_=rowSnpList;
	//	columnSnpList_=rowSnpList;
		
	//	rowWeights_ = new double[rowSnpList_.size()];
	//	columnWeights_ = new double[rowSnpList_.size()];	
	//	Arrays.fill(rowWeights_, 1);
	//	Arrays.fill(columnWeights_, 1);		
	//	weightedCrossCorrelationMatrix_=LinkageDisequilibrium.computeCorrelationMatrixMTJ(rowSnpList_);
	//}
	
	/**
	public LdMat (ArrayList<Snp> rowSnpList, ArrayList<Snp> columnSnpList){
		
		rowSnpList_=rowSnpList;
		columnSnpList_=columnSnpList;				
		rowWeights_ = new double[rowSnpList_.size()];
		columnWeights_ = new double[columnSnpList_.size()];							
		Arrays.fill(rowWeights_, 1);
		Arrays.fill(columnWeights_, 1);	
		rowGenotypes_ = LinkageDisequilibrium.loadGenotypeIntoDenseMatrix(rowSnpList_);
		columnGenotypes_ = LinkageDisequilibrium.loadGenotypeIntoDenseMatrix(columnSnpList_);		
		
	}
	*/
	public LdMat (SnpWeightPairs rowSnpWeightPairs, ArrayList<Snp> columnSnpList){		
		rowSnpWeightPairs_=rowSnpWeightPairs;
		columnSnpList_=columnSnpList;
	}
		//LinkageDisequilibrium.computeCrossCorrelationMatrixMTJ(rowSnpList_,columnSnpList_);	
	public void processData(){
		rowSnpList_ =  rowSnpWeightPairs_.getSnps();
		columnWeights_ = new double[columnSnpList_.size()];
		rowWeights_ = rowSnpWeightPairs_.getWeightsNotList();
		Arrays.fill(columnWeights_, 1);
		rowGenotypes_ = LinkageDisequilibrium.loadGenotypeIntoDenseMatrix(rowSnpList_);
		columnGenotypes_ = LinkageDisequilibrium.loadGenotypeIntoDenseMatrix(columnSnpList_);
	
	}
	
	
	
	public boolean isCorrelationMat(){
		return rowSnpList_.equals(columnSnpList_);
	}
	public boolean isUnWeighted(){
		boolean myBool = true;
		for (double weight : rowWeights_)
			if (weight!=1) myBool=false;
		for (double weight : columnWeights_)
			if (weight!=1) myBool=false;		
		return myBool;
	}
	public DenseMatrix getWeightedCorrelationMatrix(){		
		DenseMatrix weightMatColWise = MTJConvenienceMethods.repMTJ(columnWeights_,columnWeights_.length, false);
		DenseMatrix weightMatRowWise = MTJConvenienceMethods.repMTJ(columnWeights_,columnWeights_.length, true);
		DenseMatrix Corr = MTJConvenienceMethods.calculateCovarianceMat(columnGenotypes_);			
		Corr = MTJConvenienceMethods.hadamardProduct(weightMatColWise, Corr);
		Corr = MTJConvenienceMethods.hadamardProduct(weightMatRowWise, Corr);   
		return Corr;
	}
		
	public DenseMatrix getWeighedCrossCorrelationMatrix(){		
		DenseMatrix weightMatRowWise = MTJConvenienceMethods.repMTJ(rowWeights_,columnWeights_.length, false);
		DenseMatrix crossCorr = MTJConvenienceMethods.calculateCrossCovarianceMat(rowGenotypes_, columnGenotypes_);		
		crossCorr = MTJConvenienceMethods.hadamardProduct(weightMatRowWise, crossCorr);		
		return crossCorr;
	}
	public DenseMatrix getCrossCorrelationMatrix(){		
		DenseMatrix crossCorr = MTJConvenienceMethods.calculateCrossCovarianceMat(rowGenotypes_, columnGenotypes_);				
		return crossCorr;
	}
	/** in case of normalized genotypes this is proportianal to the prior variance of the vegas statistic  
	 * because the trace of a matrix is equal to the sum of eigenvalues.
	 * */
	public double getSumOfWeightsSquared(boolean rowWeights){	
		double out = 0;
		double[] weights;
		if(rowWeights)
			weights = rowWeights_;
		else
			weights = columnWeights_;
		for (int i=0; i<weights.length; i++){
			out+=(weights[i]*weights[i]);
		}
			return out;
		
	}
	
	
	public void setRowSnpList(ArrayList<Snp> rowSnpList){rowSnpList_=rowSnpList;}
	public void setColSnpList(ArrayList<Snp> colSnpList){rowSnpList_=colSnpList;}	
	public void setColumnWeights(double[] columnWeights){columnWeights_ = columnWeights;}	
	public void setRowWeights(double[] rowWeights){rowWeights_ = rowWeights;}
	public double[] getColumnWeights(){return columnWeights_;}
	public double[] getRowWeights(){return rowWeights_;}

}