package ch.unil.genescore.vegas;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;


public class ProjectionVegas {
	
	DenseMatrix crossCovariance_ = null;
	UpperSymmDenseMatrix covariance_ = null;
	DenseMatrix regularizedInvCovariance_ =null;
	DenseVector transformedZscores_ = null;
	DenseMatrix matToDecomposeMTJ_ = null;
	
	double fractionToBeExplained_;
	double conditionFraction_;
	//double conditionFractionCross_;
	double priorVar_;
	double posteriorVar_;
	double[] snpwiseVar_;
	private ArrayList<Double> snpScores_;
	private int numSnps_;
	/** crossLd has */
	public ProjectionVegas(ArrayList<Double> snpScores, DenseMatrix ld,DenseMatrix crossLd, double priorVar) {
		if (ld.numColumns()!=ld.numRows()){
			throw new RuntimeException("covariance mat not square");			
		}
		if (crossLd.numColumns()!=ld.numRows()){
			throw new RuntimeException("crossCovariance not compatible with covariance");			
		}
		
		fractionToBeExplained_=Settings.fractionToBeExplained_;
		conditionFraction_=Settings.conditionFraction_;
		priorVar_=priorVar;
		
		//conditionFractionCross_=Settings.conditionFractionCross_;
		setSnpScores(snpScores);
		setCovariance(ld);
		setCrossCovariance(crossLd);
	}
	public void processData(){
		//regularizeCrossCovariance();
		constructRegularizedInverseCovariance();
		constructMatToDecomposeMTJ();
		calculatePosteriorVar();
		calculateTransformedZscores();				
	}
	/**input: eigenvals: eigenvalues sorted in ascending order (MTJ-convention)
	 * get lowest index of eigenvalues that explain at least 'percentageExplained'
	 * important: eigenvalues have to be sorted in ascending order
	 * */
	private static int getEigenValsIndexAboveCutoff(double[] eigenvals, double percentageToBeExplained){
		double totalVariance = eigenvals.length;
		double varianceToBeExplained = totalVariance * percentageToBeExplained; 
		double explainedVariance = 0.0;
		for (int i=(eigenvals.length-1);i >= 0;i--){
			explainedVariance += eigenvals[i]; 
			if  (explainedVariance >= varianceToBeExplained)
				return i;
		}
		return -1;//should never be called;				
	}

	private void constructRegularizedInverseCovariance() {
		//TODO: remove writing files
//		WritingMethods.writeMTJ(covariance_, "run_Cov", Settings.writeCorFiles_);		
		UpperSymmDenseMatrix regularizedCovariance = MTJConvenienceMethods.regularizeMat(covariance_, conditionFraction_);
		SymmDenseEVD evd = new SymmDenseEVD(regularizedCovariance.numColumns(), true);		
		try {
			evd = SymmDenseEVD.factorize(regularizedCovariance);
		} catch (NotConvergedException e) {
			Main.println("Error: Eigendecomposition failed to converge");
		}
		DenseMatrix eigenVects = evd.getEigenvectors();	
		double[] eigenVals = evd.getEigenvalues();			
		// find how many eigenvectors should be used.
		int index = getEigenValsIndexAboveCutoff(eigenVals, fractionToBeExplained_);
		//int index = 0;
		DenseMatrix subEigenVects=MTJConvenienceMethods.getSubDenseMatrixMTJ(eigenVects,0,index);
		double[] subInvEigenVals = new double[eigenVals.length-index];
		for (int i=0 ; i < subInvEigenVals.length ; i++)
			subInvEigenVals[i]=1/eigenVals[index+i];
		DenseMatrix subInvEigenValsMat = MTJConvenienceMethods.diagMTJ(subInvEigenVals);		
		regularizedInvCovariance_ = MTJConvenienceMethods.doGammaLambdaGammaTMTJ(subEigenVects,subInvEigenValsMat); 		
		//TODO: remove writing files
//		WritingMethods.writeMTJ(regularizedInvCovariance_, "run_regInvCov", Settings.writeCorFiles_);
	}
private void constructMatToDecomposeMTJ() {
	//TODO: remove writing files
	//WritingMethods.writeMTJ(regularizedInvCovariance_, "run_regInvCov", Settings.writeCorFiles_);
		DenseMatrix interim= null;
		interim = MTJConvenienceMethods.doGammaLambdaGammaTMTJ(crossCovariance_,regularizedInvCovariance_);
//		WritingMethods.writeMTJ(interim, "run_matToDecompose", Settings.writeCorFiles_);
		matToDecomposeMTJ_ = interim;				
	}
private DenseVector calculateTransformedZscores() {
	double[] zscores = new double[snpScores_.size()]; 
	for (int i=0 ; i< snpScores_.size(); i++){
				zscores[i]=snpScores_.get(i);
	}
	DenseVector Zscores = new DenseVector(zscores);
	DenseMatrix interim= new DenseMatrix(crossCovariance_.numRows(),crossCovariance_.numColumns());
	
	crossCovariance_.mult(regularizedInvCovariance_,interim);			
	transformedZscores_ = new DenseVector(crossCovariance_.numRows());
	interim.mult(Zscores, transformedZscores_);
	return transformedZscores_;
}

	
private void calculatePosteriorVar(){
	double projectionVar = MTJConvenienceMethods.traceMTJ(matToDecomposeMTJ_);
	posteriorVar_ = priorVar_ - projectionVar;
}
public void setSnpScores(ArrayList<Double> snpScores) {
	
	snpScores_ = snpScores;
	numSnps_ = snpScores_.size();
}

private void setCrossCovariance(DenseMatrix crossCovariance){	
	//TODO: just for debuggin
	//crossCovariance_ = new DenseMatrix(regularizeMat(tmp, conditionFraction_));
	//crossCovariance_ = new DenseMatrix(MTJConvenienceMethods.regularizeMat(crossCovariance, 0.0));
	crossCovariance_ = new DenseMatrix(crossCovariance);
}	
//TODO: unneccessary?
//private void regularizeCrossCovariance(){		
//	crossCovariance_.scale(1-conditionFractionCross_);
//}	

//TODO: change setting in AnalyticVegas
public void setCovariance(DenseMatrix covariance){		
	covariance_ = new UpperSymmDenseMatrix(covariance);
	//covariance_ = new UpperSymmDenseMatrix(regularizeMat(tmp, conditionFraction_));	
	//System.out.println("blub");
}	

public void setFractionToBeExplained(double fractionToBeExplained) {fractionToBeExplained_=fractionToBeExplained;}
public void setConditionFraction(double conditionFraction) {conditionFraction_=conditionFraction;}
public DenseMatrix getRegularizedInvCovariance(){return regularizedInvCovariance_;}
public DenseVector getTransformedZscores(){return transformedZscores_;}
public DenseMatrix getMatToDecomposeMTJ(){return matToDecomposeMTJ_;}

}