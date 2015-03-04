package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.DistributionMethods;
import ch.unil.genescore.vegas.GeneData;
import ch.unil.genescore.vegas.GeneDataEqtlZscoreHolder;
import ch.unil.genescore.vegas.GeneScoreEvaluator;
import ch.unil.genescore.vegas.MTJConvenienceMethods;

public class EqtlEvaluator extends GeneScoreEvaluator {

	//GeneDataWithEqtl  geneData_ = null;
	protected DenseMatrix LdMat_ = null;
	protected double[] eqtlZscores_ = null;
	protected double[] gwasZscores_ = null;
	
	double geneScore_ ;
	public EqtlEvaluator(){}
	
	public EqtlEvaluator(GeneDataEqtlZscoreHolder geneData){
		LdMat_ = geneData.getCorr();
		eqtlZscores_ = geneData.getEqtlZscores();		
		gwasZscores_ = ConvenienceMethods.arListToDoubleAr(geneData.getScores());
	}
	

	
	@Override
	public boolean computeScore() {
		
		DenseMatrix mat = LdMat_;
		mat = MTJConvenienceMethods.scaleToDiagOne(MTJConvenienceMethods.regularizeMat(mat, 0.05));
		DenseMatrix beta = new DenseMatrix(new DenseVector(eqtlZscores_));
		double[] xD = gwasZscores_;
		DenseMatrix x = new DenseMatrix(new DenseVector(xD));
		DenseMatrix invSigmaBeta = new DenseMatrix(beta.numRows(),beta.numColumns());
		mat.solve(beta, invSigmaBeta);
		DenseMatrix betaInvSigmaBeta = new DenseMatrix(1,1);
		DenseMatrix xInvSigmaBeta = new DenseMatrix(1,1);
		invSigmaBeta.transAmult(beta, betaInvSigmaBeta);
		invSigmaBeta.transAmult(x, xInvSigmaBeta);
		double zScore = xInvSigmaBeta.get(0,0)/Math.sqrt(betaInvSigmaBeta.get(0,0));
		// TODO Auto-generated method stub
		geneScore_ = 2*DistributionMethods.normalCumulativeProbability((-1)*Math.abs(zScore));
		return true;
	}

	@Override
	public double[] getScore() {
		
			double[] result = { geneScore_ };
			return result;
	}
	public void setScore(double myScore) {
		
		  geneScore_ = myScore;		
}

	@Override
	/** Get string representation of result that will be written to the output file */
	  public String getResultsAsString() {
		  
		  String line = "\t" + Utils.toStringScientific10(geneScore_);		 
	        return line;
	    }	

	    /** Get a header line corresponding to getResultsAsString() */
	    public String getResultsAsStringHeader() {

	        //String header = "\tpvalue\tnumSnps";
	    	String header = "\tpvalue";
	  //      if (Settings.writeDetailedOutput_)
	   //         header += "\tavgSnpCorrelation\tstatus";
	        return header;
	    }

	@Override
	public void setDataFromGeneData(GeneData geneData) {			
		//LdMat_ = geneData.getCorr();
		//eqtlZscores_ = geneData.getEqtlZscores();		
		//gwasZscores_ = ConvenienceMethods.arListToDoubleAr(geneData.getScores());
	}
	    
	public void setDataFromGeneData(GeneDataEqtlZscoreHolder geneData) {			
		LdMat_ = geneData.getCorr();
		eqtlZscores_ = geneData.getEqtlZscores();		
		gwasZscores_ = ConvenienceMethods.arListToDoubleAr(geneData.getScores());
	}
}
