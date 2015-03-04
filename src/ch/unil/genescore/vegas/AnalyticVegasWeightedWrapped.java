package ch.unil.genescore.vegas;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import ch.unil.genescore.main.Utils;

public class AnalyticVegasWeightedWrapped extends GeneScoreEvaluator{

	double[] allGeneScores_ = null;
	ArrayList<Double> deltaVect_ = null;
	int deltaPos_ = -1;
	AnalyticVegas analyticVegasWeighted_ = null;


	// ============================================================================
	// PUBLIC METHODS

	/** Constructor */
	public AnalyticVegasWeightedWrapped(ArrayList<Double> snpScores, DenseMatrix ld, ArrayList<Double> deltaVect) {
		double[] emptyWeights = new double[snpScores.size()];
		analyticVegasWeighted_ = new AnalyticVegas(snpScores, ld, emptyWeights);		
		deltaVect_ = deltaVect;
	}	
	
	
	// ----------------------------------------------------------------------------

	/** Compute gene score for each delta */
	public boolean computeScore() {

		boolean allSuccess = true;
		allGeneScores_ = new double[deltaVect_.size()];
		
		for (int i=0; i<deltaVect_.size(); i++){
			moveDelta();
			//long t0 = System.currentTimeMillis();
			boolean success =  analyticVegasWeighted_.computeScore();
			allSuccess &= success;
			
			//long t1 = System.currentTimeMillis();
			//		Main.print("\t" + Utils.chronometer(t1 - t0));			
			double aa=analyticVegasWeighted_.getScore()[0];
			//		System.out.println(aa[0]);
			allGeneScores_[i] = aa; 
			//allGeneScores_.add(AnalyticVegasWeighted_.getScore()); 
		}
		return allSuccess;
	}

	
	// ----------------------------------------------------------------------------

	/** Get the score(s) */
	public double[] getScore() { return allGeneScores_;	}
	
	
	// ----------------------------------------------------------------------------

	/** Get string representation of result that will be written to the output file */
	public String getResultsAsString() {

		// Tab-separated list of gene scores
		String line = "";
		for (double score : allGeneScores_)
			line += "\t" + Utils.toStringScientific10(score);
		return line;
	}
	
	// TODO
	public String getResultsAsStringHeader() {
		
		return "TBD";
	}

	
	/** Get output to be printed on console after computing score for a gene */
	// TODO
	public String getConsoleOutput() {
		
		return "";
	}

	
	// ============================================================================
	// PRIVATE METHODS

	/** Initialize AnalyticVegasWeighted with the next delta */
	private void moveDelta() {
		
		deltaPos_++;		
		double[] weightVect=analyticVegasWeighted_.getSnpWeights();
		double[] newWeightVect = new double[weightVect.length];
		for (int i=0 ; i<weightVect.length; i++){
			newWeightVect[i]=deltaVect_.get(deltaPos_);
		}
		analyticVegasWeighted_.initialize(newWeightVect);
	}


	@Override
	public void setDataFromGeneData(GeneData Dat) {
		//setCovariance(Dat.getCorr());
    	//setWeights(Dat.getWeights());
		
		
	}

}


