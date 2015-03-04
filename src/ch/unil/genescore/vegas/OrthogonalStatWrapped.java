
package ch.unil.genescore.vegas;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import ch.unil.genescore.main.Utils;



public class OrthogonalStatWrapped extends GeneScoreEvaluator {

	OrthogonalStat OrthogonalStat_ = null;
	int[] allNrComponents_ = null;
	double[] allGeneScores_ = null;
	public OrthogonalStatWrapped(ArrayList<Double> snpScores, DenseMatrix ld, int[] allNrComponents) {
		allNrComponents_ = allNrComponents;
		allGeneScores_ = new double[allNrComponents.length];
		int firstNrComponents = allNrComponents[0];
		OrthogonalStat_ = new OrthogonalStat(snpScores, ld, firstNrComponents);	
	} 
	@Override
	public String getResultsAsString() {
		// TODO Auto-generated method stub
		String line = "";
		for (int i=0; i < allGeneScores_.length; ++i){
			double geneScore_ = allGeneScores_[i];		
			line += "\t" + Utils.toStringScientific10(geneScore_);
		}
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

	@Override
	public boolean computeScore() {
		for (int i = 0; i < allNrComponents_.length; ++i){
			OrthogonalStat_.setNrComponents(allNrComponents_[i]);
			OrthogonalStat_.computeScore();
			allGeneScores_[i] = OrthogonalStat_.getGeneScore(); 
		}
		
		// TODO return true if success, false if fail
		return true;
		
	}

	
	/** Get gene score(s) */
	public double[] getScore() { return allGeneScores_; }
	@Override
	public void setDataFromGeneData(GeneData Dat) {
		// TODO Auto-generated method stub
		
	}

}
