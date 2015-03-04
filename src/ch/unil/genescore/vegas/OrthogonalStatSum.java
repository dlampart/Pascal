package ch.unil.genescore.vegas;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import ch.unil.genescore.main.Main;

public class OrthogonalStatSum extends OrthogonalStat {
	
	public OrthogonalStatSum(){
		super();
	}
	public OrthogonalStatSum(ArrayList<Double> snpScores, DenseMatrix ld, int nrComponents) {
		super(snpScores,ld, nrComponents);
	}
	public OrthogonalStatSum(ArrayList<Double> snpScores, DenseMatrix ld) {
		super(snpScores,ld, 1);
	}
	
	
	public boolean computeScore(){
		double sumStatistic = 0;
		for (int i=0; i< indepScaledStats_.length ;i++){
				sumStatistic += indepScaledStats_[i]*indepScaledStats_[i];
		}		
		if (indepScaledStats_.length < 1){						
			geneScore_ = 1.0;
			return true;
		}
		ChiSquaredDistribution chiSquareDistN = new ChiSquaredDistribution(indepScaledStats_.length);
		geneScore_= 1.0 - chiSquareDistN.cumulativeProbability(sumStatistic);			
		return true;
	}
}
