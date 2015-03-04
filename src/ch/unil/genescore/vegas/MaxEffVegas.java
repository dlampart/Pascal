package ch.unil.genescore.vegas;

import java.util.ArrayList;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.MaxVegas.Status;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

public class MaxEffVegas extends MaxVegas {
	
	public MaxEffVegas(){
		super();
		pruningCutoff_=Settings.maxPruningCutoff_;	
	}
	public MaxEffVegas(ArrayList<Double> snpScores, UpperSymmDenseMatrix ld, double[] weights) {			
		super(snpScores, ld, weights);
		}

	public MaxEffVegas(ArrayList<Double> snpScores, DenseMatrix ld,	double pruningCutoff) {
		super(snpScores, ld, pruningCutoff);
		
	}
	
	
	@Override
	protected void calculateGenescore(){	
	geneScore_=computeEffectivePval();	
	status_=Status.EFF_NR_OF_TEST_ADJ;
	
	}	
}
