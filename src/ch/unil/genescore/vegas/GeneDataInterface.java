package ch.unil.genescore.vegas;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;

public interface GeneDataInterface {

	public DenseMatrix getCorr();
	public ArrayList<Double> getScores();
	public double[] getWeights();
	public void processData();	
}
