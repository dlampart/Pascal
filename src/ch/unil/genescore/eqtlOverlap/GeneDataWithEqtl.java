package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;

import ch.unil.genescore.vegas.GeneDataEqtlZscoreHolder;
import ch.unil.genescore.vegas.GeneData;
import ch.unil.genescore.vegas.Snp;

public class GeneDataWithEqtl extends GeneData implements GeneDataEqtlZscoreHolder {
	
	double[] eqtlZscores_ = null;
	public GeneDataWithEqtl(ArrayList<Snp> snpList, double[] eqtlZscores) {
		super(snpList);
		eqtlZscores_= eqtlZscores;
		
	}
	@Override
	public void setEqtlZscores(double[] zScores) {		
		eqtlZscores_= zScores;
	}
	@Override
	public double[] getEqtlZscores() {
		// TODO Auto-generated method stub
		return eqtlZscores_;
	}
	
	
}
