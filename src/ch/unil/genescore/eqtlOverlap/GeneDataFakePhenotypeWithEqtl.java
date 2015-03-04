package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;

import ch.unil.genescore.vegas.GeneDataEqtlZscoreHolder;
import ch.unil.genescore.vegas.GeneDataFakePhenotype;
import ch.unil.genescore.vegas.Snp;

public class GeneDataFakePhenotypeWithEqtl extends GeneDataFakePhenotype implements GeneDataEqtlZscoreHolder {

	double[] eqtlZscores_ = null;
	public GeneDataFakePhenotypeWithEqtl(ArrayList<Snp> snpList,double[] FakeSignal_, double[] eqtlZscores) {
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
