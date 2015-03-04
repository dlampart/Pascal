package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.List;

import ch.unil.genescore.vegas.Snp;

public class TopSnpEqtlGenePair extends EqtlGenePair {
	
	/** create overlap list that only consists of one snp where 
	 *  which is the top p-value in the eqtl set.*/
	@Override
	public void overlapSnps(List<Snp> geneSnps){
		
		overlappedSnpList_ = new ArrayList<Snp>();
		Snp myTopSnp= null;	
		Snp gwasSnp = null;
		double myTopSnpPval=1.0;
		int len=geneSnps.size();
		for (int i=0; i < len;i++){
			gwasSnp = geneSnps.get(i);
			String currentSnpId = gwasSnp.id_;
			if(snpList_.containsKey(currentSnpId)){
				Snp eqtlSnp = snpList_.get(currentSnpId);
				if(myTopSnpPval>eqtlSnp.getPval()){
					myTopSnp=gwasSnp;
					myTopSnpPval=eqtlSnp.getPval();
				}
					
			}				
		}
		if(myTopSnp!=null)
			overlappedSnpList_.add(myTopSnp);	
	setEqtlZscoresFromSnpList();
	}
}
