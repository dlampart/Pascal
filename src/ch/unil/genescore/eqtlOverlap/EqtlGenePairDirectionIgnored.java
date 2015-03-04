package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.List;

import ch.unil.genescore.vegas.Snp;

public class EqtlGenePairDirectionIgnored extends EqtlGenePair {
	
	/** gets zscores of snps that overlap the two Snp sets (in snpList_ and geneSnps).*/
	public void overlapSnp(List<Snp> List){
		
		overlappedSnpList_ = new ArrayList<Snp>();
		int len=List.size();
		for (int i=0; i < len;i++){
			Snp gwasSnp = List.get(i);
			String currentSnpId = gwasSnp.id_;
			if(snpList_.containsKey(currentSnpId)){
				Snp eqtlSnp = snpList_.get(currentSnpId);
					overlappedSnpList_.add(gwasSnp);											
			}			
		}
		setEqtlZscoresFromSnpList();
	}

}
