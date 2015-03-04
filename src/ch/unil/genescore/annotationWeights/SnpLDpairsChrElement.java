package ch.unil.genescore.annotationWeights;
import java.util.LinkedList;

import ch.unil.genescore.gene.GenomicElement;
public class SnpLDpairsChrElement extends GenomicElement {

	public String snpId_;
	public LinkedList<SnpLdPair>  ldValues_ = null;
	public SnpLDpairsChrElement(SnpLdPairs el){
		super(el.getSnp().getId());
		String snpChr = el.getSnp().chr_;
		int snpPosStart = el.getSnp().start_;
		int snpPosEnd = el.getSnp().end_;
		boolean snpStrand = el.getSnp().posStrand_;
		this.setPosition(snpChr, snpPosStart, snpPosEnd,snpStrand);
		ldValues_= el.getList();
	}
}
