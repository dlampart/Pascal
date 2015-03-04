package  ch.unil.genescore.annotationWeights;
import java.util.LinkedList;

import ch.unil.genescore.vegas.Snp;
//import  ch.unil.genescore.annotationWeights.SnpLdPair;
public class SnpLdPairs {

	Snp Snp_;
	LinkedList<SnpLdPair>  ldValues_ = new LinkedList<SnpLdPair>();
	public SnpLdPairs(Snp mySnp) {
		
		Snp_ = mySnp;
		
		// TODO Auto-generated constructor stub
	}
	void appendToLDList(String snpId, double ldValue){
			ldValues_.add(new SnpLdPair(snpId, ldValue));
	}
	public Snp getSnp(){
		
		return Snp_;
	}
	public LinkedList<SnpLdPair>  getList(){
		
		return ldValues_;
	}
	public double getSnpPos(){
		return Snp_.getStart();
	}
	public String getSnpId(){
		return Snp_.getId();
	}

}
