package ch.unil.genescore.annotationWeights;


public class SnpLdPair {
	
	public String snpId_;
	public double ldValue_;
	public SnpLdPair(String snpId,double ldValue){
		
		snpId_ = snpId;
		ldValue_ = ldValue;
	}
	public String getId(){
		return snpId_;
	}
}
