package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.vegas.DataInconsistencyException;
import ch.unil.genescore.vegas.Snp;

public class EqtlGenePair {
	
	String geneName_ = null; 
	HashMap<String, Snp> snpList_ = null;
	ArrayList<Snp> overlappedSnpList_ = null;	
	double[] overlappedEqtlZscores_ = null;
	
	
	public EqtlGenePair(){				
	}
	
	public EqtlGenePair(String geneName, HashMap<String, Snp> snpList){
		
		geneName_=geneName;
		snpList_=snpList;
	}
	
	public EqtlGenePair(String geneName, ArrayList<Snp> snpAr){
		
		geneName_=geneName;
		snpList_  = new HashMap<String, Snp>();
		
		for (int i=0;i<snpAr.size();i++){
			Snp curSnp = snpAr.get(i);
			snpList_.put(curSnp.id_, curSnp);
		}		
	}
	
	/** gets zscores of snps that overlap the two Snp sets (in snpList_ and geneSnps).*/
	public void overlapSnps(List<Snp> geneSnps){
		ArrayList<Double> zScores = new ArrayList<Double>();
		overlappedSnpList_ = new ArrayList<Snp>();
		int len=geneSnps.size();
		for (int i=0; i < len;i++){
			Snp gwasSnp = geneSnps.get(i);
			String currentSnpId = gwasSnp.id_;
			if(snpList_.containsKey(currentSnpId)){
				Snp eqtlSnp = snpList_.get(currentSnpId);
				try{
					compareAlleleDirection(gwasSnp,eqtlSnp);
					zScores.add(eqtlSnp.getZscore());
					overlappedSnpList_.add(gwasSnp);					
				}
				catch(DataInconsistencyException e){									
			        System.err.println("Caught DataInconsistencyException: " + e.getMessage());			       			     
					continue;						
				}				
			}		
		}
		
		overlappedEqtlZscores_ = ConvenienceMethods.arListToDoubleAr(zScores);		
	}
	
	
	/** gets zscores of snps that overlap the two Snp sets (in snpList_ and geneSnps).*/
	public void overlapSnpWithoutCheck(List<Snp> List){
		
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
	
	/** create overlap list that only consists of one snp where 
	 *  which is the top p-value in the eqtl set.*/
	public void getSnpWithTopEqtlValue(List<Snp> geneSnps){
		
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
	}
	
	
	private void compareAlleleDirection(Snp gwasSnp,Snp eqtlSnp) throws DataInconsistencyException{
		
		if(gwasSnp.getMajorAllele()==eqtlSnp.getMajorAllele() && gwasSnp.getMinorAllele()==eqtlSnp.getMinorAllele())
			return;
		else if(gwasSnp.getMajorAllele()==eqtlSnp.getMinorAllele() && gwasSnp.getMinorAllele()==eqtlSnp.getMajorAllele())
			eqtlSnp.setZscore((-1)*eqtlSnp.getZscore());
		else{
			throw new DataInconsistencyException("Major and Minor allele conflict between eqtl data and reference population");
			}
	}
	public double[] getOverlappedEqtlZscores(){
		return(overlappedEqtlZscores_);		
	}
	public ArrayList<Snp> getOverlappedSnpList(){
		return(overlappedSnpList_);
	}
	public HashMap<String, Snp> getSnpList(){
		return(snpList_);
	}
	
	public void setGeneName(String str){
		geneName_ = str;
	}		
	public  void setSnpList(HashMap<String, Snp> l){
		snpList_ = l;
	}
	
	
	protected void setEqtlZscoresFromSnpList(){
		ArrayList<Double> zScores = new ArrayList<Double>();
		for (Snp eqtlSnp : overlappedSnpList_){					
			zScores.add(eqtlSnp.getZscore());			
		}
		overlappedEqtlZscores_ = ConvenienceMethods.arListToDoubleAr(zScores);	
	}
	
	public  void setSnpList(ArrayList<Snp> snpAr){
		
		snpList_  = new HashMap<String, Snp>();	
		for (int i=0;i<snpAr.size();i++){
			Snp curSnp = snpAr.get(i);
			snpList_.put(curSnp.id_, curSnp);
		}
	}
}
