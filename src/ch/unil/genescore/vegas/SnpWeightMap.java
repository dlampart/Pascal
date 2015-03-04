package ch.unil.genescore.vegas;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Settings;
/** */
public class SnpWeightMap {
	HashMap<String, Double> genewiseWeightMap_  = null; 
	private ArrayList<Snp> overlappedSnps_ = new ArrayList<Snp>();
	private ArrayList<Double> overlappedWeights_ = new ArrayList<Double>();	
	private SnpWeightPairs overlappedSnpsWithWeights_ = null;
	
	public SnpWeightMap(){
		genewiseWeightMap_  = new HashMap<String, Double>();
	}
	/** set by hand ! set not only for snp class!*/ 
	public SnpWeightMap(List<? extends GenomicElement> myElements, double weight){
		genewiseWeightMap_  = new HashMap<String, Double>();
		for (GenomicElement currentElement : myElements){
			genewiseWeightMap_.put(currentElement.id_, weight);
		}
	}
	public void findOverlappedSnps(List<Snp> snpList){
		ArrayList<Snp> snpsInMap = new ArrayList<Snp>();		 
		ArrayList<Double> WeightsOfSnpsInMap = new ArrayList<Double>();		 
		//for (String inMap : genewiseWeightMap_.keySet())
			//System.out.println(inMap);
		//System.out.println("number of relevant snps");
		System.out.println(genewiseWeightMap_.keySet().size());
		for (Snp currentSnp : snpList){	
			if (currentSnp.getMaf() < 0){
				throw new RuntimeException("Maf hasn't been set yet!");				
			}
			if (currentSnp.getMaf() < Settings.useMafCutoffForProjection_)
				continue;			
				
			//System.out.println(currentSnp.id_);
			String currentId = currentSnp.getId();
			if (genewiseWeightMap_.containsKey(currentId)){
		//		System.out.println(currentSnp.id_);
				snpsInMap.add(currentSnp);				
				WeightsOfSnpsInMap.add(genewiseWeightMap_.get(currentId));				
			}								
		}
		//System.out.println("number of relevant snps projected upon");
		//System.out.println(snpsInMap.size());	
		overlappedSnpsWithWeights_ = new SnpWeightPairs(snpsInMap,WeightsOfSnpsInMap);
		overlappedSnps_ = snpsInMap;		
		overlappedWeights_ = WeightsOfSnpsInMap;
	}
		
	/**reads in SnpWeightPairList from a split line. format: each String has format snpId;weight .*/
	public void readMapFromSplitLine(String[] inputStrings){
		String snpId = null;
		double value;		
		for (String inputString : inputStrings){
			String[] bothValues = inputString.split(";");
			snpId = bothValues[0];
			value = Double.parseDouble(bothValues[1]);
			if (!genewiseWeightMap_.containsKey(snpId)){
				genewiseWeightMap_.put(snpId, value);	
			}
			else {
				double priorWeight = genewiseWeightMap_.get(snpId);
				double newWeight = priorWeight * value;
				genewiseWeightMap_.put(snpId, newWeight);	
			}
		}			
	}
	public void updateWeights(SnpWeightMap other){
		Double newVal = null;		
		Iterator<String> iter = genewiseWeightMap_.keySet().iterator();		
		while(iter.hasNext()){
			String key = iter.next();
			if (other.checkExistence(key)){		
				double priorVal = genewiseWeightMap_.get(key);
				double otherVal = other.getWeight(key);
				newVal = new Double(otherVal * priorVal);
				genewiseWeightMap_.put(key, newVal);
			}
		}		
	}
	
	public boolean checkExistence(String id){
		return genewiseWeightMap_.containsKey(id);
	}
	public double getWeight(String id){
		return genewiseWeightMap_.get(id);
	}
	public ArrayList<Snp> getOverlappedSnps(){return overlappedSnps_;}
	public ArrayList<Double> getOverlappedWeights(){return overlappedWeights_;}
	public SnpWeightPairs getOverlappedSnpsWithWeights(){return overlappedSnpsWithWeights_;}
	public double[] getOverlappedWeightsNotList(){
		int n=overlappedWeights_.size();
		double[] outAr  = new double[n];
		for (int i=0 ; i < n ; i++)
			outAr[i]=overlappedWeights_.get(i);		
		return outAr;
	}
	public ArrayList<String> getSnpIds(){
		Iterator<String> iter = genewiseWeightMap_.keySet().iterator();		
		ArrayList<String> myAr = new ArrayList<String>(genewiseWeightMap_.size()); 				
		while (iter.hasNext()){
			myAr.add(iter.next());		
		}
		return myAr;
	} 
}
	
