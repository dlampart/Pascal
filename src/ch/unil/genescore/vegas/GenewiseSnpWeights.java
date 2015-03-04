package ch.unil.genescore.vegas;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Settings;

public class GenewiseSnpWeights implements AccessGenewiseSnpWeights {

	
	protected  HashMap<String, SnpWeightMap> snpWeightsForEachGene_ = null;
	
	public GenewiseSnpWeights(){
		snpWeightsForEachGene_ = new HashMap<String, SnpWeightMap>();
	}
	/** set GenewiseSnpWeights Directly */
	public GenewiseSnpWeights(HashMap<String, SnpWeightMap> snpWeightsForEachGene){
		snpWeightsForEachGene_=snpWeightsForEachGene;
	}
    
	/** load snp weights from File;path to file: in Settings.geneWiseSnpWeightsFile_;*/
	public void loadSnpWeightsForEachGeneFromFile(Collection<Gene> geneList){
    	snpWeightsForEachGene_ = new HashMap<String, SnpWeightMap>();
    	String filename = Settings.geneWiseSnpWeightsFile_;
    	FileParser parser = new FileParser(filename,"::");
    	while (true){
            String[] nextLine = parser.readLine();
            if (nextLine == null)
                break;
            if (nextLine.length > 2){
            	 throw new RuntimeException("snpWeightsFile has format problem, separator :: detected more than once on a line.");
            }
            String geneId = nextLine[0];
            String[] SnpIdsAndWeights = nextLine[1].split("\t");
            SnpWeightMap currentMap = new SnpWeightMap();
            currentMap.readMapFromSplitLine(SnpIdsAndWeights);
            //currentMap.findOverlappedSnps(snpSet);          
            snpWeightsForEachGene_.put(geneId, currentMap);
    	}
	}
	/** for a list of snps, returns the ones that have a weighting for the given gene */
	@Override
	public SnpWeightPairs getOverlappedSnpsWithWeights(String geneId, List<Snp> SnpList){		
		
		SnpWeightPairs relevantSnpsAndWeights;
		if(!snpWeightsForEachGene_.containsKey(geneId)){
    		return null;
    	}  	
		else {
			SnpWeightMap currentMap = snpWeightsForEachGene_.get(geneId);   
	    	currentMap.findOverlappedSnps(SnpList);
	    	relevantSnpsAndWeights = currentMap.getOverlappedSnpsWithWeights();
		}
		return relevantSnpsAndWeights;				
	}
	/**returns all snpIds that are in the snpWeightsForEachGene_ (basically this is all snps that are loaded
	 * from the weight file)*/
	public HashSet<String> returnAllSnpIds(){
		HashSet<String> allSnpIds = new HashSet<String>();
		for (SnpWeightMap currentSnpWeightMap : snpWeightsForEachGene_.values()){
			for (String relevantSnpId : currentSnpWeightMap.genewiseWeightMap_.keySet()){
				allSnpIds.add(relevantSnpId);
			}			
		}			
		return allSnpIds;
	}
	public void updateWeights(GenewiseSnpWeights other){
	
		Iterator<String> iter = snpWeightsForEachGene_.keySet().iterator();		
		while(iter.hasNext()){
			String key = iter.next();
			SnpWeightMap val = snpWeightsForEachGene_.get(key);
			if (other.checkExistence(key)){				
				SnpWeightMap otherVal = other.getSnpWeightMap(key); 
				val.updateWeights(otherVal);
						
			}
		}		
		
	}
	public boolean checkExistence(String id){
		return snpWeightsForEachGene_.containsKey(id);
	}
	
	
	public HashMap<String, SnpWeightMap> getHash(){
		return snpWeightsForEachGene_;
	}	
	public void combineWithOtherGenewiseSnpWeights (GenewiseSnpWeights other){
		HashMap<String, SnpWeightMap> otherHash = other.getHash();
		snpWeightsForEachGene_.putAll(otherHash);
	}
	public SnpWeightMap getSnpWeightMap(String id){
		return snpWeightsForEachGene_.get(id);
	}
}
