package ch.unil.genescore.vegas;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MetaGenewiseSnpWeights extends GenewiseSnpWeights {

	private String separator_ = null;
	public MetaGenewiseSnpWeights(){
		super();
	}
	public MetaGenewiseSnpWeights(HashMap<String, SnpWeightMap> snpWeightsForEachGene){
		super(snpWeightsForEachGene);
	}
	
	public String[] getMetaIds(String metaGeneId){
		String[] splitVals = metaGeneId.split(separator_, -1);
		String[] stringArrayOut = null;
		int startPos;
		if (splitVals[0].equals("meta"))
			startPos=1;
		else
			startPos=0;
		int count=0;
		stringArrayOut= new String[splitVals.length-startPos];
		for ( int i=startPos; i<splitVals.length; i++){						
			stringArrayOut[count]=splitVals[i];
			count++;
		}
		return stringArrayOut;
	}
	
	public SnpWeightPairs getOverlappedSnpsWithWeights(String metaGeneId, List<Snp> SnpList){
		
		String[] separatedGeneIds = getMetaIds(metaGeneId);
		SnpWeightPairs fusedSnpWeightPairs = new SnpWeightPairs();
		for (String geneId : separatedGeneIds){
			 fusedSnpWeightPairs.fuseSnpWeightPairs(super.getOverlappedSnpsWithWeights(geneId, SnpList));					
		}
		return fusedSnpWeightPairs;		
	}
}
