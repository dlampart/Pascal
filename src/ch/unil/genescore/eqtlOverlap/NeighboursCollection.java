package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;

public class NeighboursCollection {

	HashMap<String, ArrayList<String>> geneNeighbours_ = new HashMap<String, ArrayList<String>>();
	
	public HashSet<String> getAllNeighbourName(ArrayList<String> genes){
		HashSet<String> allNames = new HashSet<String>();
		for (String geneSymbol : genes){			
			allNames.addAll(geneNeighbours_.get(geneSymbol));
		}
		allNames.removeAll(genes);
		return allNames;	
		}
	
	
	public void findNeighbours(Collection<Gene> genes){
		TreeSet<Gene> geneTree = new TreeSet<Gene>(genes);
		Iterator<Gene> lowerIt = geneTree.iterator();
		Iterator<Gene> upperIt = geneTree.iterator();
		int cutoff=1000000;
		Gene lowerGene = lowerIt.next();
		Gene upperGene = upperIt.next();		
		for (Gene gene : geneTree){
			ArrayList<String> geneNeighbourNames = new ArrayList<String>();
			if(lowerIt.hasNext()){					
				if (!lowerGene.extensionPartiallyOverlapsElement(gene,0, cutoff)){
					lowerGene=lowerIt.next();
				}
			}
			if(upperIt.hasNext()){					
				if (upperGene.extensionPartiallyOverlapsElement(gene,-cutoff, 0)){
					upperGene=upperIt.next();
				}
			}								
			NavigableSet<Gene> naviSet = geneTree.subSet(lowerGene, true, upperGene, false);
			for(Gene gene2 :naviSet){				
				geneNeighbourNames.add(gene2.symbol_);
			}	
			geneNeighbours_.put(gene.symbol_, geneNeighbourNames);
		}	
	}
	public void clearMap(){
		geneNeighbours_ = new HashMap<String, ArrayList<String>>();
	}
}
