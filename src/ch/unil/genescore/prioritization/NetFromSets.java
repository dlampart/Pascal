package ch.unil.genescore.prioritization;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.pathway.GeneSet;
import ch.unil.genescore.pathway.GeneSetLibrary;

import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;

public class NetFromSets implements ConnectionTreeGetter {

	//GeneSetLibrary geneSets_ = null; 
	ArrayList<GeneSet> geneSets_ = null;
	protected FileParser parser_ = null;
	Collection<Gene> geneLists = null;
	
	@Override
	public void loadNetworkData() {
		
		
		//geneSets_.getGenes()
		
		
		
		
		
	}
	
	
	/**
	 * input Gene queryGene:
	 * output TreeMap<Double,Gene> sortedByValue
	 * returns a Tree of genes that are connected to currentGene that share a pathway sorted by connection strength.
	 * connection strength is defined as 1/pathway_size
	
	 * returns empty tree if gene is not found in network; 
	 * 
	 * first entry: <1.1,queryGene>
	 * */
	@Override
	public TreeMultimap<Double, Gene> returnReverseOrderedConnectionTree(
			Gene queryGene) {
		TreeMultimap<Double,Gene> sortedByValue = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
		for(GeneSet set : geneSets_){
			HashSet<Gene> gset = set.getGenes();
			if(gset.contains(queryGene)){
				Double weight = 1/ (double) (gset.size());
				for (Gene g : gset){
					if (g != queryGene)
						sortedByValue.put(weight,g);
				}
			}
			sortedByValue.put(1.1, queryGene);	
		}
		return sortedByValue;
	}

	@Override
	public void setGenes(Collection<Gene> Genes) {
		HashMap<String, Gene> geneHash  = new HashMap<String, Gene>();
		for(Gene g: Genes){
			geneHash.put(g.id_, g);
		}
		geneSets_ = new GeneSetLibrary(Settings.netPath_, geneHash).getGeneSets();
	}

}
