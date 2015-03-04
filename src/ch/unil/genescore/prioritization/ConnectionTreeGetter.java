package ch.unil.genescore.prioritization;

import java.util.Collection;
import com.google.common.collect.TreeMultimap;
import ch.unil.genescore.gene.Gene;

/**
 * input Gene currentGene:
 * output TreeMultiMap<Double,Gene> sortedByValue sortedByValue:
 * */
public interface ConnectionTreeGetter {
	
	/**loads the network data*/
	public void loadNetworkData();
	
	/**
	 * input Gene currentGene:
	 * 
	 * output TreeMap<Double,Gene> sortedByValue
	 * returns a Tree of genes that are connected to currentGene in the network sorted by connection strength.
	 * 
	 * bordercases
	 * returns empty tree if gene is not found in network; 
	 * 
	 * first entry: <1,currentGene>
	 * */
	public TreeMultimap<Double,Gene> returnReverseOrderedConnectionTree(Gene currentGene);
	
	/**sets Genes that are potentially returned*/
	public void setGenes(Collection<Gene> Genes);	
}
