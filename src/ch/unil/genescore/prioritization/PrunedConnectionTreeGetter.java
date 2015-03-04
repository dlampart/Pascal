package ch.unil.genescore.prioritization;

import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.gene.GenomicElementComparedByEnd;

import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;

/** wraps ConnectionTreeGetter and applies pruning to the returned Trees */

public class PrunedConnectionTreeGetter implements ConnectionTreeGetter {

	ConnectionTreeGetter myConnectionTreeGetter_ = null;
	int ext_;
	@Override
	public void loadNetworkData() {
		myConnectionTreeGetter_.loadNetworkData();

	}

	/**/
	@Override
	public TreeMultimap<Double, Gene> returnReverseOrderedConnectionTree(
			Gene currentGene) {
		TreeMultimap<Double, Gene> tmp = myConnectionTreeGetter_
				.returnReverseOrderedConnectionTree(currentGene);
		return (pruneTreeByDistance(tmp));

	}

	@Override
	public void setGenes(Collection<Gene> Genes) {
		myConnectionTreeGetter_.setGenes(Genes);

	}

	public void setConnectionTreeGetter(
			ConnectionTreeGetter myConnectionTreeGetter) {
		myConnectionTreeGetter_ = myConnectionTreeGetter;
	}

	//TODO: hack: made pruneTreeByDistance public so as to use it like "static" method in PrioritizationMain
	
	/**
	 * input:TreeMultimap<Double,Gene> sortedByValue : Genes that are connected
	 * to to the current running gene; reverse sorted!!!! Note: starting element
	 * is the gene itself value:1.1 (all other genes have maximally value 1)
	 * output: TreeMultimap<Double,Gene> prunedSortedByValue : same tree but
	 * positions are pruned based on position process: iterate over
	 * sortedByValue (iterator descending by value) and add to
	 * prunedSortedByValue unless tree contains already gene within the gene
	 * region param: ext: add ext to one gene before comparing overlap
	 * */
	public TreeMultimap<Double, Gene> pruneTreeByDistance(
			TreeMultimap<Double, Gene> sortedByValue) {


		TreeSet<Gene> startSortedSet = new TreeSet<Gene>();
		TreeSet<GenomicElementComparedByEnd> endSortedSet = new TreeSet<GenomicElementComparedByEnd>();
		TreeMultimap<Double, Gene> prunedSortedByValue = TreeMultimap.create(
				Ordering.natural().reverse(), Ordering.natural());

		boolean isFirst = true;
		for (Map.Entry<Double, Gene> entry : sortedByValue.entries()) {
			Gene myCurrentGene = entry.getValue();
			GenomicElementComparedByEnd myCurrentGeneComparedByEnd = new GenomicElementComparedByEnd(myCurrentGene);
			boolean toAdd = isFirst
				||	doesntOverlap(startSortedSet,myCurrentGene,endSortedSet,myCurrentGeneComparedByEnd,ext_);
			isFirst = false;
			if(toAdd){				
				prunedSortedByValue.put(entry.getKey(), entry.getValue());
				startSortedSet.add(myCurrentGene);
				endSortedSet.add(myCurrentGeneComparedByEnd);
			}
		}
		return (prunedSortedByValue);
	}

private boolean doesntOverlap(TreeSet<Gene> startSortedSet,Gene gene,
		TreeSet<GenomicElementComparedByEnd> endSortedSet,GenomicElementComparedByEnd geneComparedByEnd,int ext){
	boolean out = (!overlaps(startSortedSet,gene,ext)
		 	&& !overlaps(endSortedSet,geneComparedByEnd,ext));
	return out;
	
}
	
private boolean checkFloor(TreeSet<Gene> set, Gene el, int ext){
	boolean out = (set.floor(el) == null)
			|| !el.extensionPartiallyOverlapsElement(set.floor(el), -ext, ext);			
	return out;
}

private boolean checkFloor(TreeSet<GenomicElementComparedByEnd> set, GenomicElementComparedByEnd el, int ext){
	boolean out = (set.floor(el) == null)
			|| !el.extensionPartiallyOverlapsElement(set.floor(el), -ext, ext);			
	return out;
}

private boolean checkCeiling(TreeSet<GenomicElementComparedByEnd> set, GenomicElementComparedByEnd el, int ext){
	boolean out = (set.ceiling(el) == null) 
			|| !el.extensionPartiallyOverlapsElement(set.ceiling(el), -ext, ext);			
	return out;
}


private boolean checkCeiling(TreeSet<Gene> set, Gene el, int ext){
	boolean out = (set.ceiling(el) == null) 
			|| !el.extensionPartiallyOverlapsElement(set.ceiling(el), -ext, ext);			
	return out;
}


private boolean overlaps(TreeSet<Gene> set,Gene el,int ext){
	boolean out = !checkFloor(set,el,ext)
			|| !checkCeiling(set,el,ext);
			return out;
}
private boolean overlaps(TreeSet<GenomicElementComparedByEnd> set,GenomicElementComparedByEnd el,int ext){
	boolean out = !checkFloor(set,el,ext)
			|| !checkCeiling(set,el,ext);
			return out;
}
public void setExt(int ext){
	ext_=ext;
}

}
