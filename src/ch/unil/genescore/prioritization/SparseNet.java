package ch.unil.genescore.prioritization;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileParser;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;




public class SparseNet implements ConnectionTreeGetter {

	
	HashBasedTable<String, String, Double> network = HashBasedTable.create();
	//LinkedList<String>  firstGeneNameList = new LinkedList<String>();
	//LinkedList<String>  secondGeneNameList = new LinkedList<String>();
//	LinkedList<Double>  values = new LinkedList<Double>();
	
	private HashMap<String,Gene> geneList_ = new HashMap<String,Gene>();
	private FileParser parser_ = null;
	
	
	
	/**expects headerless file. tab-sep; first two rows: gene-symbol; third row double(interaction Strength)
	 * also assumes that each pair only occurs once.
	 * 
	 * */
	@Override
	public void loadNetworkData() {

		while (true) {
			String[] nextLine = parser_.readLine();
			
			if (nextLine == null)
				break;						
			if(nextLine.length != 3){throw new RuntimeException("Exception during parsing of network: not 3 columns in network file");}
			Double myVal = Double.valueOf(nextLine[2]);
			if(myVal>1){
				throw new RuntimeException("Exception during parsing of network: Value is larger than one. One is the allowed maximum. Scale network externally before trying again.");
			}
			if(!myVal.isNaN() & !myVal.isInfinite())
				network.put(nextLine[0], nextLine[1],myVal);				
		}		
		parser_.close();		
	}


	private Map<String,Double> getPairsForGene(String GeneName){
		
		//check uniqueness
		HashMap<String,Double> totset = new HashMap<String,Double>();
		Map<String,Double> colset = network.column(GeneName);
		int s1 = colset.size();
		totset.putAll(colset);
		Map<String,Double> rowset = network.row(GeneName);
		int s2 = rowset.size();
		totset.putAll(rowset);
		if(s1+s2!=totset.size()){
			throw new RuntimeException("Double Entry in parsed network or connections to itself. Gene "+ GeneName +" involved. remove double entries");
		}				
		return totset;
	}
	
	/**
	 * input Gene queryGene:
	 * output TreeMap<Double,Gene> sortedByValue
	 * returns a Tree of genes that are connected to currentGene in the network sorted by connection strength.
	 * 
	 * Mapping is done via gene symbol; 
	 * qu
	 * 
	 * bordercases
	 * returns empty tree if gene is not found in network; 
	 * 
	 * first entry: <1.1,queryGene>
	 * */
	@Override
	public TreeMultimap<Double, Gene> returnReverseOrderedConnectionTree(Gene queryGene) {
		TreeMultimap<Double,Gene> sortedByValue = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
		String geneName = queryGene.symbol_;		
		if(!geneList_.containsKey(geneName)){
			return(sortedByValue);
		}
		if(geneList_.get(geneName) !=  queryGene){
			throw new RuntimeException("Exception: query gene has same symbol as some gene in provided gene list but its not the same underlying gene object.");
			
		}
		sortedByValue.put(1.1, queryGene);		//has to be larger than 1 to be sorted completely right
		Map<String,Double> subset=getPairsForGene(geneName);
		for (Map.Entry<String,Double> entry : subset.entrySet()){
			if(geneList_.containsKey(entry.getKey())){
				Gene myGene = geneList_.get(entry.getKey());
				sortedByValue.put(entry.getValue(), myGene);
			}
		}		
		return sortedByValue;
	}
	
	public void setParser(String filename){
		parser_ = new FileParser(filename);
	}
	public void setParser(FileParser parser){
		parser_ = parser;
	}

	
	@Override
	public void setGenes(Collection<Gene> Genes) {
		for (Gene gene : Genes){
			geneList_.put(gene.symbol_,gene);
		}
		
	}
	
	

}
