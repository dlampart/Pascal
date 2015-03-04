package ch.unil.genescore.prioritization;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Utils;

import com.google.common.collect.TreeMultimap;

public class AnalyticPrioritizationNoRemoval extends AnalyticPrioritization {

	public AnalyticPrioritizationNoRemoval(
			TreeMultimap<Double, Gene> prunedValueGeneMap) {
		super(prunedValueGeneMap);
		// TODO Auto-generated constructor stub
	}

	//no removal
	protected void removeFirst(){		
	}
	
	/** Get string representation of result that will be written to the output file */
	  public String getResultsAsString() {
		  
		  String line = lambda_.length + "\t" + Utils.toStringScientific10(geneScore_) + "\t" + getStatusString();		 
	        return line;
	    }
	  
	  /** Get string representation of result that will be written to the output file */
	  public String getResultsAsStringHeader() {
		  
		  String line = "numNeighbours\tpvalue\tStatus";		 
	        return line;
	    }
}
