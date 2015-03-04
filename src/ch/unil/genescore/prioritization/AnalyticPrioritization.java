package ch.unil.genescore.prioritization;

import java.util.Map;
import java.util.TreeMap;

import com.google.common.collect.TreeMultimap;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;
import ch.unil.genescore.vegas.AnalyticVegas;
import ch.unil.genescore.vegas.DistributionMethods;

public class AnalyticPrioritization extends AnalyticVegas{
	
	TreeMultimap<Double,Gene> prunedValueGeneMap_ = null;	
	
	public AnalyticPrioritization(TreeMultimap<Double,Gene> prunedValueGeneMap){
		prunedValueGeneMap_=prunedValueGeneMap;
    }
	
	public boolean computeScore() {
		
		reinitialize();				
		
		processData();
		computeAnalyticVegasPvalues();
		return status_ != Status.DAVIES_FAIL_FAREBROTHER_FAIL;
	}
	
	//public void removeGene
	protected void removeFirst(){		
		prunedValueGeneMap_.removeAll(1.1);//remove gene itself.
	}
	
	private void processData(){				
		
		removeFirst();
		testStatisticReal_=0;
		lambda_ = new double[prunedValueGeneMap_.entries().size()];
		int count=0;
		for(Map.Entry<Double, Gene> entry : prunedValueGeneMap_.entries()){
			
		lambda_[count]=entry.getKey();		
		double currentNormalizedScore=entry.getValue().getNormalizedScore(); 
		testStatisticReal_ += lambda_[count]* DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(currentNormalizedScore);
		count++;
		}
	}
	
	
	public double[] computeLambda(){
		return lambda_;
	}
	
	/** Get string representation of result that will be written to the output file */
	  public String getResultsAsString() {
		  
		  String line = "\t" + lambda_.length + "\t" + Utils.toStringScientific10(geneScore_) + "\t" + getStatusString();		 
	        return line;
	    }

	    /** Get a header line corresponding to getResultsAsString() */
	    public String getResultsAsStringHeader() {

	        //String header = "\tpvalue\tnumSnps";
	    	String header = "\tnumNeighbours\tpvalue\tStatus";
	  //      if (Settings.writeDetailedOutput_)
	   //         header += "\tavgSnpCorrelation\tstatus";
	        return header;
	    }
	
	// ----------------------------------------------------------------------------

	/** Get output to be printed on console after computing score for a gene */
	public String getConsoleOutput() { 
		
		if (Settings.writeDetailedOutput_) {
			return getStatusString();
		
		}
		else return "";
	}
}
