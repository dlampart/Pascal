package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.vegas.DistributionMethods;
import ch.unil.genescore.vegas.Snp;

/** 
 *  
 * interface:qtlGenePair  getEqtlGenePair(String GeneName);
 *  loadEqtlFile(String filename); 
*/
public class EqtlResults {

	
	HashMap<String, EqtlGenePair> eqtlResultList = null;
	boolean hasHeader_ = false;
	HashMap<String, ArrayList<Snp>> snpsInList_ =null;
	FileParser parser_ = null;
	
	public void loadEqtlFile(String filename){
		
		snpsInList_ = new LinkedHashMap<String, ArrayList<Snp>>();
        // Open the file
        parser_  = createParser(filename);
       
        String[] nextLine=null;
        int numNA=0;
        if (hasHeader_)
        	nextLine = parser_.readLine();
        // Check that all subsequent lines have the same number of cols
        nextLine = parser_.readLine();
        int numCols = nextLine.length;
        while (nextLine != null) {
            // Check number of columns
            if (nextLine.length != numCols)
                throw new RuntimeException("Line " + parser_.getLineCounter() + ": Inconsistent number of columns");
            String id = nextLine[0];
            //  zval = Double.parseDouble(nextLine[4]); 
             // pval = Double.parseDouble(nextLine[5]); 
              if (nextLine[5].equalsIgnoreCase("NA") || nextLine[5].equalsIgnoreCase("NaN") || nextLine[5].equalsIgnoreCase("")
              		||
              		nextLine[4].equalsIgnoreCase("NA") || nextLine[4].equalsIgnoreCase("NaN") || nextLine[4].equalsIgnoreCase("")) {
  				numNA++;
  				nextLine = parser_.readLine();
  				continue;
  			}
            double zval = 0;
            try{
				zval = Double.parseDouble(nextLine[4]);
				if (Double.isNaN(zval)){					
					String msg = "snp " + nextLine[0] + " couldn't be processed properly. (Maybe it's a NaN with additional white-spaces?)";
					Main.warning(msg);
					nextLine = parser_.readLine();
					continue;
				}
				if (Double.isInfinite(zval)){					
					String msg = "snp " + nextLine[0] + " was read as infinity. Whas removed set to reasonable value yourself";
					Main.warning(msg);
					nextLine = parser_.readLine();
					continue;
				}
				
			}catch(NumberFormatException e){
				String msg = "snp " + nextLine[0] + "has no valid z-value (use NA or NaN for omitted values).";
				Main.warning(msg);
				nextLine = parser_.readLine();
				continue;
				
			}  
            double pval=1;
            try{
				pval = Double.parseDouble(nextLine[5]);
				if (Double.isNaN(pval)){					
					String msg = "snp " + nextLine[0] + " couldn't be processed properly. (Maybe it's a NaN with additional white-spaces?)";
					Main.warning(msg);
					nextLine = parser_.readLine();
					continue;
				}
			}catch(NumberFormatException e){
				String msg = "snp " + nextLine[0] + "has no valid p-value (use NA or NaN for omitted values).";
				Main.warning(msg);
				nextLine = parser_.readLine();
				continue;
				
			}            
            if (pval < 1e-300) {            
            	pval = 1e-300;
            	zval = Math.signum(zval)*(-1)*DistributionMethods.normalInverseCumulativeProbability(0.5e-300);
				
			}
            char majorAllele = Character.toUpperCase(nextLine[1].charAt(0));
            char minorAllele = Character.toUpperCase(nextLine[2].charAt(0));
            Snp newSnp =  new Snp(id, pval, zval);
            newSnp.setMajorAllele(majorAllele);
            newSnp.setMinorAllele(minorAllele);            
            String geneName = nextLine[6];
            addSnpToSnpsInList(newSnp, geneName);            
            nextLine = parser_.readLine();
        }
        parser_.close();
        if (numNA > 0)
			Main.println("   - " + numNA + " SNPs were ignored because they have 'NA', 'NaN' or '' as p-value");	      
    }	        
	
	private void addSnpToSnpsInList(Snp snpToAdd, String geneName){
		
		if (!snpsInList_.containsKey(geneName))
			snpsInList_.put(geneName, new ArrayList<Snp>());
		
		snpsInList_.get(geneName).add(snpToAdd);
	}
	
	public boolean eqtlGenePairExists(String geneId){		
		return(snpsInList_.containsKey(geneId));
	}
	
	public EqtlGenePair  getEqtlGenePair(String geneId){
		if (!snpsInList_.containsKey(geneId))
			throw new RuntimeException("calling for gene that hasn't been loaded from eqtl results");
		ArrayList<Snp> snpAr = snpsInList_.get(geneId);		
		EqtlGenePair myPair = new EqtlGenePair(geneId, snpAr);		
		return(myPair);
	}
	
	public void  fillEqtlGenePair(String geneId, EqtlGenePair myPair){
		if (!snpsInList_.containsKey(geneId))
			throw new RuntimeException("calling for gene that hasn't been loaded from eqtl results");
		ArrayList<Snp> snpAr = snpsInList_.get(geneId);		
		myPair.setSnpList(snpAr);
		myPair.setGeneName(geneId);
		//return(myPair);
	}
	
	protected FileParser createParser(String filename){
		return(new FileParser(filename));
	}
	
}
