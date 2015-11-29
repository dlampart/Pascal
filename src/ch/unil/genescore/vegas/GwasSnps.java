/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.vegas;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;

import ch.unil.genescore.gene.GeneAnnotationGencode;
import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileParser;

public class GwasSnps {


	private HashMap<String, Snp>  snpsInList_ = null;			
	boolean hasHeader_;
	//TODO: sort out header situation
	
	
	public void loadSnpPvalZval(File file){
        
        LinkedHashMap<String, Snp> snps = new LinkedHashMap<String, Snp>();
        // Open the file
        FileParser parser = new FileParser(Pascal.log, file);
       
        String[] nextLine=null;
        int numNA=0;
        if (hasHeader_)
        	nextLine = parser.readLine();
        // Check that all subsequent lines have the same number of cols
        nextLine = parser.readLine();
        int numCols = nextLine.length;
        while (nextLine != null) {
            // Check number of columns
            if (nextLine.length != numCols)
                throw new RuntimeException("Line " + parser.getLineCounter() + ": Inconsistent number of columns");
            
            // Parse  p-value and z-value
           
            // TODO: check file format
            String id = nextLine[0];
          //  zval = Double.parseDouble(nextLine[4]); 
           // pval = Double.parseDouble(nextLine[5]); 
            if(id.equals("rs11177594")){
            	System.out.println("");
            }
            if (nextLine[5].equalsIgnoreCase("NA") || nextLine[5].equalsIgnoreCase("NaN") || nextLine[5].equalsIgnoreCase("")
            		||
            		nextLine[4].equalsIgnoreCase("NA") || nextLine[4].equalsIgnoreCase("NaN") || nextLine[4].equalsIgnoreCase("")) {
				numNA++;
				nextLine = parser.readLine();
				continue;
			}
            double zval = 0;
            try{
				zval = Double.parseDouble(nextLine[4]);
				if (Double.isNaN(zval)){					
					String msg = "snp " + nextLine[0] + " couldn't be processed properly. (Maybe it's a NaN with additional white-spaces?)";
					Pascal.warning(msg);
					nextLine = parser.readLine();
					continue;
				}
				if (Double.isInfinite(zval)){					
					String msg = "snp " + nextLine[0] + " was read as infinity. Whas removed set to reasonable value yourself";
					Pascal.warning(msg);
					nextLine = parser.readLine();
					continue;
				}
				
				
			}catch(NumberFormatException e){
				String msg = "snp " + nextLine[0] + " has no valid z-value (use NA or NaN for omitted values).";
				Pascal.warning(msg);
				nextLine = parser.readLine();
				continue;				
			}        
            double pval=1;
            try{
				pval = Double.parseDouble(nextLine[5]);
			}catch(NumberFormatException e){
				String msg = "snp " + nextLine[0] + " has no valid p-value (use NA or NaN for omitted values).";
				Pascal.warning(msg);
				nextLine = parser.readLine();
				continue;				
			}            
            if (Double.isNaN(pval)){					
				String msg = "snp " + nextLine[0] + " couldn't be processed properly. (Maybe it's a NaN with additional white-spaces?)";
				Pascal.warning(msg);
				nextLine = parser.readLine();
				continue;
			}
            if (pval < 1e-300) {            
            	pval = 1e-300;
            	zval = Math.signum(zval)*(-1)*DistributionMethods.normalInverseCumulativeProbability(0.5e-300);
				
			}
            // this is random we don't really know if its true but we still need to impute the relative direction
            char majorAllele = Character.toUpperCase(nextLine[1].charAt(0));
            char minorAllele = Character.toUpperCase(nextLine[2].charAt(0));
            Snp newSnp =  new Snp(id, pval, zval);
            newSnp.setMajorAllele(majorAllele);
            newSnp.setMinorAllele(minorAllele);
            snps.put(id, newSnp);
            nextLine = parser.readLine();
        }
        parser.close();
        if (numNA > 0)
			Pascal.println("   - " + numNA + " SNPs were ignored because they have 'NA', 'NaN' or '' as p-value");
	
        snpsInList_ = snps;
    }	
	
	/** Load SNPs with p-values */
	public void loadSnpPvals(File file) {
		
		LinkedHashMap<String, Snp> snps = new LinkedHashMap<String, Snp>();
		FileParser parser = new FileParser(Pascal.log, file);		
		 String[] nextLine=null;
	        if (hasHeader_)
	        	nextLine = parser.readLine();
		// Open the file
		// Check that list is sorted
//		double previousPval = 0.0;
		
		// Read first line
		nextLine = parser.readLine();
		// Check that all subsequent lines have the same number of cols
		int numCols = nextLine.length;
		// Check that there are enough columns
		int pvalCol = Pascal.set.pvalCol_ - 1; // Here we start at 0
		if (pvalCol >= numCols)
			parser.error("Not enough columns, expected p-value in column specified in settings file snpPvalCol=" + Pascal.set.pvalCol_);
		
		// Snps that are listed multiple times in the input file
		LinkedHashMap<String, Integer> multiSnps = new LinkedHashMap<String, Integer>();
		// Snps with NA values
		int numNA = 0;
		// Snps with p-values below minPvalue
		int numBelowMin = 0;
		
		while (nextLine != null) {
			
			// Check number of columns
			if (nextLine.length != numCols)
				throw new RuntimeException("Line " + parser.getLineCounter() + ": Inconsistent number of columns");

			// Parse the p-value
			if (nextLine[pvalCol].equalsIgnoreCase("NA") || nextLine[pvalCol].equalsIgnoreCase("NaN") || nextLine[pvalCol].equalsIgnoreCase("")) {
				numNA++;
				nextLine = parser.readLine();
				continue;
			}
			double pvals=1;
			try{
				pvals = Double.parseDouble(nextLine[pvalCol]);
			}catch(NumberFormatException e){
				String msg = "snp " + nextLine[0] + "has no valid p-value (use NA or NaN for omitted values).";
				nextLine = parser.readLine();
				Pascal.warning(msg);
				continue;
				
			}
			 if (Double.isNaN(pvals)){					
					String msg = "snp" + nextLine[0] + "couldn't be processed properly. (Maybe it's a NaN with additional white-spaces?)";
					Pascal.warning(msg);
					continue;
			 }
			// Check that p-value is in [0,1]
			if (pvals < 0.0 || pvals > 1.0)
				parser.error("P-value not in [0,1]");
			// Some p-values become zero because they are smaller than the min double value
			if (pvals < 1e-300) {
				String msg = "snp " + nextLine[0] + "has p-value below 1E-300. Set to 1E-300.";
				Pascal.warning(msg);
				pvals = 1e-300;
				numBelowMin++;
			}
						
			// Add the SNP... Some pval files have the same snp listed multiple times
			String id = nextLine[0];
			
			// If it was already listed multiple times, increase it's count
			if (multiSnps.containsKey(id)) {
				multiSnps.put(id, multiSnps.get(id) + 1);
				
			// If it is already in the snp list, remove and add it to the ambigous set
			} else if (snps.containsKey(id)) {
				//throw new RuntimeException("The following SNP is listed twice: " + nextLine[0]);
				multiSnps.put(id, 1);
				snps.remove(id);
			
			// If it occurs the first time, add it (ideally, we should only have this case, I don't know why they list the same snp multiple times in some gwas)
			} else {
				snps.put(id, new Snp(id, pvals));
			}
			
			// Read next line
			nextLine = parser.readLine();
		}
		parser.close();
		
		if (multiSnps.size() > 0)
			Pascal.warning(multiSnps.size() + " SNPs were removed because they are listed multiple times in the p-value file");
		if (numNA > 0)
			Pascal.println("   - " + numNA + " SNPs were ignored because they have 'NA', 'NaN' or '' as p-value");
		if (numBelowMin > 0)
			Pascal.println("   - " + numBelowMin + " SNPs were set to the minimum allowed p-value (1e-300)");
		
		// Sort the list
		//Collections.sort(list_);
		snpsInList_ = snps;
	}

	/** Load the positions of the given snps */ 
	public void loadCodingSnps(File file, HashMap<String, Snp> snps) {
		
		// Open file
		FileParser parser = new FileParser(Pascal.log, file);
		
		while (true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			
			// Check number of columns
			if (nextLine.length != 2)
				parser.error("Line does not have same number of columns as first line");

			// Get the snp		
			Snp snp = snps.get(nextLine[0]);
			// Skip snps that are not part of the initial set
			if (snp == null)
				continue;

			// Get the gene
			String gene = GeneAnnotationGencode.removeEnsemblVersion(nextLine[1]);
			snp.addCoding(gene);
		}
		parser.close();
	}
	/** keep only snps form list that are in snpFilterFile */
	public void keepOnlySnpsInFilterListFile(File snpFilterFile){

		FileParser parser = new FileParser(Pascal.log, snpFilterFile);
		HashMap<String, Snp> filteredSnps = new HashMap<String, Snp>();
		String snpId = null;
		
		Snp presentSnp = null; 
		while (true){
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			snpId = nextLine[0];					
			presentSnp = snpsInList_.get(snpId);
			if (presentSnp != null){
				
				filteredSnps.put(snpId, presentSnp);
			}			
					
		}
		snpsInList_ = filteredSnps;
	}		
	public HashMap<String, Snp> getSnpsInList(){return snpsInList_;}
	
	public void setHeader(boolean hasHeader){
		hasHeader_=hasHeader;
	}
}

