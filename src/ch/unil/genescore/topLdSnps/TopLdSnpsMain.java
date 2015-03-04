package ch.unil.genescore.topLdSnps;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeSet;

import no.uib.cipr.matrix.DenseMatrix;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.LinkageDisequilibrium;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.genescore.vegas.Snp;

public class TopLdSnpsMain extends ReferencePopulation {

	TreeSet<Snp> chr_ = null;
	
	public void run(){
		
		setSnps();
		initializeSnps();
		removeSnpsWithoutPos();
		setupSnpTree();
		goThroughTree();
		
	}
	private void removeSnpsWithoutPos(){
		removeSnpsWithoutPos(gwasSnps_);
		removeSnpsWithoutPos(relevantSnps_);
	}
	
	private void setupSnpTree(){
		
		chr_ =  new TreeSet<Snp>();
		for (Entry<String, Snp> entry : relevantSnps_.entrySet()){
			Snp currentSnp = entry.getValue();
			chr_.add(currentSnp);
		}		
	}
	
	private void goThroughTree(){
		
		
		ArrayList<String> stringsToWrite = new ArrayList<String>();
		for(Iterator<Snp> it = chr_.iterator(); it.hasNext();) {
			  //Integer key = entry.getKey();
			  Snp snp = it.next();
			  updateLoadedGenotypes(snp);
			  ArrayList<Snp> SnpsWithinRegion=getSnpsWithGenotypeWithinRegion(snp);			  
			  if(SnpsWithinRegion.size()==0){
				  System.out.println(snp.id_ + "no snps in region");
				  continue;
			  }			  
			  ArrayList<Snp> singleSnp = new ArrayList<Snp>();
			  singleSnp.add(snp);			  
			  DenseMatrix mat = LinkageDisequilibrium.computeCrossCorrelationMatrixMTJ(singleSnp, SnpsWithinRegion);
			  int ind= whichMax(mat);
			  String str = snp.id_ + "\t" + SnpsWithinRegion.get(ind).id_ + "\t" + mat.get(0,ind);
			  stringsToWrite.add(str);
			  
		}
		writeToFile(stringsToWrite,Settings.topLdSnpOutFile_);
		
	}
	private void writeToFile(ArrayList<String> strs, String filename){
		FileExport writer = new FileExport(filename);
		for (String str : strs){
			writer.println(str);			
		}
		writer.close();
	}
	
	private int whichMax(DenseMatrix mat){
		double val = 0;
		int index = -1;
		for (int i=0; i<mat.numColumns();i++){
			double curVal=Math.abs(mat.get(0,i));
			if (curVal>val){
				index=i;
				val=curVal;
			}
		}
		return(index);
	}
	
	private  ArrayList<Snp> getSnpsWithGenotypeWithinRegion(GenomicElement snp){
				
		ArrayList<Snp> myAr = new ArrayList<Snp>();
		for (Snp snpWithGenotype : loadedGenotypes_.getSnpsWithGenotypes()){
			if (snp.extensionPartiallyOverlapsElement(snpWithGenotype, -1000000,1000000)){
				if (gwasSnps_.containsKey(snpWithGenotype.id_)){
					myAr.add(snpWithGenotype);
				}
			}
		}
	 return myAr;
	}
		
	
	
	private void setSnps(){
		
		String list1Name = Settings.firstList_;
		String list2Name = Settings.secondList_;
		
		//FileParser parser2= new FileParser(list2Name);
		ArrayList<String> gwasList =  getSnpIds(list2Name);	
		ArrayList<String> relevantList =  getSnpIds(list1Name);	
		setGwasSnps(gwasList);
		setRelevantSnps(relevantList);	
		
		
		
	}	
	private void removeSnpsWithoutPos(HashMap<String,Snp> myMap){
		
		for(Iterator<Map.Entry<String,Snp>> it = myMap.entrySet().iterator(); it.hasNext();) {
			
			Map.Entry<String,Snp> entry = it.next();
			Snp mySnp = entry.getValue();
			if(mySnp.chr_==null){
				it.remove();
			}
			//String key = entry.getKey();
			  //Snp snp = (Snp) entry.getValue();
			 // entry.
			  
		}
	}
	
	private void setGwasSnps(ArrayList<String> snpList){
		
		
		 gwasSnps_ = new HashMap<String, Snp>();
		 for (int i=0; i<snpList.size(); i++){
			 String snpId=snpList.get(i);			 
			 Snp mySnp = new Snp(snpId);
			 gwasSnps_.put(snpId, mySnp);
		 }
		 
	}
	private void setRelevantSnps(ArrayList<String> snpList){				
				
			if (gwasSnps_==null){
				throw new RuntimeException("gwasSnps_ not yet set");
			}
				
			 relevantSnps_ = new HashMap<String, Snp>();
			 for (int i=0; i<snpList.size(); i++){
				 String snpId=snpList.get(i);			 
				 if(gwasSnps_.containsKey(snpId))
					 continue;
				 Snp mySnp = new Snp(snpId);				 
				 relevantSnps_.put(snpId, mySnp);
			 }		 
	}
	
	
	
	private ArrayList<String> getSnpIds(String listName){
		FileParser parser= new FileParser(listName);
		ArrayList<String> snpNamesList = new  ArrayList<String>();
			
			String[] line = parser.readLine();
			while (line!=null){
				snpNamesList.add(line[0]);
				line = parser.readLine();
			}
			parser.close();
			return (snpNamesList);
	}
	

}
	
	