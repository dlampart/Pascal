package ch.unil.genescore.vegas;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Settings;

public class AllOverlappedElements{
 
	TreeSet<OverlappedGenomicElement> myLargeList_ = null; 
	OverlappedGenomicElementStream myStream_ = null;
	int mod_start_ = 0;
	int mod_end_ = 0;
	
	public AllOverlappedElements(OverlappedGenomicElementStream myStream, Collection<? extends GenomicElement> elementList){
		// has to be base list
		myLargeList_ = new TreeSet<OverlappedGenomicElement>();
		for (GenomicElement currentElement : elementList){
			myLargeList_.add(new OverlappedGenomicElement(currentElement));
		}
		myStream_ = myStream;
		
	}
	
	public AllOverlappedElements(OverlappedGenomicElementStream myStream, OverlappedGenomicElementStream myRefStream){
			//Collection<? extends GenomicElement> elements){
		if (!myStream.streamOpen() || !myRefStream.streamOpen()){
			throw new RuntimeException("one of the streams not open");
		}
		myLargeList_ = new TreeSet<OverlappedGenomicElement>();
		OverlappedGenomicElement currentLoaded = null;
		while (myRefStream.streamOpen()){
			currentLoaded=myRefStream.getNextAsOverlappedGenomicElement();
			myLargeList_.add(currentLoaded);
		}
		myStream_ = myStream;
	
	}
	public void fillTreeSet(){	
		ArrayDeque<OverlappedGenomicElement> myDeque = new ArrayDeque<OverlappedGenomicElement>();
		OverlappedGenomicElement currentStreamedElement= null;
		
		Iterator<OverlappedGenomicElement>  myIt = myLargeList_.iterator();		
		while(true){
			//load StreamedElement
			if (myStream_.streamOpen())
				currentStreamedElement = myStream_.getNextAsOverlappedGenomicElement();	
			else 
				break;
			//load TreeElements
			while (!myDeque.isEmpty() && currentStreamedElement.isPastElement(myDeque.getLast(), -mod_end_)){
				myDeque.removeLast();
			}			
			while (myIt.hasNext() && (myDeque.isEmpty() || !myDeque.getFirst().isPastElement(currentStreamedElement, mod_start_))){
				myDeque.addFirst(myIt.next());
			}
			for (OverlappedGenomicElement dequeElement : myDeque){
				if (dequeElement.completelyOverlapsElement(currentStreamedElement, mod_start_, mod_end_)){
					dequeElement.addToList(currentStreamedElement);
				}
			}		
		}
		
	}
	public TreeSet<OverlappedGenomicElement> getMyLargeList(){
		return  myLargeList_;
	}
	public GenewiseSnpWeights getAsGenewiseSnpWeightDataStructure(double weight, int recursionLevel){
		HashMap<String, SnpWeightMap> snpWeightsForEachGene = new HashMap<String, SnpWeightMap>();
		LinkedList<GenomicElement> myGenomicElements = new LinkedList<GenomicElement>(); 
		for (OverlappedGenomicElement currentElement : myLargeList_){			
			myGenomicElements=currentElement.getRecursiveOverlappedElementAtLevelN(recursionLevel);			
			String ElementId = currentElement.getMainElement().id_;
			snpWeightsForEachGene.put(ElementId, new SnpWeightMap(myGenomicElements, weight));			 
		}
		if (Settings.runPathwayAnalysis_)
			return new MetaGenewiseSnpWeights(snpWeightsForEachGene);				
		else
			return new GenewiseSnpWeights(snpWeightsForEachGene);				
	}
	/**removes elements in lists is custom made for gene_symbol ids from encode-map*/
	public void filterTree(){
		HashMap<String,String> myIdMap = new HashMap<String, String>();
		Gene asGene = null;
		for(OverlappedGenomicElement el : myLargeList_){
			asGene = (Gene) el.mainElement_;
			myIdMap.put(asGene.id_, asGene.symbol_);			
		}
	
		String mappedId = null;
		String currentId = null;
		for(OverlappedGenomicElement el : myLargeList_){
			currentId=el.mainElement_.id_;
			if(myIdMap.containsKey(currentId)){
				mappedId=myIdMap.get(currentId);
			}
			else {mappedId = "";}
			el.filterMembers(mappedId);
		}
	}
	public void setMods(int mod_start, int mod_end){
		mod_start_ = mod_start;
		mod_end_ = mod_end;
		
		
	}
	

	
		
}
