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

import java.util.Iterator;
import java.util.LinkedList;

import ch.unil.genescore.gene.Chromosome;
import ch.unil.genescore.gene.GenomicElement;
/**this is a recursive DataStructure*/
public class OverlappedGenomicElement implements Comparable<OverlappedGenomicElement> {

	GenomicElement mainElement_ = null;
	LinkedList<OverlappedGenomicElement> overlappedElements_ = null;	
	
	public OverlappedGenomicElement(GenomicElement mainElement){
		mainElement_ = mainElement;
		overlappedElements_ = new LinkedList<OverlappedGenomicElement>();
	}
	
	public void  addToList(OverlappedGenomicElement element){
		overlappedElements_.add(element);
	}
	public GenomicElement getMainElement(){
		return mainElement_;
	}
	public LinkedList<GenomicElement> getAllOverlappedElements(){		
		LinkedList<GenomicElement> elementList = new LinkedList<GenomicElement>(); 
		for (OverlappedGenomicElement overlappedElement :  overlappedElements_){
			elementList.add(overlappedElement.getMainElement());
		}
		return elementList;
	}
	public LinkedList<OverlappedGenomicElement> addAllOverlappedElements(LinkedList<OverlappedGenomicElement> priorList){		
		LinkedList<OverlappedGenomicElement> elementList = new LinkedList<OverlappedGenomicElement>(); 
		for (OverlappedGenomicElement overlappedElement :  overlappedElements_){
			elementList.add(overlappedElement);
		}
		priorList.addAll(elementList);
		return elementList;
	}
	public void filterMembers(String id){
		Iterator<OverlappedGenomicElement> iter = overlappedElements_.iterator();
		OverlappedGenomicElement current = null;
		String currentId = null;
		while(iter.hasNext()){
			current=iter.next();
			currentId = current.mainElement_.id_;
			if (!currentId.equals(id)){
				iter.remove();
			}
		}
	}
	
	public LinkedList<GenomicElement> getRecursiveOverlappedElementAtLevelN(int n){
		LinkedList<GenomicElement> outList= new LinkedList<GenomicElement>();
		if (n>0){
			for (OverlappedGenomicElement overlappedElement : overlappedElements_){				
				outList.addAll(overlappedElement.getRecursiveOverlappedElementAtLevelN(n-1));
			}
		}
		if (n==0){
			outList=getAllOverlappedElements();
		}
		return outList;
	}
	//TODO: potentially interface properly
	// now just double def in both classes.
	/** Check whether main element completely overlaps second main element
    */
    public boolean completelyOverlapsElement(OverlappedGenomicElement other, int start_mod, int end_mod) {          	
    	return mainElement_.completelyOverlapsElement(other.getMainElement(), start_mod, end_mod);
    }
	/** Check whether main element partially overlaps second main element
     */
     public boolean partiallyOverlapsElement(OverlappedGenomicElement other) {          	
     	return mainElement_.partiallyOverlapsElement(other.getMainElement());
     }
     
     /** Check whether main element partially overlaps second main element
      * * conceptually move first element borders by start_mod and end_mod first.start+start_mod, first.end+end_mod
      */
      public boolean extensionPartiallyOverlapsElement(OverlappedGenomicElement other, int start_mod, int end_mod) {          	
      	return mainElement_.extensionPartiallyOverlapsElement(other.getMainElement(), start_mod, end_mod);
      }
      
      /** Check whether element is completely downStream of second
       * Compare first by chromosome
       * then check wether first.start + start_mod is higher than sec.end  
       * conceptually first move  first.start by start_mod*/
    public boolean isPastElement(OverlappedGenomicElement other, int start_mod) {      
    	return mainElement_.isPastElement(other.getMainElement(),  start_mod);    	    	
    }
    /** Compare first by chromosome, second by start position, third by id */
    @Override
    public int compareTo(OverlappedGenomicElement other) { 
    	return mainElement_.compareTo(other.getMainElement());
    }
    
    /** 
     * Override hashCode to be consistent with equals() and compareTo().
     * Automatically generated using Eclipse
     */
	@Override
	public int hashCode() {
		return mainElement_.hashCode();
	}
	/** 
	 * Override equals comparing by chromosome, start and id (consistent with compareTo() and hashCode()).
	 * Automatically generated using Eclipse.
	 */
	@Override
	public boolean equals(Object obj) {
		return mainElement_.equals(obj);
	
	}
}
	
