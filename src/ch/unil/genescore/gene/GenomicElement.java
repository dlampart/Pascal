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
package ch.unil.genescore.gene;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import ch.unil.genescore.main.Pascal;


/**
 * A position on the genome (e.g., SNP)
 */
public class GenomicElement implements Comparable<GenomicElement> {

	/** The ID of the element */
	public String id_ = null;
	
	/** The chromosome (chr1, ..., chr22, chrX, chrY, chrM) */
	public String chr_ = "none";
	/** The start position */
	public int start_ = -1;
	/** The end position */
	public int end_ = -1;
	/** The strand (true for +, false for -) */
	public boolean posStrand_ = true;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GenomicElement(String id) {

		id_ = id;
	}

	
	// ----------------------------------------------------------------------------
//TODO:@Daniel changed compareTo for chromosome to lexical compare. (the reason being that the bedops tool sorts things the same way).
	//question: do we need to generate the hashfunction again?
	
	/** Compare first by chromosome, second by start position, third by end position, fourth by id */
    @Override
    public int compareTo(GenomicElement other) {      

    	// Compare chromosome
    	int chrDiff = chr_.compareTo(other.chr_);
    	if (chrDiff != 0)
    		return chrDiff;
    	
    	// Compare start position
    	int startDiff = start_ - other.start_;
    	if (startDiff != 0)
    		return startDiff;
    	int endDiff = end_ - other.end_;
    	if (endDiff != 0)
    		return endDiff;
    	// Compare id
        int idDiff = id_.compareToIgnoreCase(other.id_);
        return idDiff;
    }
    
    /** Check whether element completely overlaps second
      */
    public boolean completelyOverlapsElement(GenomicElement other,int start_mod, int end_mod) {      

    	
    	// Compare chromosome
    	int chrDiff = chr_.compareTo(other.chr_);
    	if (chrDiff != 0)
    		return false;
    	
    	if (start_ + start_mod <= other.start_ && end_ + end_mod >= other.end_ ) 
    		return true;
    	else
    		return false;
    }
    
    /** Check whether element partially overlaps second
     */
   public boolean partiallyOverlapsElement(GenomicElement other) {      

   	   	// Compare chromosome
   	int chrDiff = chr_.compareTo(other.chr_);
   	if (chrDiff != 0)
   		return false;
   	
   	if (start_ <= other.start_ && end_ >= other.start_ ) 
   		return true;
   	else if (start_ <= other.end_ && end_ >= other.end_ ) 
   		return true;
   	if (other.start_ <= start_ && other.end_ >= start_ ) 
   		return true;
   	else if (other.start_ <= end_ && other.end_ >= end_ ) 
   		return true;
   	else
   		return false;
   }
    
   /** Check whether element partially overlaps second
    * conceptually move first element borders by start_mod and end_mod first.start+start_mod, first.end+end_mod
    */
  public boolean extensionPartiallyOverlapsElement(GenomicElement other, int start_mod, int end_mod) {      
	  int start = start_ + start_mod;
	  int end = end_ + end_mod;
  	   	// Compare chromosome
  	int chrDiff = chr_.compareTo(other.chr_);
  	if (chrDiff != 0)
  		return false;
  	
  	if (start <= other.start_ && end >= other.start_ ) 
  		return true;
  	else if (start <= other.end_ && end >= other.end_ ) 
  		return true;
  	if (other.start_ <= start && other.end_ >= start ) 
  		return true;
  	else if (other.start_ <= end && other.end_ >= end ) 
  		return true;
  	else
  		return false;
  }
   
    /** Check whether element is completely downStream of second
     * Compare first by chromosome
     * then check wether first.start + start_mod is higher than sec.end  
     * conceptually first move  first.start by start_mod*/
    public boolean isPastElement(GenomicElement other, int start_mod) {      

    	// Compare chromosome
    	int chrDiff = chr_.compareTo(other.chr_);   	
    	if (chrDiff > 0)
    		return true;
    	if (chrDiff < 0)
    		return false;    	
    	if (start_ + start_mod >= other.end_)
    		return true;
    	else
    		return false;
    	
    }
	// ----------------------------------------------------------------------------

    /** 
     * Override hashCode to be consistent with equals() and compareTo().
     * Automatically generated using Eclipse
     */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chr_ == null) ? 0 : chr_.hashCode());
		result = prime * result + ((id_ == null) ? 0 : id_.hashCode());
		result = prime * result + start_;
		return result;
	}


	// ----------------------------------------------------------------------------

	/** 
	 * Override equals comparing by chromosome, start and id (consistent with compareTo() and hashCode()).
	 * Automatically generated using Eclipse.
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GenomicElement other = (GenomicElement) obj;
		if (chr_ == null) {
			if (other.chr_ != null)
				return false;
		} else if (!chr_.equals(other.chr_))
			return false;
		if (id_ == null) {
			if (other.id_ != null)
				return false;
		} else if (!id_.equals(other.id_))
			return false;
		if (start_ != other.start_)
			return false;
		return true;
	}
	
	/** Get the elements that are in the window around the given Elements */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public ArrayList<GenomicElement> findElements(Genome allElements) {
		
		// Define the window around the gene (depends on orientation)
		int start;
		int end;
		if (posStrand_) {
			start = start_ - Pascal.set.geneWindowUpstream_;
			end = end_ + Pascal.set.geneWindowDownstream_;
		} else {
			start = start_ - Pascal.set.geneWindowDownstream_;
			end = end_ + Pascal.set.geneWindowUpstream_;
		}
		
		// Get the set of snps assigned to this gene as an array
		ArrayList<GenomicElement> overlappedElements = (ArrayList) allElements.getElementsIn(chr_, start, end);
		return overlappedElements;
	}
	
	//@SuppressWarnings({ "rawtypes", "unchecked" })
	public HashSet<GenomicElement> findElementsHash(Genome allElements) {
		
		HashSet<GenomicElement> hashset = new HashSet<GenomicElement>();
		ArrayList<GenomicElement> elements = findElements(allElements);
		for(GenomicElement el : elements)
			hashset.add(el);
		return(hashset);
	}


	// ----------------------------------------------------------------------------

	/** Get a string representation in UCSC-BED format (0-1 coding)*/
	public String UCSCbedFormatString() {
		
		int start0 = start_ - 1 ;
		return chr_ + "\t" + start0 + "\t" + end_ + "\t" + id_;
	}


	//	----------------------------------------------------------------------------

	/** 
	 * Write this genomic element to the binary file. Format is:
	 * 1. id_
	 * 2. chr_
	 * 3. start_
	 * 4. end_
	 * 5. posStrand_
	 * @throws IOException 
	 */
	public void writePos(DataOutputStream os) throws IOException {
		
		// NOTE: ALSO CHANGE readBinary() IF YOU CHANGE THIS
		os.writeUTF(id_);
		os.writeUTF(chr_);
		os.writeInt(start_);
		os.writeInt(end_);
		os.writeBoolean(posStrand_);
	}


	//	----------------------------------------------------------------------------

	/** 
	 * Read this genomic element from a binary file. Format is:
	 * 1. id_
	 * 2. chr_
	 * 3. start_
	 * 4. end_
	 * 5. posStrand_
	 * @throws IOException 
	 */
	public void readPos(DataInputStream is) throws IOException {
		
		// id_ is already read 
		chr_ = is.readUTF();
		start_ = is.readInt();
		end_ = is.readInt();
		posStrand_ = is.readBoolean();
		
	}

	
	//	----------------------------------------------------------------------------

	/** Return '+' / '-' for positive / negative strand, respectively */
	public String getPosStrandStr() {
		
		return posStrand_ ? "+" : "-";
	}
	
	
	//	----------------------------------------------------------------------------

	/** Return true for '+', false for '-' and error otherwise */
	static public boolean isPosStrand(String ch) {
	
		if (ch.equals("+"))
			return true;
		else if (ch.equals("-"))
			return false;
		else 
			throw new RuntimeException("Strand has to be '+' or '-'");
	}

	public void copyPosition(GenomicElement el){
		chr_ = el.getChr();
		start_ = el.getStart();
		end_ = el.getEnd();
		//if (curChr != null || curStart != -1 || curEnd != -1){
			//	if (!chr_.equals(curChr) || start_ != curStart || end_ != curEnd){
				//	throw new RuntimeException("snp seems to have been set before to another value");
				//}
		//}
		posStrand_ = el.getPosStrand();
	}
	
	
	// ============================================================================
	// GETTERS AND SETTERS

	public String getId() { return id_; }
	public int getStart() { return start_; }
	private int getEnd() { return end_; }
	private String getChr() { return chr_; }
	public int getChrInt() { return Integer.parseInt(chr_.substring(3)); }
	private boolean getPosStrand() { return posStrand_;}
	
	//public void setPosition(String chr, int start) {
	//	chr_ = chr;
	//	start_ = start;
	//	end_ = start;
	//}

	public void setPosition(String chr, int start, int end, boolean posStrand) {
		chr_ = chr;
		start_ = start;
		end_ = end;
		posStrand_ = posStrand;
	}

}
