/*
Copyright (c) 2013 Daniel Marbach

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper available at:
http://networkinference.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
package ch.unil.genescore.gene;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.Snp;


/**
 * A set of genomic elements organized by chromosome.
 * Note: In the current implementation we cannot have two distinct elements on the same
 * chromosome with the same start_, keep this in mind... (see comment in Chromosome)
 */
public class Genome {

	// TODO make this an arraylist, access chromosomes by their index throughout the project. Hm, maybe it's not worth changing
	/** The chromosomes (chr1, ..., chr22, chrX, chrY, chrM) */
	private LinkedHashMap<String, Chromosome> chromosomes_ = null;

	/** The total number of elements */
	private int numElements_ = 0;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Genome() {

		// initialize
		chromosomes_ = new LinkedHashMap<String, Chromosome>();
		numElements_ = 0;
		
		for (int i=1; i<=22; i++)
			chromosomes_.put("chr" + i, new Chromosome());
		if (!Settings.ignoreAllosomes_) {
			chromosomes_.put("chrX", new Chromosome());
			chromosomes_.put("chrY", new Chromosome());
		}
		if (!Settings.ignoreChrM_)
			chromosomes_.put("chrM", new Chromosome());
	}

	
	/** Constructor */
	@SuppressWarnings("rawtypes")
	public Genome(Collection elements) {
		
		this();
		addElements(elements);
	}

	// ----------------------------------------------------------------------------

	/** Add the elements of the specified chromosome to the genome (add all elements if chromosome is null or empty) */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void addElements(Collection elements) {
		String chromosome=Settings.chromosome_;
		boolean addAll = (chromosome == null) || (chromosome.length() == 0);
		Iterator<GenomicElement> iter = elements.iterator();
		
		while (iter.hasNext()) {
			GenomicElement nextElement = iter.next();
			
			// Skip if not specified chromosome
			if (nextElement.chr_ == "none" || 
					(!addAll && !nextElement.chr_.equals(chromosome)))
				continue;
			
			addElement(nextElement);
		}
		
		// Count number of elements
		if (addAll)
			numElements_ = elements.size();
		else
			numElements_ = chromosomes_.get(chromosome).getNumElements();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Add a genomic element */
	public void addElement(GenomicElement element) {

		// Add element to the corresponding chromosome
		chromosomes_.get(element.chr_).addElement(element);
	}

	
	// ----------------------------------------------------------------------------

	/** Get the elements in the given window */
	public ArrayList<GenomicElement> getElementsIn(String chr, int start, int end) {
		
		return chromosomes_.get(chr).getElementsIn(start, end);
	}

	
	// ----------------------------------------------------------------------------

	/** Get all elements of the given chromosome */
	public Collection<GenomicElement> getElements(String chr) {
		
		return chromosomes_.get(chr).getElements();
	}

	
	// ----------------------------------------------------------------------------

	/** Get all elements hashed by id and sorted by chromosome and start position */
	public LinkedHashMap<String, GenomicElement> asLinkedHashMap() {
		
		LinkedHashMap<String, GenomicElement> map = new LinkedHashMap<String, GenomicElement>();
		for (Chromosome chr : chromosomes_.values()) {
			for (GenomicElement el : chr.getElements())
				map.put(el.id_, el);
		}
		return map;
	}

	
	// ----------------------------------------------------------------------------

	/** Get the nearest element (in either direction) */
	public GenomicElement getNearestElement(String chr, int pos) {
		
		return chromosomes_.get(chr).getNearestElement(pos);
	}

	
	// ----------------------------------------------------------------------------

	/** Write a BED file with all the elements */
	public void writeBedFile(String filename) {
		
		FileExport writer = new FileExport(filename);
		
		for (Chromosome chr : chromosomes_.values()) {
			for (GenomicElement el : chr.getElements()) {
				String nextLine = el.UCSCbedFormatString();
				writer.println(nextLine);
			}
		}
		writer.close();
	}
	public void writeTpedPlinkFile(String filename) {		
		String filenameTped = filename + ".tped";
		FileExport writer = new FileExport(filenameTped);
		
		for (Chromosome chr : chromosomes_.values()) {
			for (GenomicElement el : chr.getElements()) {
				if(!(el instanceof Snp)){//TODO : replace with exception
					System.out.println("genome not filled only with snps.");
					writer.close();
					return;
				}										
				String nextLine = ((Snp) el).tpedString();
				writer.println(nextLine);
			}
		}
		writer.close();		
	}
	public void writePseudoTfamPlinkFile(String filename) {		
		String filenameTfam = filename + ".tfam";
		FileExport writer = new FileExport(filenameTfam);
		
		for (Chromosome chr : chromosomes_.values()) {	
			if (chr.getNumElements() == 0){
				continue;
			}
			GenomicElement el =  chr.getNearestElement(1);
				if(!(el instanceof Snp)){//TODO : replace with exception
					System.out.println("first in genome not filled with snp.");
					writer.close();
				return;
				}
				int genotypeLength = ((Snp) el).getGenotypes().length;
				String nextLine = "";
				for (int i=0; i < genotypeLength ; ++i){
				nextLine = 	i + " " + "NA" + i + " " +"0 0 1 -9";
				writer.println(nextLine);
			}
		}
		writer.close();		
	}

	
	// ============================================================================
	// PRIVATE METHODS
		

	// ============================================================================
	// GETTERS AND SETTERS

	public int getNumElements() { return numElements_; }
	public int getNumChromosomes() { return chromosomes_.size(); }
	
	public HashMap<String, Chromosome> getChromosomes() { return chromosomes_; }
	public Chromosome getChromosome(String chrId) { return chromosomes_.get(chrId); }
	
}
