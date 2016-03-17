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
package ch.unil.genescore.pathway;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.Genome;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.Snp;


/**
 * 
 */
public class MetaGene extends Gene {

	/** The genes that are merged in this meta-gene */
	private TreeSet<Gene> genes_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public MetaGene(TreeSet<Gene> genes) {

		super(null);
		initialize(genes);
	}

	
	// ----------------------------------------------------------------------------

	/** Get the snps that are in the gene windows merged in this meta-gene */
	@Override
	public ArrayList<Snp> findSnps(Genome snps) {
		
		// Use a set to not include the same snp multiple times
		TreeSet<Snp> geneSnps = new TreeSet<Snp>();
		for (Gene gene : genes_)
			geneSnps.addAll(gene.findSnps(snps));
		
		return new ArrayList<Snp>(geneSnps);
	}

	
	// ============================================================================
	// PRIVATE METHODS

	/** Initialize this meta-gene with the given genes */
	private void initialize(TreeSet<Gene> genes) {
		
		if (genes == null || genes.size() < 2)
			throw new RuntimeException("Cannot create meta-gene with less than two genes");
		
		genes_ = genes;
		
		Iterator<Gene> iter = genes.iterator();
		Gene first = iter.next();
		
		id_ = "meta_" + first.id_;
		symbol_ = "meta_" + first.symbol_;
		chr_ = first.chr_;
		start_ = first.start_;
		end_ = first.end_;
		// Note, posStrand_ is not defined for meta-genes
		
		while (iter.hasNext()) {
			Gene gene = iter.next();
			
			// Check ordering
			if (!gene.chr_.equals(chr_))
				throw new RuntimeException("Genes of a meta-gene must be on same chromosome");
			if (gene.start_ < start_)
				throw new RuntimeException("Genes of a meta-gene must be ordered by start position");
			
			// Update end
			if (gene.end_ > end_)
				end_ = gene.end_;
			
			// Update id and symb
			id_ += "_" + gene.id_;
			symbol_ += "_" + gene.symbol_;
		
		}
	}
	
	/** Remove Snps in gene body region */
	public void removeSnpsInGeneBody(ArrayList<Snp> snpList) {
		Iterator<Snp> it = snpList.iterator();
		while (it.hasNext()){
			Snp snp = it.next();		
			for(Gene gene : genes_){
				if(snp.chr_.equals(gene.chr_) &&  snp.start_ >= gene.start_ &&  snp.end_ <= gene.end_){
					it.remove();
					break;
				}	
			}				
		}
	}
	
//	public void removeHighMafSnps(Collection<Snp> snpList){
//		if(Settings.useMaxMafCutoff_ < 0.5){
//			Iterator<Snp> it = snpList.iterator();
//			while (it.hasNext()){
//				Snp snp = it.next();
//				if (snp.getMaf() > Settings.useMaxMafCutoff_ || snp.getAlleleSd()==0)
//					it.remove();
//			}
//		}
//	}
	
	// ============================================================================
	// GETTERS AND SETTERS
		
	public TreeSet<Gene> getGenes() { return genes_; }
	
	@Override
	public ArrayList<String> getSymbolList(){
		ArrayList<String> myList = new ArrayList<String>();
		for (Gene gene : genes_){
			myList.add(gene.symbol_);
		}
			return(myList);
	}
}
