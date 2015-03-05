/*******************************************************************************
 * Copyright (c) 2015 IBM Corporation and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     IBM Corporation - initial API and implementation
 *******************************************************************************/
package ch.unil.genescore.pathway;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.Genome;
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
