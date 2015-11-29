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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Pascal;


/**
 * Compute enrichment scores for the given gene sets
 */
public class GeneSet {

	/** The name of this gene set */
	private String id_ = null;
	
	/** The ordered set of genes (some may be meta-genes) */
	 HashSet<Gene> genes_ = null;
	/** The set of meta-genes (this is a subset of genes_) */
	private HashSet<MetaGene> metaGenes_ = null;
	/** The p-value computed using a simulation test */
	private double simulPvalue_ = Double.NaN;
	/** The p-value computed using a simulation test Weighted */
	private double simulPvalueWeighted_ = Double.NaN;
	/** The p-value computed using chi-squared distribution */
	private double chi2Pvalue_ = Double.NaN;
	/** The p-value computed using rank-sum test*/
	private double rankSumPvalue_ = Double.NaN;
	/** The p-value computed using the hypergeometric distribution */
	private double[] hypGeomPvalues_ = null;
	/** The p-value computed using the gamma distribution */
	private double[] gammaPvalues_ = null;
	/** Indicates if there is depletion of this gene set */
	private boolean depletion_ = false;

	private double[] expHypPvalues_;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneSet(String id) {
		
		id_ = id;
		genes_ = new HashSet<Gene>();
	}

	
	/** Constructor */
	public GeneSet(String id, Collection<Gene> genes) {
		
		this(id);
		genes_.addAll(genes);
	}

	
	// ----------------------------------------------------------------------------

	/** Add the given gene to the set */
	public void add(Gene gene) {
		
		genes_.add(gene);
	}
	

	// ----------------------------------------------------------------------------

	/** Create meta-genes for neighboring genes closer than Settings.mergeGenesDistance_ */
	public void createMetaGenes(HashMap<String, MetaGene> allMetaGenes) {
		
		setMetaGenes(new HashSet<MetaGene>());
		final int mergeDist = (int) Math.floor(1000000 * Pascal.set.mergeGenesDistance_);
		
		if (genes_.size() < 1)
			return;
		
		// Sort genes by position
		ArrayList<Gene> genesSorted = Gene.sortByPosition(genes_);
		
		// Iterate over genes
		Iterator<Gene> iter = genesSorted.iterator();
		// The first gene (above we made sure there is at least one element)
		Gene prevGene = iter.next();
		// The next meta-gene
		TreeSet<Gene> merge = null; 
		
		while (iter.hasNext()) {
			Gene curGene = iter.next();
			
			// If this gene's window overlaps with the previous gene
			if (curGene.start_ - mergeDist <= prevGene.end_ && curGene.chr_.equals(prevGene.chr_)) {
				if (merge == null) {
					// Create new meta-gene from the two genes 
					merge = new TreeSet<Gene>();
					merge.add(prevGene);
					merge.add(curGene);
				} else {
					// Add to growing meta-gene
					merge.add(curGene);
				}
			
			} else {
				// No overlap => finish previous meta-gene
				if (merge != null) {
					MetaGene metaGene = new MetaGene(merge);
					String metaGeneId = metaGene.getId();
					// if metagene has been seen before add it add it from list,else add new.
					if(allMetaGenes.containsKey(metaGeneId)){						
						getMetaGenes().add(allMetaGenes.get(metaGeneId));
					}
					else{
						getMetaGenes().add(metaGene);
						}
					
					// Reset genes to be merged
					merge = null;
				}
			}
			prevGene = curGene;
		}

		// If the last genes were merged, finish up
		if (merge != null) {
			MetaGene metaGene = new MetaGene(merge);
			String metaGeneId = metaGene.getId();
			// if metagene has been seen before add it add it from list,else add new.
			if(allMetaGenes.containsKey(metaGeneId)){
				getMetaGenes().add(allMetaGenes.get(metaGeneId));
			}
			else{	
				getMetaGenes().add(metaGene);
			}
		}
			
		// Add the meta-genes and remove the corresponding individual genes (has to be done here to not screw up the iterator)
		for (MetaGene meta : getMetaGenes()) {
			genes_.add(meta);
			for (Gene g : meta.getGenes())
				genes_.remove(g);
		}
	}
	
	/** removes metaGenes that didn't have a score. has to be removed 
	 * from metaGenes_ and genes_.
	 * */
	public HashSet<MetaGene> removeMetaGenesWithoutScore(){
		
		HashSet<MetaGene> metaGenesToBeRemoved = new HashSet<MetaGene>();
		if(getMetaGenes()==null){
			return metaGenesToBeRemoved;
		}
		for (MetaGene currentGene : getMetaGenes()){
			
			if (currentGene.getChi2Stat()==-1){
				metaGenesToBeRemoved.add(currentGene);
			}			
		}
		for (MetaGene geneToBeRemoved : metaGenesToBeRemoved){
			genes_.remove(geneToBeRemoved);
			getMetaGenes().remove(geneToBeRemoved);			
		}
		return(metaGenesToBeRemoved);
		
		
	}
	
	// ----------------------------------------------------------------------------

	public void countMergedGenes(int[] minMaxTot) {

		for (MetaGene meta : getMetaGenes()) {
			int numMerged = meta.getGenes().size();
			minMaxTot[0] = Math.min(minMaxTot[0], numMerged);
			minMaxTot[1] = Math.max(minMaxTot[1], numMerged);
			minMaxTot[2] += numMerged;
		}
	}

	
	// ----------------------------------------------------------------------------

	/** String representation of this set (tab-separated list of gene IDs) */
	@Override
	public String toString() {
		
		if (genes_.size() == 0)
			return "empty_set";
		
		// Sort genes by position
		ArrayList<Gene> genesSorted = Gene.sortByPosition(genes_);

		Iterator<Gene> iter = genesSorted.iterator();
		String str = iter.next().id_;
		while (iter.hasNext())
			str += "\t" + iter.next().id_;
		
		return str;
	}

	

	// ----------------------------------------------------------------------------

	/** Write enrichment to text file */
	public String getResultsAsString() {

		// Id
		String str = id_;

		// Enrichment pvals
		if (Pascal.set.useChi2_) {
			str += "\t" + Pascal.utils.toStringScientific10(chi2Pvalue_);
		}
		if (Pascal.set.useSimulation_)
			str += "\t" + Pascal.utils.toStringScientific10(simulPvalue_);
		
		if (Pascal.set.useSimulationWeightedSampling_)
			str += "\t" + Pascal.utils.toStringScientific10(simulPvalueWeighted_);
		
		
		if (Pascal.set.useRankSum_)
			str += "\t" + Pascal.utils.toStringScientific10(rankSumPvalue_);
		
		if (Pascal.set.useHypGeom_)
			for (int i=0;i<hypGeomPvalues_.length; i++)
				str += "\t" + Pascal.utils.toStringScientific10(hypGeomPvalues_[i]);
		
		if (Pascal.set.useGamma_)
			for (int i=0;i<gammaPvalues_.length; i++)
				str += "\t" + Pascal.utils.toStringScientific10(gammaPvalues_[i]);
		
		if (Pascal.set.useExpHyp_)
			for (int i=0;i<expHypPvalues_.length; i++)
				str += "\t" + Pascal.utils.toStringScientific10(expHypPvalues_[i]);
		// Sort genes by score
		
		// Comma-separated list of gene ids
		//TODO: this option is brocken atm
		if (Pascal.set.writeDetailedOutput_) {
			ArrayList<Gene> genesSorted = new GeneScoreList(genes_, true).getGenes();
			
		str += "\t";
		if (genesSorted.size() > 0) {
			Iterator<Gene> iter = genesSorted.iterator();
			str += iter.next().id_;
			while (iter.hasNext())
				str += "," + iter.next().id_;
		}
		
		// Comma-separated list of gene scores (kind of redundant with the gene score files, but could be convenient)
		
			str += "\t";
			if (genesSorted.size() > 0) {
				Iterator<Gene> iter = genesSorted.iterator();
				str += Pascal.utils.toStringScientific10(iter.next().getScore(0));
				while (iter.hasNext())
					str += "," + Pascal.utils.toStringScientific10(iter.next().getScore(0));
			}
		}
		return str;
	}

	
	// ----------------------------------------------------------------------------

	/** Get the header corresponding to getResultsAsString() */
	public static String getResultsAsStringHeader() {
		
		String str = "Name";
		if (Pascal.set.useChi2_)
			str += "\t" + "chi2Pvalue";
		if (Pascal.set.useSimulation_)
			str += "\t" + "empPvalue";
		if (Pascal.set.useSimulationWeightedSampling_)
			str += "\t" + "simulWeightedSamplingPvalue";
		if (Pascal.set.useRankSum_)
			str += "\t" + "ranksumPvalue";
		if(Pascal.set.useHypGeom_){
			for (int i=0;i< Pascal.set.hypGeomQuantiles_.size(); i++)
				str += "\t" + "hypgeomPvalue" + Pascal.set.hypGeomQuantiles_.get(i);  		
		}
		if(Pascal.set.useGamma_){
			for (int i=0;i< Pascal.set.gammaShapeParameters_.size(); i++)
				str += "\t" + "gammaPvalue" + Pascal.set.gammaShapeParameters_.get(i);
		}
		if(Pascal.set.useExpHyp_){
			for (int i=0;i< Pascal.set.expHypParameters_.size(); i++)
				str += "\t" + "expHypPvalue" + Pascal.set.expHypParameters_.get(i);
		}
		
		if (Pascal.set.writeDetailedOutput_){
			str += "\tGenes";
			str += "\tGeneScores";
		}
		return str;
	}
	
	public ArrayList<String> getGeneIds(){
		ArrayList<String> out = new ArrayList<String>(genes_.size());
		for(Gene g : genes_){
			out.add(g.id_);
		}	
		return(out);
	}
	
	
	// ============================================================================
	// GETTERS AND SETTERS
		
	public HashSet<Gene> getGenes() { return genes_; }
	public HashSet<MetaGene> getMetaGenes() { return metaGenes_; }
	public String getId() { return id_; }
	

	public void setSimulPvalue(double x) { simulPvalue_ = x; }
	public double getSimulPvalue() { return simulPvalue_; }
	
	public void setSimulPvalueWeighted(double x) { simulPvalueWeighted_ = x; }
	public double getSimulPvalueWeighted() { return simulPvalueWeighted_; }
	
	public void setHypGeomPvalues(double[] x) { hypGeomPvalues_ = x; }
	public double[] getHypGeomPvalues() { return hypGeomPvalues_; }
	
	public void setGammaPvalues(double[] x) { gammaPvalues_ = x; }
	public double[] getGammaPvalues() { return gammaPvalues_; }
	
	public void setExpHypPvalues(double[] x) { expHypPvalues_ = x; }
	public double[] getExpHypPvalues() { return expHypPvalues_; }
	
	public void setChi2Pvalue(double x) { chi2Pvalue_ = x; }
	public double getChi2Pvalue() { return chi2Pvalue_; }

	public void setRankSumPvalue(double x) { rankSumPvalue_ = x; }
	public double getRankSumPvalue() { return rankSumPvalue_; }

	public void setDepletion(boolean b) { depletion_ = b; }
	public boolean getDepletion() { return depletion_; }


	public void setMetaGenes(HashSet<MetaGene> metaGenes_) {
		this.metaGenes_ = metaGenes_;
	}


	
	
}
