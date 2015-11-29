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

import java.util.ArrayList;
import java.util.Iterator;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.main.Pascal;

public class GeneWithItsSnps {

	private Gene gene_;
	private ArrayList<Snp> snpList_;

	public GeneWithItsSnps(Gene gene,ArrayList<Snp> snpList){
		gene_=gene;
		snpList_=snpList;
		checkAndExecuteCodingRemoval();
	}
	
	public ArrayList<Snp> getSnpList(){
		return snpList_;
		
	}
	public Gene getGene(){
		return gene_;
	}
	public String getGeneId(){
		return gene_.getId();
	}
		
	private void checkAndExecuteCodingRemoval(){
		if(Pascal.set.removeCodingSnpsOfOtherGenes_){
			removeCodingSnpsOfOtherGenes();
		}
	}
	
	public int getNrOfSnps(){
		return(snpList_.size());
	}
	/** Remove coding snps of other genes in that window */
	private void removeCodingSnpsOfOtherGenes() {
		
		Iterator<Snp> iter = snpList_.iterator();
		while (iter.hasNext()) {
			Snp snp = iter.next();
			if (snp.isCodingForOtherGene(gene_.getId()))
				iter.remove();
		}
	}
	
	/** Check if gene has greater than zero and smaller than maxSnpsPerGene SNPs */
	public boolean checkNrOfSnps() {
		if (snpList_.size() == 0 || (Pascal.set.maxSnpsPerGene_ > 0 && snpList_.size() > Pascal.set.maxSnpsPerGene_)) {
			if (Pascal.set.verbose_) {
				Pascal.print("\t" + ConvenienceMethods.padRight("0h 0min 0s 0ms", 22));
				if (snpList_.size() == 0)
					Pascal.print("Gene has no SNPs\n");
				else
					Pascal.print("Gene exceeds max number of SNPs defined in settings file\n");			
			
			}
			return false;
		}
		return true;
	}		
	
	
	/** Print the gene name and number of SNPs to the console */
	public void printGeneNameAndNrOfSnps(){
		String symb = ConvenienceMethods.padRight(((gene_.symbol_ == null) ? "NA" : gene_.symbol_), 16);
		Pascal.print(ConvenienceMethods.padRight(gene_.id_, 18) + symb);
		Pascal.print(String.format("%6d", snpList_.size()));
	}
}

