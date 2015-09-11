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
import java.util.Collection;
import java.util.Iterator;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.main.Utils;

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
		if(Settings.removeCodingSnpsOfOtherGenes_){
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
	
	public boolean checkNrOfSnps(GeneResultsSnpsOutOfBounds GeneResultsZeroOrAboveLimit_){
		if (getNrOfSnps() == 0 || (Settings.maxSnpsPerGene_ > 0 && getNrOfSnps() > Settings.maxSnpsPerGene_)) {
		
			GeneResultsZeroOrAboveLimit_.addToMap(this);
					
			if (Settings.verbose_) {
				Main.print("\t" + Utils.padRight("0h 0min 0s 0ms", 22));
				if (getNrOfSnps() == 0)
					Main.print("Gene has no SNPs\n");
				else
					Main.print("Gene exceeds max number of SNPs defined in settings file\n");			
			
			}
			return false;
		}
		return true;
	}		
	
	public void printGeneNameAndNrOfSnps(){
		String symb = Utils.padRight(((gene_.symbol_ == null) ? "NA" : gene_.symbol_), 16);
		Main.print(Utils.padRight(gene_.id_, 18) + symb);
		Main.print(String.format("%6d", snpList_.size()));
	}
}

