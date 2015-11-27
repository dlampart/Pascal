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

import java.io.IOException;
import java.util.LinkedList;

import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Pascal;
import ch.unil.genescore.main.Settings;

public class LoadedGenotypes {
	
	/** A stack with the SNPs for which genotypes are currently loaded, newly loaded snps are put on top */
	private LinkedList<Snp> snpsWithGenotypes_ = null;

	/** The current chromosome */
	private String loadedChr_ = "not_initialized";	
	/** The start and end coordinates of the region for which genotypes are currently loaded */
	private int loadedStart_ = -1;
	/** The start and end coordinates of the region for which genotypes are currently loaded */
	private int loadedEnd_ = -1;
	/** Flag set when end of file is reached */
	protected boolean endOfFile_ = false;
	private ReferencePopulation refpop_ =null;
	public LoadedGenotypes(){
		snpsWithGenotypes_ = new LinkedList<Snp>();
	}
	
	public LoadedGenotypes(ReferencePopulation refpop){
		snpsWithGenotypes_ = new LinkedList<Snp>();
		refpop_ = refpop;
	}
	
	

	/** 
	 * Load genotypes for SNPs in the window around the given gene, delete genotypes of SNPs
	 * that are not needed anymore. 
	 */
	public void updateLoadedGenotypes(GenomicElement gene) {

		// If the gene is not on the current chromosome, delete all loaded genotypes, reinitialize and open binary file
		if (!gene.chr_.equals(loadedChr_))
			refpop_.initialize(gene.chr_);
		else if (endOfFile_)
			return;
			
		// Boundaries of region for which genotypes should be loaded (note, with meta-genes we don't really know what's up and down, that's why we use max())
		int d = Math.max(Settings.geneWindowDownstream_, Settings.geneWindowUpstream_);
		int laxityFactor=3000000;
		int newStart = gene.start_ - d - laxityFactor;
		int newEnd = gene.end_ + d + laxityFactor;;		
			
		// Remove snps before newStart from the front of the list and delete their genotypes
		while (getSnpsWithGenotypes().size() > 0 && getSnpsWithGenotypes().getFirst().start_ < newStart)
			getSnpsWithGenotypes().poll().setGenotypes(null);

		loadedStart_ = newStart;		
		while (loadedEnd_ < newEnd && !endOfFile_)
			refpop_.loadNextGenotype(); // updates loadedEnd_

		// The first snp of the list should be the first gwas snp > loadedStart
		assert getSnpsWithGenotypes().size() == 0 || getSnpsWithGenotypes().getFirst().start_ >= loadedStart_;
		// The last snp of the list should be the first gwas snp > loadedEnd (we load one too many), unless the end of the file has been reached
		assert getSnpsWithGenotypes().size() == 0 || endOfFile_ || getSnpsWithGenotypes().getLast().start_ >= loadedEnd_;
	}
	


	
	public LinkedList<Snp> getSnpsWithGenotypes() {
		return snpsWithGenotypes_;
	}

	public void setSnpsWithGenotypes(LinkedList<Snp> snpsWithGenotypes) {
		this.snpsWithGenotypes_ = snpsWithGenotypes;
	}

}
