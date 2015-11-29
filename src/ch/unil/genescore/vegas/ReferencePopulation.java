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

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.Genome;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Pascal;


/**
 * 
 */
public class ReferencePopulation {

	/** The binary file with the reference population */
	protected ArrayList<String> chromosomes_ = null;
	
	/** The complete set of SNPs for which we have p-values */
	protected HashMap<String, Snp> gwasSnps_ = null;
	/** A stack with the SNPs for which genotypes are currently loaded, newly loaded snps are put on top */
	
	protected LoadedGenotypes loadedGenotypes_ = null;
	/** The set of SNPs that are weighed*/
	protected HashMap<String, Snp> relevantSnps_ = null;
		
	/** A genome instance populated with the gwas-snps */
	private Genome snps_ = null;

	/** The current chromosome */
	private String loadedChr_ = "not_initialized";	
	/** The start and end coordinates of the region for which genotypes are currently loaded */
	private int loadedStart_ = -1;
	/** The start and end coordinates of the region for which genotypes are currently loaded */
	private int loadedEnd_ = -1;
	/** Flag set when end of file is reached */
	private boolean endOfFile_ = false;
	
	/** The path and prefix of the ref pop files */
	//protected String genotypeFileHandler_.getFilePrefix() = null;
	/** The binary file being read */
	private DataInputStream inStream_ = null;
	
	private GenotypeFileHandler genotypeFileHandler_ = null;			
	
	
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public ReferencePopulation() {
		
		setGenotypeFileHandler(new GenotypeFileHandler());		
		//genotypeFileHandler_.getFilePrefix() = Settings.refPopDirectory_ + "/" + Settings.refPopFilePrefix_ + ".";
		
		// Initialize list of chromosomes to be considered
		chromosomes_ = new ArrayList<String>();
		if (Pascal.set.chromosome_ == null || Pascal.set.chromosome_.length() == 0)
			for (int i=1; i<=22; i++)
				chromosomes_.add("chr" + i);
		else
			chromosomes_.add(Pascal.set.chromosome_);
		
		// Create missing binary files if text or tped files were specified in settings
		for (String chr : chromosomes_)
			getGenotypeFileHandler().writeBinaryFiles(chr);
		
		
		// Initialize empty set in case no file with relevant snps is loaded
		relevantSnps_ = new HashMap<String, Snp>();
	}
	
	// ----------------------------------------------------------------------------

	/** load snp positions from the binary ref pop files */
	public void initializeSnps() { 
				
		for (String chr : chromosomes_)
			loadGwasAndRelevantSnpsPos(chr);
	}


	// ----------------------------------------------------------------------------

	/** 
	 * Load genotypes for SNPs in the window around the given gene, delete genotypes of SNPs
	 * that are not needed anymore. 
	 */
	public void updateLoadedGenotypes(GenomicElement gene) {

		// If the gene is not on the current chromosome, delete all loaded genotypes, reinitialize and open binary file
		if (!gene.chr_.equals(loadedChr_))
			initialize(gene.chr_);
		else if (isEndOfFile())
			return;
			
		// Boundaries of region for which genotypes should be loaded (note, with meta-genes we don't really know what's up and down, that's why we use max())
		int d = Math.max(Pascal.set.geneWindowDownstream_, Pascal.set.geneWindowUpstream_);
		int laxityFactor=3000000;
		int newStart = gene.start_ - d - laxityFactor;
		int newEnd = gene.end_ + d + laxityFactor;;		
			
		// Remove snps before newStart from the front of the list and delete their genotypes
		while (loadedGenotypes_.getSnpsWithGenotypes().size() > 0 && loadedGenotypes_.getSnpsWithGenotypes().getFirst().start_ < newStart)
			loadedGenotypes_.getSnpsWithGenotypes().poll().setGenotypes(null);

		loadedStart_ = newStart;		
		while (loadedEnd_ < newEnd && !isEndOfFile())
			loadNextGenotype(); // updates loadedEnd_

		// The first snp of the list should be the first gwas snp > loadedStart
		assert loadedGenotypes_.getSnpsWithGenotypes().size() == 0 || loadedGenotypes_.getSnpsWithGenotypes().getFirst().start_ >= loadedStart_;
		// The last snp of the list should be the first gwas snp > loadedEnd (we load one too many), unless the end of the file has been reached
		assert loadedGenotypes_.getSnpsWithGenotypes().size() == 0 || isEndOfFile() || loadedGenotypes_.getSnpsWithGenotypes().getLast().start_ >= loadedEnd_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Delete loaded genotypes and re-initialize with the given chromosome */
	public void initialize(String chr) {

		// Delete loaded genotypes
		if (loadedGenotypes_ !=null && loadedGenotypes_.getSnpsWithGenotypes() != null) {
			for (Snp snp : loadedGenotypes_.getSnpsWithGenotypes())
				snp.setGenotypes(null);
		}

		// Close the previous file (if there is one open)
		closeDataInputStream();

		// Reset the treemap and loaded region	
		loadedGenotypes_ = new LoadedGenotypes();
		loadedChr_ = chr;
		loadedStart_ = -1;
		loadedEnd_ = -1;
		
		// Open input binary file
		String filename = getGenotypeFileHandler().getFilePrefix() + chr + ".gnt.ser.gz";
		openDataInputStream(filename);
		
		// Read genotype vector length and phased flag
		try {
			int length = getInStream().readInt();
			boolean phased = getInStream().readBoolean();
			
			Snp.setGenotypeLength(length);
			Snp.setGenotypeIsPhased(phased);
			
		} catch (IOException e) {
			Pascal.error(e, "Error reading genotype length");
		}
	}


	// ----------------------------------------------------------------------------

	/** Load the next snp from the input stream */
	public void loadNextGenotype() {
		
		try {
			String snpId = getInStream().readUTF();
			
			// End of file is marked with the version ID
			if (snpId.equals(getGenotypeFileHandler().getBinaryFileVersionID())) {
				setEndOfFile(true);
				return;
			}
			
			Snp snp = gwasSnps_.get(snpId);
			if (snp == null)
				snp = relevantSnps_.get(snpId);	
			
			// If this is not a gwas snp, skip
			if (snp == null || snp.start_ < loadedStart_) {
				Snp.skip(getInStream());
				return;
			}			
			// Load the genotype of this snp
			snp.readGenotype(getInStream());			
			// Update
			loadedGenotypes_.getSnpsWithGenotypes().add(snp);
			loadedEnd_ = snp.start_;
		
		} catch (Exception e) {
			Pascal.error(e, "Error loading genotypes");
		}

	}

	
	// ----------------------------------------------------------------------------

	/** Close the currently open binary filem, reset endOfFile_ */
	protected void closeDataInputStream() {
		
		try {
			setEndOfFile(false);
			if (getInStream() != null)
				getInStream().close();
		} catch (IOException e) {
			Pascal.error(e, "Error closing binary input genotype file");
		}
	}
			
	
	/**get inputStream-Object that that iterates over snp positions*/
	 public SnpPositionStream getSnpPositionStream(String chr){		 
		 return(getGenotypeFileHandler().getSnpPositionStream(chr));
			

	 }
	 /**get inputStream-Object that that iterates over snp positions*/
	 public WrappedSnpPositionStream getSnpPositionStream(String[] chrs){		 		
		 return(getGenotypeFileHandler().getSnpPositionStream(chrs));
	 }
	 
	 
	// ----------------------------------------------------------------------------
	 /** Load the genomic coordinates for gwasSnps_ of the given chromosome */
		private void loadGwasAndRelevantSnpsPos(String chr) {
			
			int counter=0;
			SnpSerializedPositionStream str = getGenotypeFileHandler().getSnpSerializedPositionStream(chr);
			while (str.lineAvailable()) {			
				Snp newSnp = str.readLine();
				String snpId = newSnp.getId();
				Snp snp = null;
				try{
					if ((gwasSnps_.containsKey(snpId) || relevantSnps_.containsKey(snpId))){						
						if (gwasSnps_.containsKey(snpId))
							snp = gwasSnps_.get(snpId);
						else 
							snp = relevantSnps_.get(snpId);
							
						if (snp.chr_ != "none" || snp.start_ != -1 || snp.end_ != -1){
							throw new RuntimeException("snp seems to have been set before");
						}						
						snp.copyPosAndAllele(newSnp);
					}
				}
				catch(DataInconsistencyException e){					
					counter++;
					System.err.println("Caught DataInconsistencyException: " + e.getMessage() + counter);
					continue;			            
				}	
		}
	}
	 
	 
	public void addToRelevantSnps(Set<String> relevantSnpIds){
		
		for (String relevantSnpId : relevantSnpIds)
			relevantSnps_.put(relevantSnpId, new Snp(relevantSnpId));	
		System.out.println(relevantSnps_.size());
		
	}
	
	
	// ----------------------------------------------------------------------------

	/** Load the snp pvals and positions */
	public void loadGwasAndRelevantSnps(){
		// Initialize reference population instance
			
		// Load pvals (all chromosomes)
		GwasSnps snps = new GwasSnps();
		if (!Pascal.set.withZScore_){
			snps.setHeader(false);
			snps.loadSnpPvals(Pascal.set.snpPvalFile_);
		}
		else{
			snps.setHeader(true);
			snps.loadSnpPvalZval(Pascal.set.snpPvalFile_);
		}
		// Remove SNPs not in snpPruneFile
		File snpFilterFile = Pascal.set.snpFilterFile_;
		if (snpFilterFile != null)
			snps.keepOnlySnpsInFilterListFile(snpFilterFile);

//TODO: process again		if (Settings.removeCodingSnpsOfOtherGenes_)
//			snps.loadCodingSnps(Settings.codingSnpsFile_);
		// Load snp positions from ref pop files for chromosome specified in Settings.chromosome				
		
		setGwasSnps(snps.getSnpsInList());
		
		
		initializeSnps();	
		if (Pascal.set.useOnlyGwasSnps_){
			 useOnlyGwasSnps();
		}
        //TODO: this has to be somehow taken care of Don't freak
//		if (Settings.populationFormat_=="tped")
//           loadReferenceBaseCoding(refPopBasesFile,snps);

		// Load coding snps
		Genome gSnps = new Genome();
		gSnps.addElements(snps.getSnpsInList().values());		
		
		// Print location of SNPs of this study as bed file
		String filename = Pascal.set.gwasName_ + ".snps" + Pascal.set.chromFileExtension_;
		if (Pascal.set.writeSnpBedFile_) {
			File bedFile = new File(Pascal.set.outputDirectory_, filename + ".txt");
			gSnps.writeBedFile(bedFile);
		}
		if (Pascal.set.writeTpedFile_){
			File tfamFile = new File(Pascal.set.outputDirectory_, filename + ".tfam");
            gSnps.writePseudoTfamPlinkFile(tfamFile);
            File tpedFile = new File(Pascal.set.outputDirectory_, filename + ".tped");
            gSnps.writeTpedPlinkFile(tpedFile);			
		}
		setGenomeSnps(gSnps);
		
		// Add snps to genome, filters for specified chromosome
	}
	

	
	// ----------------------------------------------------------------------------

	/** Open a data input stream */
	protected void openDataInputStream(String filename)  {
		
		if (Pascal.set.verbose_)
			Pascal.println("Reading file: " + filename);

		try {
			FileInputStream infile = new FileInputStream(filename);
			GZIPInputStream gzip = new GZIPInputStream(infile);
			BufferedInputStream buf = new BufferedInputStream(gzip);
			DataInputStream inStream = new DataInputStream(buf);
			
			// Check that the file starts with the right versionId
			String versionId = inStream.readUTF();
			if (!versionId.equals(getGenotypeFileHandler().getBinaryFileVersionID()))
				throw new RuntimeException("Incompatible version ID of binary file, delete the file to create a new one");
			
			inStream_=inStream;
			
		} catch (Exception e) {
			throw new RuntimeException("Could not open binary input file: " + filename);
		}
	}
	
	public ArrayList<Snp> findSnps(Gene el){
		return(el.findSnps(snps_));
	}
	public HashSet<Snp> findSnpsHash(Gene el){
		return(el.findSnpsHash(snps_));
	}
	
	public void useOnlyGwasSnps(){
		relevantSnps_ = gwasSnps_;
	}
	
	// ============================================================================
	// GETTERS AND SETTERS
	public void  setGenomeSnps(Genome snps){snps_=snps;}
	public void  setGwasSnps(HashMap<String, Snp> gwasSnps){gwasSnps_=gwasSnps;}
	public void  setRelevantSnps(HashMap<String, Snp> relevantSnps){relevantSnps_=relevantSnps;}
	public LinkedList<Snp> getSnpsWithGenotypes(){return loadedGenotypes_.getSnpsWithGenotypes();}

	protected GenotypeFileHandler getGenotypeFileHandler() {
		return genotypeFileHandler_;
	}

	protected void setGenotypeFileHandler(GenotypeFileHandler genotypeFileHandler_) {
		this.genotypeFileHandler_ = genotypeFileHandler_;
	}

	protected DataInputStream getInStream() {
		return inStream_;
	}

	protected void setInStream(DataInputStream inStream_) {
		this.inStream_ = inStream_;
	}

	protected boolean isEndOfFile() {
		return endOfFile_;
	}

	protected void setEndOfFile(boolean endOfFile_) {
		this.endOfFile_ = endOfFile_;
	}
}
