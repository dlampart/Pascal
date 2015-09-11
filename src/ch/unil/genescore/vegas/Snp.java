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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.HashSet;


import org.apache.commons.math3.distribution.TDistribution;

import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.DistributionMethods;

/**
 * A SNP with its p-value.
 */
public class Snp extends GenomicElement {

	/** Length of genotype vector (written at beginning of binary file) */
	static private int genotypeLength_ = -1;
	/** True if genotype is phased (written at beginning of binary file) */
	static private boolean genotypeIsPhased_ = true;

	/** True if genotype is phased (written at beginning of binary file) */
	private boolean genotypeIsPhased2_ = true;
	
	/** Number of bytes of one snp in the binary genotype file (not including the ID) */
	static private int numBytes_ = -1; 

	
	/** The real p-value, optionally followed by p-values from phenotype permutations */
	private double pval_ = -1;
	/** The pvals transformed to 1 df chi squared statistic (VEGAS) */
	private double chi2Stat_ = Double.NaN;
	/** The z-score */
	private double zscore_ = -1;	
	/** The genotypes coded as 0, 1 */
	private byte[] genotypes_ = null;
	/** Minor allele frequency */
	private double maf_ = -1;
	/** Allele mean */
	private double alleleMean_ = -1;
	/** Allele standard deviation */
	private double alleleSd_ = -1;
	
	private char majorAllele_ = 'N';
	private char minorAllele_ = 'N';
	
	/** Number of alleles that are not NA in the genotype (not written to binary file) */
	private transient int alleleNonMissing_ = -1;
	
	/** Genes for which this snp is a coding mutation */
	private HashSet<String> coding_ = null;
	
	private static HashMap<Byte, String> fromDephasedToPlinkMap_ = new HashMap<Byte, String>();
	private static HashMap<Byte, String> fromPhasedToPlinkMap_ = new HashMap<Byte, String>();
	private static HashMap<String, Byte> fromPlinkToPhasedMap_ = new HashMap<String, Byte>(2);
	private static HashMap<String, Byte> fromPlinkToDephasedMap_ = new HashMap<String, Byte>(3);
		
		
	// ============================================================================
	// PUBLIC METHODS
	
	/** Initialize plink maps */
	static {
		fromDephasedToPlinkMap_.put(new Byte((byte) 0), "1 1");
		fromDephasedToPlinkMap_.put(new Byte((byte) 1), "1 2");
		fromDephasedToPlinkMap_.put(new Byte((byte) 2), "1 2");
		fromDephasedToPlinkMap_.put(new Byte((byte) -9), "0 0");	
		
		fromPhasedToPlinkMap_.put(new Byte((byte) 0), "1");
		fromPhasedToPlinkMap_.put(new Byte((byte) 1), "2");
		fromPhasedToPlinkMap_.put(new Byte((byte) -9), "0");
				
		fromPlinkToPhasedMap_.put("1", new Byte((byte) 0));
		fromPlinkToPhasedMap_.put("2", new Byte((byte) 1));
		fromPlinkToPhasedMap_.put("0", new Byte((byte) -9));
		
		fromPlinkToDephasedMap_.put("1 1", new Byte((byte) 0));
		fromPlinkToDephasedMap_.put("1 2", new Byte((byte) 1));
		fromPlinkToDephasedMap_.put("2 2", new Byte((byte) 2));
	}
	
	
	/** Constructor */
	public Snp(String id) {
		super(id);
	}

	/** Constructor */
	public Snp(String id, double pval) {
		super(id);
		pval_ = pval;
		computeChiSquaredStatistics();		
	}
	
	/** Constructor */
	public Snp(String id, double pval, double zscore) {
		super(id);		
		pval_ = pval;
		zscore_ = zscore;	
	}


	/** calculate new p-values and zscores via linear regression equation. */
	public void recalcPvalsAndZscores(double[] phenotype, double phenotypeMean, double phenotypeSd, TDistribution myT) {
		
		if (phenotype.length!=genotypes_.length)
			throw new RuntimeException("phenotype vector not same length as genotype vector.");
		//calculate allele mean
		double alleleMean=0;
		for (int i=0; i<genotypes_.length; i++){
			alleleMean=alleleMean+genotypes_[i];
		}
		alleleMean=alleleMean/genotypes_.length;
				
		//all formulas are from Regression Script by Kuensch (ETHZ)
		//calculating alpha and beta of regression equation.
		double betaUpper=0;
		double betaLower=0;
		for (int i=0; i<genotypes_.length; i++){
			betaUpper=betaUpper + phenotype[i]*(genotypes_[i]-alleleMean_);
			betaLower=betaLower + Math.pow((genotypes_[i]-alleleMean_), 2);
		
		}
		double beta=betaUpper/betaLower;
		double alpha=phenotypeMean + beta*alleleMean_;
		//calculating stdErr
		double stdErr=0;
		for (int i=0; i<genotypes_.length; i++){
			stdErr=stdErr + Math.pow((phenotype[i]-alpha+beta*genotypes_[i]),2);					
		}
		stdErr=stdErr/(genotypes_.length-2);
		//calculating SSx
		double SSx=0;
		for (int i=0; i<genotypes_.length; i++){
			SSx=SSx + Math.pow((genotypes_[i]-alleleMean_), 2);			
		}		
		double tStat=beta*Math.sqrt(SSx)/Math.sqrt(stdErr);
		double priorP = myT.cumulativeProbability(tStat);
		zscore_ = DistributionMethods.normalInverseCumulativeProbability(priorP);
		
		chi2Stat_=zscore_*zscore_;
		//pval_=1-chiSquared1df_.cumulativeProbability(chi2Stat_);
		pval_=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(chi2Stat_);
		
	}
	/** calculate new p-values and zscores via direct statistics. Only works if phenotype is normal(0,1) distributed
	 * Formula for zScore:=\sum((x-\sum(x)/n_x)*y)/sqrt(sum(x-sum(x)/n_x)^2).
	 * */
	public void recalcPvalsAndZscoresDirect(double[] phenotype) {
		
		if (phenotype.length!=genotypes_.length)
			throw new RuntimeException("phenotype vector not same length as genotype vector.");
		//calculate allele mean
		double alleleMean=0;
		for (int i=0; i<genotypes_.length; i++){
			alleleMean=alleleMean+genotypes_[i];
		}
		alleleMean=alleleMean/genotypes_.length;
				
		double zUnsc=0;		
		for (int i=0; i<genotypes_.length; i++){
			zUnsc=zUnsc + phenotype[i]*(genotypes_[i]-alleleMean_);		
		}
		//calculating SSx
		double SSx=0;
		for (int i=0; i<genotypes_.length; i++){
			SSx=SSx + Math.pow((genotypes_[i]-alleleMean_), 2);			
		}		
		zscore_=zUnsc/Math.sqrt(SSx);						
		chi2Stat_=zscore_*zscore_;
		//pval_=1-chiSquared1df_.cumulativeProbability(chi2Stat_);
		pval_=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(chi2Stat_);
		
	}
					
	
	// ----------------------------------------------------------------------------

	/** Set this snp as coding for the given gene */
	public void addCoding(String geneId) {
		
		if (coding_ == null)
			coding_ = new HashSet<String>();
		coding_.add(geneId);
	}

		
	// ----------------------------------------------------------------------------

	/** Return true if this snp is a coding mutation for another gene than the given gene */
	public boolean isCodingForOtherGene(String thisGene) {
		
		if (coding_ != null && !coding_.contains(thisGene)) {
			assert coding_.size() > 0;
			return true;
		} else
			return false;
	}


	// ----------------------------------------------------------------------------

	/** Get a string representation in BED format */
	public String bedFormatString() {
		
		String bed = super.UCSCbedFormatString();
		bed += "\t" + pval_;
		
		NumberFormat df = new DecimalFormat("0.######");
		
		// Minor allele frequency
		if (genotypes_ != null) {
			double maf = getMaf();
			assert(maf <= 1.0 && maf >= 0.0);
			bed += "\t" + df.format(maf);
		}
		return bed;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get a string representation in UCSC BED format */
	public String UCSCbedFormatString() {
		
		int start0 = start_ - 1 ;
		return chr_ + "\t" + start0 + "\t" + end_ + "\t" + id_ + "\t" + zscore_;
	}

	
	// ----------------------------------------------------------------------------

	/** Compute allele statistics (mean, standard deviation, number of non-missing genotypes, maf) */
	public void computeAlleleStats() {
		
		computeAlleleMeanAndNonMissing();	
		computeAlleleSd();	
		computeMinorAlleleFreq();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Parse the string array to initialize snp position and genotypes, return true if success 
	 * @throws DataInconsistencyException */
	public boolean parseGenotypeAndChrPos(String[] splittedStr, boolean tpedFormat) {		
		
		if (tpedFormat)
			setGenotypeAndChrPosFromPlinkTPed(splittedStr);
		else
			setGenotypeAndChrPosFromHomegrown(splittedStr);
		
		return genotypes_ != null;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Transform p-values to chi-squared statistics if it hasen't been done yet */ 
	public void computeChiSquaredStatistics() {
		
		// If chi2 stat has already been computed return
		if (!Double.isNaN(chi2Stat_))
			return;
		
		// This is only consistent with R down to p-values 1e-14 (see SnpTest), hence the Settings.minSnpPvalue=1e-14
				chi2Stat_ = DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(pval_);
		//chi2Stat_ = chiSquared1df_.inverseCumulativeProbability(1-pval_);
		// The inverse cumulative prob for p-val=1e-14 is 59.89...
		//f (chi2Stat_ > 59.90) // this should never happen if pval > minPval=1e-14
		//	throw new RuntimeException("P-values < 1e-14 are not supported (chi2 stat > 59.896088)");
	}
	
	
	// ----------------------------------------------------------------------------

	/** Convert from phased to unphased genotype format */
	public void dephase() {
		
		// Note, do not check for genotypePhasedFormat_ here, because when dephasing it's now set to
		// false before dephase() is called the first time
		
		byte[] unphasedGenotypes = new byte[genotypes_.length/2];
		//assert(!Math.mod(((double) genotypes_.length),2));
		
		for (int i=0; i < unphasedGenotypes.length; ++i){
			assert (!(genotypes_[i*2] != -9  && genotypes_[i*2 + 1] == -9));
			assert (!(genotypes_[i*2] == -9  && genotypes_[i*2 + 1] != -9));
			unphasedGenotypes[i] = (byte) (genotypes_[i*2] + genotypes_[i*2 + 1]);
			if (genotypes_[i*2] == -9)
				unphasedGenotypes[i] = -9;
		}
		genotypes_ = unphasedGenotypes;
	}

	
	// ----------------------------------------------------------------------------

	/** 
	 * Write genotype information of this snp to the binary file. Format is:
	 * 0. id_
	 * 1. genotypes_.length
	 * 2. genotypes_[]
	 * 3. maf_
	 * 4. alleleMean_
	 * 5. alleleSd_
	 * @throws IOException 
	 */
	public void writeGenotype(DataOutputStream os) throws IOException {
				
		// NOTE: CHANGE ALSO readBinary() AND skip() WHEN CHANGING ANYTHING HERE!
		os.writeUTF(id_);

		for (int i=0; i<genotypes_.length; i++)
			os.writeByte(genotypes_[i]);
		
		os.writeDouble(maf_);
		os.writeDouble(alleleMean_);
		os.writeDouble(alleleSd_);
	}


	// ----------------------------------------------------------------------------

	/** 
	 * Read this snp from the binary file. Format is:
	 * 1. genotypes_.length
	 * 2. genotypes_
	 * 3. maf_
	 * 4. alleleSd_
	 * 5. alleleMean_
	 * @throws IOException 
	 */
	public void readGenotype(DataInputStream is) throws IOException {
		
		// id_ is already read 
		genotypes_ = new byte[genotypeLength_];
		for (int i=0; i<genotypeLength_; i++)
			genotypes_[i] = is.readByte();
		
		maf_ = is.readDouble();
		alleleMean_ = is.readDouble();
		alleleSd_ = is.readDouble();
	}

	
	// ----------------------------------------------------------------------------

	/** Skip a SNP in the binary file */
	public static void skip(DataInputStream is) throws IOException {
		
		// Skip bytes
		int numSkipped = is.skipBytes(numBytes_);
		assert numSkipped == numBytes_;
	}
	
	// ----------------------------------------------------------------------------

	/** String representation in tped format */
	public String tpedString(){
		
		// Chromosome number
		String chrNr = chr_.substring(3);

		HashMap<Byte, String> toPlinkMap;
		if (genotypeIsPhased_)
			toPlinkMap = fromPhasedToPlinkMap_;
		else 
			toPlinkMap = fromDephasedToPlinkMap_;							

		String tpedGenotype = "";
		for (int genotype:genotypes_)						
			tpedGenotype += " " + toPlinkMap.get(genotype);	

		String line = chrNr + " " + id_ + " 0 " + start_ + tpedGenotype; 	
		return line;
	}		

	
	// ============================================================================
	// PRIVATE METHODS

	/** Compute allele mean and number of non-missing values */
	private void computeAlleleMeanAndNonMissing(){		
		
		int alleleSum = 0;
		int nonMissing = 0;
		for (int i = 0; i < genotypes_.length; i++) {
			if (genotypes_[i] != -9){
				alleleSum += genotypes_[i];
				++nonMissing;
			}
		}
		alleleMean_ = ((double) alleleSum)/nonMissing;
		alleleNonMissing_ = nonMissing;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Compute allele standard deviation */
	private void computeAlleleSd() {
		
		double alleleSumOfSquares = 0;
		int numberOfMissing = genotypes_.length - alleleNonMissing_;
		for (int i = 0; i < genotypes_.length; i++) 
			if (genotypes_[i] != -9)
				alleleSumOfSquares += genotypes_[i] * genotypes_[i];
		alleleSumOfSquares += Math.pow(alleleMean_, 2)*numberOfMissing;				
		double alleleVar = ((double) alleleSumOfSquares)/genotypes_.length - Math.pow(alleleMean_, 2);
		alleleSd_ = Math.sqrt(alleleVar);
	}	
	
	
	// ----------------------------------------------------------------------------

	/** Compute MAF */
	private void computeMinorAlleleFreq() {	
		
		double maf = 0;
		for (int i=0; i<genotypes_.length; i++) 
			if (genotypes_[i]!= -9)
				maf += genotypes_[i];
			//if (genotypes_[i] > 1)	
			//	Main.error("Check, not yet implemented");
		
		if (genotypeIsPhased_)
			maf /= alleleNonMissing_;
		else 
			maf /= (2*alleleNonMissing_);
		// Make sure it's the minor allele frequency
		if (maf > 0.5)
			maf = 1.0 - maf;
		assert maf <= 0.5 && maf >=0;
		
		maf_ = maf;
	}

		
	// ----------------------------------------------------------------------------

	/** Parse tped format */
	private void setGenotypeAndChrPosFromPlinkTPed(String[] splittedTpedStr) {
		
		assert(splittedTpedStr[1].equals(this.id_));
		
		byte[] genotypes = new byte[splittedTpedStr.length-4];		
		for (int i=0; i < genotypes.length ; ++i){		
			genotypes[i]=fromPlinkToPhasedMap_.get(splittedTpedStr[i+4]);			
		}	

		// Do not initialize position and genotype for snps with a constant genotype
		if (isConstant(genotypes))
			return;
		
		genotypes_ = genotypes;
		String chr = "chr" + splittedTpedStr[0];
		int start = Integer.parseInt(splittedTpedStr[3]);
		this.setPosition(chr, start, start+1, true);
		assert(Settings.dePhase_);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Parse homegrown format */
	private void setGenotypeAndChrPosFromHomegrown(String[] nextLine){
		
		assert(nextLine[0].equals(this.id_));
		
		String chr = nextLine[1];
		int start = Integer.parseInt(nextLine[2]);
		char minorAllele ='N';
		if (nextLine[4].length()==0){
			minorAllele='D';
		}
		else if (nextLine[4].length()>1){
			minorAllele='I';
		}				
		else {
			minorAllele =  nextLine[4].charAt(0);					
		}
		String[] genotypesStr = nextLine[nextLine.length-1].split(";");
		byte[] genotypes = new byte[genotypesStr.length];		
		for (int i=0; i < genotypesStr.length ; ++i){		
			genotypes[i]=Byte.parseByte(genotypesStr[i]);			
		}
		
		// Do not initialize position and genotype for snps with a constant genotype
		if (isConstant(genotypes))
			return;
		
		genotypes_ = genotypes;
		this.setPosition(chr, start, start+1, true);
		this.setMinorAllele(minorAllele);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Return true if the given vector is constant or length zero */
	private boolean isConstant(byte[] x) {

		boolean isConstant = true;
		for (int i=1; i<x.length; i++) {
			if (x[i] != x[0]) {
				isConstant = false;
				break;
			}
		}
		return isConstant;
	}

	public void copyPosAndAllele(Snp snp) throws DataInconsistencyException {
				
		copyPosition(snp);
		copyMinorMajorAllele(snp);
				
	}
	
	private void copyMinorMajorAllele(Snp snp) throws DataInconsistencyException {
		
		boolean snpHasBeenSeenInGWAS = false;
		if (minorAllele_!='N' || majorAllele_!='N'){
			snpHasBeenSeenInGWAS = true;
		}
		if (Settings.withZScore_ && minorAllele_!=snp.getMinorAllele()) {
			if(snp.getMinorAllele() !=majorAllele_ && snpHasBeenSeenInGWAS){
				throw new DataInconsistencyException("different minor allele of reference population not found GWAS data. Snp left out.");
		     }
			zscore_*=-1;
			char minorAlleleSummaryFile=minorAllele_;
			char majorAlleleSummaryFile=majorAllele_;
			minorAllele_=majorAlleleSummaryFile;
			majorAllele_=minorAlleleSummaryFile;	            
		}	
		
	}
	
	public void readPosAndMinorAllele(DataInputStream is) throws IOException,DataInconsistencyException {
		//TODO: Only solves homegrown case atm;
		// @David I changed this to an IllegalArgumentException because the other one was unknown on my system
			// id_ is already read 
			
		//String curChr = chr_;
			//int curStart = start_;
			//int curEnd = end_;
			
			chr_ = is.readUTF();		
			start_ = is.readInt();
			end_ = is.readInt();
			//if (curChr != null || curStart != -1 || curEnd != -1){
				//	if (!chr_.equals(curChr) || start_ != curStart || end_ != curEnd){
					//	throw new RuntimeException("snp seems to have been set before to another value");
					//}
			//}
			posStrand_ = is.readBoolean();
			minorAllele_ = is.readChar();
	}
	
public void readPosAndAllele(DataInputStream is) throws IOException,DataInconsistencyException {
	//TODO: Only solves homegrown case atm;
	// @David I changed this to an IllegalArgumentException because the other one was unknown on my system
		// id_ is already read 
		
	//String curChr = chr_;
		//int curStart = start_;
		//int curEnd = end_;
		
		chr_ = is.readUTF();		
		start_ = is.readInt();
		end_ = is.readInt();
		//if (curChr != null || curStart != -1 || curEnd != -1){
			//	if (!chr_.equals(curChr) || start_ != curStart || end_ != curEnd){
				//	throw new RuntimeException("snp seems to have been set before to another value");
				//}
		//}
		posStrand_ = is.readBoolean();
		char minorAllele = is.readChar();
		
		
		boolean snpHasBeenSeenInGWAS = false;
		if (minorAllele_!='N' || majorAllele_!='N'){
			snpHasBeenSeenInGWAS = true;
		}
		if (Settings.withZScore_ && minorAllele_!=minorAllele) {
			if(minorAllele!=majorAllele_ && snpHasBeenSeenInGWAS){
				throw new DataInconsistencyException("different minor allele of reference population not found GWAS data. Snp left out.");
		     }
			zscore_*=-1;
			char minorAlleleSummaryFile=minorAllele_;
			char majorAlleleSummaryFile=majorAllele_;
			minorAllele_=majorAlleleSummaryFile;
			majorAllele_=minorAlleleSummaryFile;	            
		}	
}
public void writePosAndAllele(DataOutputStream os) throws IOException {
	
	// NOTE: ALSO CHANGE readBinary() IF YOU CHANGE THIS
	os.writeUTF(id_);
	os.writeUTF(chr_);
	os.writeInt(start_);
	os.writeInt(end_);
	os.writeBoolean(posStrand_);
	os.writeChar(minorAllele_);
}
	
	// ============================================================================
	// GETTERS AND SETTERS

	public static void setGenotypeLength(int length) {
		genotypeLength_ = length;
		numBytes_ = 3*8 + genotypeLength_; // Three doubles taking 8 bytes each
	}
	public static int getGenotypeLength() { return genotypeLength_; }
	
	public static void setGenotypeIsPhased(boolean b) { genotypeIsPhased_ = b; }
	public static boolean getGenotypeIsPhased() { return genotypeIsPhased_; }
	public void setGenotypeIsPhased2(boolean b) { genotypeIsPhased2_ = b; }
	public boolean getGenotypeIsPhased2() { return genotypeIsPhased2_; }
	
	public double getPval() { return pval_; }
	public void setPval(double pval) { pval_ = pval; }
    public char getMinorAllele(){ return minorAllele_;}
    public char getMajorAllele(){ return majorAllele_;}
    
    public void setMinorAllele(char minorAllele){ minorAllele_ = minorAllele; }
    public void setMajorAllele(char majorAllele){ majorAllele_ = majorAllele; }

	public double getChi2Stat() {
		if (Double.isNaN(chi2Stat_)){
			throw new RuntimeException("asking for chi2Stat but seems to be not computed yet");			
		}
		return chi2Stat_; }
	public void setChi2Stat(double s) { chi2Stat_ = s; }

	public byte[] getGenotypes() { return genotypes_; }
	public void setGenotypes(byte[] g) { genotypes_ = g; }	

	public void setMaf(double maf) {maf_ = maf;}
	public double getMaf() { return maf_;}
	
	public double getZscore() {return zscore_;}
	public void setZscore(double zscore) {zscore_=zscore;}
	
	public double getAlleleMean(){ return alleleMean_;}
	public double getAlleleSd(){ return alleleSd_;}
	
			
}
