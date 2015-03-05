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
package ch.unil.genescore.vegas;

import java.io.DataInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;


public class ReferencePopulationFakeSignal  {
//	public class ReferencePopulationFakeSignal {
		
	private ReferencePopulation refpop_ = new ReferencePopulation();		
	
	Snp currentSnp_ = null;
	Random rnd = new Random(Settings.randomSeed_);
	ArrayList<String> snpNames_= new ArrayList<String>();
	double[] signal_ = null;
	
	
	private void initializeChr(String chr) {
		
		refpop_.closeDataInputStream();
				
		// Open input binary file
		String filename = refpop_.getGenotypeFileHandler().getFilePrefix() + chr + ".gnt.ser.gz";
		refpop_.openDataInputStream(filename);
		
		// Read genotype vector length and phased flag
		try {
			int length = refpop_.getInStream().readInt();
			boolean phased = refpop_.getInStream().readBoolean();
			
			Snp.setGenotypeLength(length);
			Snp.setGenotypeIsPhased(phased);
			if (signal_ == null){				
				signal_ = new double[length];
			}
			
		} catch (IOException e) {
			Main.error(e, "Error reading genotype length");
		}
		if (Snp.getGenotypeIsPhased()==true){
			
			throw new RuntimeException("problem in method generating fake phenotype with signal: saved genotype is phased but has not been treated yet");
		}
		
	}
	public void runFakeSignal(){
		
		for (String chr : refpop_.chromosomes_){			
			runChr(chr);
		}
			
	}
	
	private void runChr(String chr){
		
		initializeChr(chr);
		while(!refpop_.isEndOfFile()){
			loadSnpGenotype();
			if(!refpop_.isEndOfFile())
				processSnp();
		}
	}
	
	/** Load the next snp from the input stream */
	private void loadSnpGenotype() {
		
		try {
			String snpId = refpop_.getInStream().readUTF();			
			// End of file is marked with the version ID
			if (snpId.equals(refpop_.getGenotypeFileHandler().getBinaryFileVersionID())) {
				refpop_.setEndOfFile(true);
				return;
			}						
			currentSnp_ = new Snp(snpId);
			currentSnp_.readGenotype(refpop_.getInStream());					
		} catch (Exception e) {
			Main.error(e, "Error loading genotypes");
		}		
	}
	
	/** Load the next snp from the input stream */
	private void processSnp() {
		double maf;
		double rndNr;
		double beta;
		double sd;
		byte[] genotypes = null;
		currentSnp_.computeAlleleStats();		
		maf = currentSnp_.getMaf();
		sd = currentSnp_.getAlleleSd();
		if (maf > Settings.useMafCutoff_ || sd==0){
			rndNr = rnd.nextDouble();
			if (rndNr < Settings.chanceOfSignal_){
				beta = 10;//rnd.nextGaussian()*10;
				genotypes = currentSnp_.getGenotypes();
				assert(genotypes.length==signal_.length);
				for (int i=0; i<signal_.length; i++){
					signal_[i]=signal_[i]+beta*genotypes[i];
				}
				snpNames_.add(currentSnp_.id_);
			}
		}				
	}
	
	public double[] getSignal(){
		return signal_;
	}
	
	
	
}