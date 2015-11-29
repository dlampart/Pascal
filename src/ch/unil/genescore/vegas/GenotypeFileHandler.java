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

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.zip.GZIPOutputStream;

import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileParser;

public class GenotypeFileHandler {

	/** The path and prefix of the ref pop files */
	protected String filePrefix_ = null;
	/** Version of the binary file format */
	protected static final String binaryFileVersionID_ = "refpop_v1";
	/** Flag set true if at least one file had to be converted to binary (useful for console output) */
	protected boolean missingBinaryFilesCreated_ = false;


	/** Check if binary file for the given chromosome exist, if not create it */
	protected void writeBinaryFiles(String chr) {		

		String prefix = getFilePrefix() + chr + ".";
		File binaryPosFile = new File(prefix + "pos.ser.gz");
		File binaryGntFile = new File(prefix + "gnt.ser.gz");

		// If the binary files don't exist yet
		if (!binaryPosFile.exists() || !binaryGntFile.exists()) {
			if (Pascal.set.refPopFileExtension_.equals("ser.gz"))
				throw new RuntimeException("The binary reference population files specified in settings don't exist: " + prefix + "*");

			// Print info to console
			if (!missingBinaryFilesCreated_) {
				Pascal.println("Creating missing binary reference population files. Hang on, it has only to be done this one time...");
				missingBinaryFilesCreated_ = true;
			}

			// Create binary files
			File settingsFile = new File(prefix + Pascal.set.refPopFileExtension_);			
			writeBinaryFiles(settingsFile, binaryPosFile, binaryGntFile);
		}
	}
	static class UnserializedSnpFile {


		int snpIdIndex = -1;
		boolean tpedFormat;
		FileParser parser_=null;
		boolean dephase_;
		/** responsible for loading the input genotype file that hasn't been serialized yet.*/
		public UnserializedSnpFile(){
		}
		public void setupFileToRead(File inputFile){
			String columnSeparator=null;
			// Detect format based on file extension
			if (inputFile.getName().endsWith("txt.gz")) {
				columnSeparator = "\t";
				snpIdIndex = 0;
				// TODO autodetect
				Snp.setGenotypeIsPhased(true);
				tpedFormat = false;

			} else if (inputFile.getName().endsWith("tped.gz")) {
				columnSeparator = "\\s+";
				snpIdIndex = 1;
				Snp.setGenotypeIsPhased(true);
				tpedFormat = true;

			} else {
				throw new RuntimeException("Reference population text files must have extension '.txt.gz' or '.tped.gz'");
			}

			Snp.setGenotypeIsPhased(Pascal.set.dePhase_ || getTpedFormat());			 				
			parser_ = new FileParser(Pascal.log, inputFile);
			parser_.setSeparator(columnSeparator);
			dephase_=setDephase();
		}

		private boolean setDephase(){
			boolean dephase = Pascal.set.dePhase_ || getTpedFormat();
			if (dephase) {				
				// Note, snps will be phased when loaded, but are then dephased before any computation is done (see below),
				// so in the end they will be in unphased format. This has to be set false here, otherwise snp.computeAlleleStats()
				// below will think the genotypes are still phased
				Snp.setGenotypeIsPhased(false);			
			}
			return dephase;
		}


		public boolean snpAvailable(){
			return parser_.getCurrentLine() != null;
		}
		/**returns next snp from handler. Will return null if genotype is constant*/
		public Snp getNextSnp() {			
			if(!snpAvailable()){
				throw new RuntimeException("snp not available.");
			}
			String[] nextLine = parser_.readLine();			
			if(nextLine[0]==null){
				System.out.println("asdf");
			}
			String snpId = nextLine[getSnpIdIndex()];
			Snp snp = new Snp(snpId);
			if (nextLine.length < 4)
				parser_.error("Expected at least 4 columns");	
			// Initialize position and genotype, skip if genotype is constant 
			if (!snp.parseGenotypeAndChrPos(nextLine, getTpedFormat()))
				return null;				
			if (dephase_){
				Snp.setGenotypeIsPhased(false);	
				snp.setGenotypeIsPhased2(false);	
				snp.dephase();
			}
			return snp;		 
		}

		public  boolean getTpedFormat(){
			return tpedFormat;
		}

		public int getSnpIdIndex(){
			return snpIdIndex;
		}

		public FileParser getParser(){
			return parser_;
		}
	}

	/** Create a binary file based on the given input text file */
	private void writeBinaryFiles(File inputFile, File binaryPosFile, File binaryGntFile) {

		UnserializedSnpFile unserial = new UnserializedSnpFile();
		unserial.setupFileToRead(inputFile);

		try {
			// Open input text file
			FileParser parser = unserial.getParser();
			// Open output files
			DataOutputStream posStream = openDataOutputStream(binaryPosFile);
			DataOutputStream gntStream = openDataOutputStream(binaryGntFile);

			boolean isFirstSnp = true;		
			//int numGenotypes = -1;
			int curChr = -1;
			int prevSnpPos = -1;

			// Check if genotypes should be dephased
			// For each snp			
			while (unserial.snpAvailable()) {

				Snp snp = unserial.getNextSnp();
				if (snp==null){
					continue;
				}									
				// Check format
				if (isFirstSnp){					

					curChr = snp.getChrInt();					
					isFirstSnp = false;

					// Write genotype length to the binary file
					int length = snp.getGenotypes().length;
					Snp.setGenotypeLength(length);
					gntStream.writeInt(length);
					gntStream.writeBoolean(Snp.getGenotypeIsPhased());

				} else {			

					if (Snp.getGenotypeLength() != snp.getGenotypes().length)
						Pascal.error("error during file-Parsing: not same number of genotypes");
					if (curChr != snp.getChrInt())
						Pascal.error("error during file-Parsing: not constant chromosome");
					if (snp.start_ < prevSnpPos)
						Pascal.error("error during file-Parsing: SNPs must be ordered by position");
				}
				prevSnpPos = snp.start_;

				// Compute allele mean, standard deviation, frequency
				snp.computeAlleleStats();

				// Write to output stream
				snp.writePosAndAllele(posStream);
				snp.writeGenotype(gntStream);
			}

			// End of file is marked with the version ID
			posStream.writeUTF(getBinaryFileVersionID());
			gntStream.writeUTF(getBinaryFileVersionID());
			// Close files
			posStream.close();
			gntStream.close();
			parser.close();

		} catch (IOException e) {
			Pascal.error(e, "IO error converting text to binary file");
		}
	}


	/**get inputStream-Object that that iterates over snp positions*/
	public SnpSerializedPositionStream getSnpSerializedPositionStream(String chr){		 

		String filename = getFilePrefix() + chr + ".pos.ser.gz";				
		return new SnpSerializedPositionStream(filename, getBinaryFileVersionID());	

	}

	/**get inputStream-Object that that iterates over snp positions*/
	public SnpPositionStream getSnpPositionStream(String chr){		 

		String filename = getFilePrefix() + chr + ".pos.ser.gz";				
		return new SnpPositionStream(filename, getBinaryFileVersionID());	

	}
	/**get inputStream-Object that that iterates over snp positions*/
	public WrappedSnpPositionStream getSnpPositionStream(String[] chrs){				 
		//String[] filenames= new String[chrs.length];
		ArrayList<String> filenames= new ArrayList<String>();
		for (int i=0; i<chrs.length;i++){
			String curString = getFilePrefix() + chrs[i] + ".pos.ser.gz";
			filenames.add(curString);
		}
		Collections.sort(filenames);
		return new WrappedSnpPositionStream(filenames, getBinaryFileVersionID());	

	}

	/** Open the data output stream for the given chromosome */
	public DataOutputStream openDataOutputStream(File file)  {

		Pascal.println("Writing file: " + file.getPath());

		try {
			FileOutputStream outfile = new FileOutputStream(file);
			GZIPOutputStream gzip = new GZIPOutputStream(outfile);
			BufferedOutputStream buf = new BufferedOutputStream(gzip);
			DataOutputStream outStream = new DataOutputStream(buf);

			// Write the version, used as a check when reading files
			outStream.writeUTF(getBinaryFileVersionID());

			return outStream;

		} catch (Exception e) {
			throw new RuntimeException("Could not open binary output file: " + file);
		}
	}

	public GenotypeFileHandler(){		
		setupFilePrefix();
	}

	private void setupFilePrefix(){
		filePrefix_ = Pascal.set.refPopDirectory_ + "/" + Pascal.set.refPopFilePrefix_ + ".";
	}

	public String getFilePrefix(){
		return filePrefix_;
	}
	public String getBinaryFileVersionID(){return binaryFileVersionID_;}




}
