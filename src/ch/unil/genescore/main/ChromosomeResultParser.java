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
package ch.unil.genescore.main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;

import ch.unil.gpsutils.FileExport;



/**
 * Concatenate result files from individual chromosomes.
 */
public class ChromosomeResultParser {

	/** The files of this directory */
	private ArrayList<File> files_ = null;
	/** The concatenated result file */
	private FileExport writer_ = null;
	/** Set true to delete individual chromosome result files after concatenating them */
	private boolean deleteOriginals_ = false;
	
	/** The current prefix of result files that are being concatenated */
	private String curPrefix_ = null;
	/** Indicates whether the current result files that are being concatenated have a header */
	private String curHeader_ = null;
	/** Indicates whether the current files are expected to be comprehensive (all chromosomes) */
	private boolean curExpectedComplete_ = false; 
	/** List of files to delete in this directory */
	private ArrayList<File> dirDelete_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public ChromosomeResultParser() {

		deleteOriginals_ = Pascal.set.deleteOriginals_;
	}

	
	// ----------------------------------------------------------------------------

	/** Concatenate all chromosome result files */
	public void concatenateChromosomeResultFiles(String path) {
		
		File root = new File(path);
		File[] list = root.listFiles();
		files_ = new ArrayList<File>();

		if (list == null) 
			return;

		// Recursively traverse subdirectories
		for (File f : list) {
			if (f.isDirectory())
				concatenateChromosomeResultFiles(f.getAbsolutePath());
			else
				files_.add(f);
		}
		
		// Concatenate files for the current directory
		Pascal.println("Directory: " + root.getName());
		if (deleteOriginals_)
			dirDelete_ = new ArrayList<File>();
		else
			dirDelete_ = null;

		while (files_.size() > 0) {
			getChromosomeResultPrefix(files_.get(0));

			// Remove the file if it's not a chromosome result
			if (curPrefix_ == null) {
				files_.remove(0);
				continue;
			}
			
			// Concatenate the files of this type
			concatenateRuns();
		}
		
		// Delete original files
		if (dirDelete_ != null)
			for (File f : dirDelete_)
				f.delete();
		
		Pascal.println();
	}

	

	// ============================================================================
	// PRIVATE FUNCTIONS

	/** Get the prefix of this chromosome result file, return null if this is not a chromosome result (*.chrX.*) */
	private void getChromosomeResultPrefix(File file) {
				
		String filename = file.getAbsolutePath();

		if (!filename.matches(".*\\.chr\\d{1,2}\\.txt$")) {
			curPrefix_ = null;
			return;
		}
		
		if (filename.contains(".genescores.")) {
			curHeader_ = getHeader(file);
			curExpectedComplete_ = true;
		} else if (filename.contains(".snps.")) {
			curHeader_ = null;
			curExpectedComplete_ = true;
		} else if (filename.contains(".aboveSnpLimit.") || filename.contains(".corNotPositiveDefinite.") || filename.contains(".scoreComputeError.")) {
			curHeader_ = getHeader(file);
			curExpectedComplete_ = false;
		} else {
			curPrefix_ = null;
			return;
		}
		
		curPrefix_ = filename.replaceFirst("\\.chr\\d{1,2}\\.txt$", "");
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the header from the given file */
	private String getHeader(File file) {
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String header = reader.readLine();
			reader.close();
			return header;
			
		} catch (Exception e) {
			Pascal.error(e);
			return null;
		}
	}

	
	// ----------------------------------------------------------------------------

	/** Concatenate result files with the given prefix */
	private void concatenateRuns() {
		
		File outfile = new File(Pascal.set.outputDirectory_, curPrefix_ + ".txt");
		writer_ = new FileExport(Pascal.log, outfile);
		
		// Write header
		if (curHeader_ != null)
			writer_.println(curHeader_);
		
		ArrayList<Integer> missingChromosomes = new ArrayList<Integer>();
		
		for (int i=1; i<=22; i++) {
			String filename = curPrefix_ + ".chr" + i + ".txt";
			File f = getFile(filename);
			
			if (f == null) {
				missingChromosomes.add(i);
			} else { 
				appendFile(f);
				if (dirDelete_ != null)
					dirDelete_.add(f);
			}
		}
		
		// Display warning with missing chromosomes
		if (curExpectedComplete_ && missingChromosomes.size() > 0) {
			String msg = "Missing chromosomes: " + missingChromosomes.get(0);
			for (int i=1; i<missingChromosomes.size(); i++)
				msg += ", " + missingChromosomes.get(i);
			Pascal.warning(msg);
			
			// Do not delete any files in this directory if there were incomplete runs
			dirDelete_ = null;
		} 
		
		writer_.close();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the file of the given name from files_, return null if not found */
	private File getFile(String filename) {
		
		Iterator<File> iter = files_.iterator();
		while (iter.hasNext()) {
			File f = iter.next();
		
			if (f.getAbsolutePath().equals(filename)) {
				iter.remove();
				return f;
			}
		}
		return null;
	}

	
	// ----------------------------------------------------------------------------

	/** Append the given file to the current file writer */
	private void appendFile(File file) {
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			System.out.println(file.getAbsolutePath());
			if (curHeader_ != null) {
				String header = reader.readLine();
				if (!header.equals(curHeader_))
					Pascal.error("Chromosome result files with inconsistent header lines");
			}
			
			// Read the first line
			String nextLine = reader.readLine();
			
			// Copy each line to the new file
			while (nextLine != null) {
				writer_.println(nextLine);
				nextLine = reader.readLine();
			}
			reader.close();
						
		} catch (Exception e) {
			Pascal.error(e);
		}
	}



}
