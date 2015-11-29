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

import java.lang.reflect.Field;
import java.util.LinkedList;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.pathway.PathwayMain;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.gpsutils.FileExport;
import ch.unil.gpsutils.Logger;
import ch.unil.gpsutils.Utils;


/**
 * Main class
 */
public class Pascal {

	/** The logger -- a different logger can be plugged in for custom logging */
	static public Logger log;
	/** The settings */
	static public PascalOptionParser set;
	/** The utilities */
	static public Utils utils;

	/** Computes gene scores */
	private GenomeWideScoring genomeWideScoring;

	
	// ============================================================================
	// STATIC METHODS

	/** Main function */
	static public void main(String[] args) {
		
		try {
			Pascal main = new Pascal(args);
			main.run();
			
		} catch (Exception e) {
			error(e);
		}
	}

	
    // ----------------------------------------------------------------------------

	/** Print the stack trace of the exception and exit */
	static public void error(Exception e) {
		e.printStackTrace();
		System.exit(-1); // return -1 in case of error
	}

	/** Print the stack trace of the exception and exit with error message */
	static public void error(Exception e, String msg) {
		System.err.println("ERROR: " + msg);
		error(e);
	}

	/** Print error message and exit */
	static public void error(String msg) {
		System.err.println("ERROR: " + msg);
		System.exit(-1);
	}

	/** Print warning message to stderr */
	static public void warning(String msg) {
		System.out.println("WARNING: " + msg);
		System.out.flush();
	}

	/** Write line to stdout */
	static public void println(String msg) {		
		System.out.println(msg);
	}

	/** Write empty line to stdout */
	static public void println() {		
		System.out.println();
	}

	/** Write string to stdout */
	static public void print(String msg) {
		System.out.print(msg);
	}

	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor with custom logger */
	public Pascal(String[] args) {
		this(args, null);
	}

	/** Constructor with custom logger */
	public Pascal() {
		this(null, null);
	}

	/** Constructor with custom logger */
	public Pascal(Logger customLog) {
		this(null, customLog);
	}

	/** Constructor, parse command-line arguments, initialize settings */
	public Pascal(String[] args, Logger customLog) {
		
		// Initialize
		if (customLog != null)
			log = customLog;
		else
			log = new Logger(); // must be first
		utils = new Utils(log);

		// Set defaults
		set = new PascalOptionParser();
		
		// Parse command-line arguments and initialize settings
		if (args != null)
			set.parse(args);
				
		// Set verbose flag in logger
		log.setVerbose(set.verbose_);
		// set jna path
		System.setProperty("jna.library.path", "jars/lib/");
		// Create output directory
		set.outputDirectory_.mkdirs();
	}


    // ----------------------------------------------------------------------------

	/** Parse the command-line arguments, read the files, perform network inference, write outputs */
	public void run() {
		
		if (!set.writeUsedSettings_.equals("")){
			writeUsedSettings();
		}
	
		if (set.runConcatenateChromosomeResults_)
			concatenateChromosomeResults();
		else if (set.runPathwayAnalysis_)
			runPathwayAnalysis();				
		else
			computeGeneScores();			
		System.out.println("Success!");
		
	}

	
    // ----------------------------------------------------------------------------

	/** Compute gene scores summarizing multiple snps */
	public void computeGeneScores() {
		
		Pascal.println("LOADING INPUT FILES");
		Pascal.println("-------------------\n");		

		genomeWideScoring = new GenomeWideScoring();
		System.out.println(set.ucscAnnotationFile_);
		// Load genes
		LinkedList<Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation(set.genesToBeLoadedFile_);
		// Load snps

		genomeWideScoring.setGenes(genes);
		ReferencePopulation myRefPop=new ReferencePopulation();				
		myRefPop.loadGwasAndRelevantSnps();
		genomeWideScoring.setReferencePopulation(myRefPop);
	
		
				
		// Print info
		Pascal.println();
		if (set.chromosome_.equals(""))
			Pascal.println("Loaded all chromosomes:");
		else
		Pascal.println("Loaded chromosome " + set.chromosome_ + ":");
		Pascal.println("- " + genes.size() + " genes ");
		Pascal.println();

		Pascal.println("COMPUTING GENE SCORES");
		Pascal.println("---------------------\n");		
		
		genomeWideScoring.computeScores();
	}


    // ----------------------------------------------------------------------------

	/** Run pathway / gene set enrichment analysis */
	public void runPathwayAnalysis() {

		PathwayMain pathways = new PathwayMain(this);
		pathways.run();
	}
	
	


	// ----------------------------------------------------------------------------

	/** Compute gene scores summarizing multiple snps */
	public void concatenateChromosomeResults() {

		// Concatenate chromosome results and exit
		Pascal.println("CONCATENATING CHROMOSOME RESULTS");
		Pascal.println("--------------------------------\n");	
	
		ChromosomeResultParser results = new ChromosomeResultParser();
		results.concatenateChromosomeResultFiles(set.concatenateChromosomeResultsDir_);
	}
	
	public void writeUsedSettings(){
		FileExport fl = new FileExport(log, set.writeUsedSettings_);

		try {
			Class<?> cls = Class.forName("ch.unil.genescore.main.Settings");
	        System.out.println("Class found = " + cls.getName());
	        System.out.println("Package = " + cls.getPackage());
	        Field f[] = cls.getFields();
	        for (int i = 0; i < f.length; i++) {
	        	 String result = String.format("%s\t%s",f[i].getName() ,f[i].get(null));
		          fl.println(result);
	        }
        } catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
        }		
		fl.close();
	}
	
	
	// ============================================================================
	// GETTERS AND SETTERS

	public GenomeWideScoring getGenomeWideScoring() {
		return genomeWideScoring;
	}
		
	
}
