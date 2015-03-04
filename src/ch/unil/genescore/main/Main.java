/*
Copyright (c) 2013 Daniel Marbach

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper available at:
http://networkinference.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
package ch.unil.genescore.main;

import java.io.File;
import java.lang.reflect.Field;
import java.util.LinkedList;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.gene.GeneAnnotation;
import ch.unil.genescore.pathway.PathwayMain;
import ch.unil.genescore.prioritization.PrioritizationMain;
import ch.unil.genescore.topLdSnps.TopLdSnpsMain;
import ch.unil.genescore.vegas.GenomeWideScoring;
import ch.unil.genescore.vegas.ReferencePopulation;


/**
 * Main class
 */
public class Main {


	// ============================================================================
	// STATIC METHODS

	/** Main function */
	static public void main(String[] args) {
		
		try {
			Main main = new Main(args);
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
	
	/** Constructor, parse command-line arguments, initialize settings */
	public Main(String[] args) {
		
		Main.println("SETTINGS FILE");
		Main.println("-------------\n");

		// Parse command-line arguments and initialize settings
		GeneScoreOptionParser optionParser = new GeneScoreOptionParser();
		optionParser.parse(args);
		optionParser.firstProcessingOfSettings();
				
		// set jna path
		System.setProperty("jna.library.path", "jars/lib/");
		// Create output directory
		File outputDir = new File(Settings.outputDirectory_);
		if (!outputDir.exists())
			outputDir.mkdirs();
	}


    // ----------------------------------------------------------------------------

	/** Parse the command-line arguments, read the files, perform network inference, write outputs */
	public void run() {
		
		if (!Settings.writeUsedSettings_.equals("")){
			writeUsedSettings();
		}
		if (Settings.runTopLdSnp_)
			runTopLdSnps();
		else if (Settings.runConcatenateChromosomeResults_)
			concatenateChromosomeResults();
		else if (Settings.runPathwayAnalysis_)
			runPathwayAnalysis();				
		else if (Settings.runPrioritizationAnalysis_)
			runPrioritizationAnalysis();
		else
			computeGeneScores();			
		System.out.println("Success!");
		
	}

	
    // ----------------------------------------------------------------------------

	/** Compute gene scores summarizing multiple snps */
	public void computeGeneScores() {
		
		Main.println("LOADING INPUT FILES");
		Main.println("-------------------\n");		

		GenomeWideScoring geneScore = new GenomeWideScoring();
		System.out.println(Settings.ucscAnnotationFile_);
		// Load genes
		LinkedList<Gene> genes = GeneAnnotation.createAnnotationInstance().loadAnnotation(Settings.genesToBeLoadedFile_);
		// Load snps

		geneScore.setGenes(genes);
		ReferencePopulation myRefPop=new ReferencePopulation();				
		myRefPop.loadGwasAndRelevantSnps();
		geneScore.setReferencePopulation(myRefPop);
	
		
				
		// Print info
		Main.println();
		if (Settings.chromosome_.equals(""))
			Main.println("Loaded all chromosomes:");
		else
		Main.println("Loaded chromosome " + Settings.chromosome_ + ":");
		Main.println("- " + genes.size() + " genes ");
		Main.println();

		Main.println("COMPUTING GENE SCORES");
		Main.println("---------------------\n");		
		
		geneScore.computeScores();
	}


    // ----------------------------------------------------------------------------

	/** Run pathway / gene set enrichment analysis */
	public void runPathwayAnalysis() {

		PathwayMain pathways = new PathwayMain();
		pathways.run();
	}
	
	public void runPrioritizationAnalysis() {

		PrioritizationMain prioritize = new PrioritizationMain();
		prioritize.run();
	}

	public void runTopLdSnps(){

		TopLdSnpsMain topLd = new TopLdSnpsMain();
		topLd.run();
	}


	// ----------------------------------------------------------------------------

	/** Compute gene scores summarizing multiple snps */
	public void concatenateChromosomeResults() {

		// Concatenate chromosome results and exit
		Main.println("CONCATENATING CHROMOSOME RESULTS");
		Main.println("--------------------------------\n");	
	
		ChromosomeResultParser results = new ChromosomeResultParser();
		results.concatenateChromosomeResultFiles(Settings.concatenateChromosomeResultsDir_);
	}
	
	public void writeUsedSettings(){
		FileExport fl = new FileExport(Settings.writeUsedSettings_, false);
		Class cls = null;
		try {
			
			cls = Class.forName("ch.unil.genescore.main.Settings");
		
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
	// PRIVATE METHODS
		
	
	}
