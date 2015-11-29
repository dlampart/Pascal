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

import java.io.File;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

/**
 * Defines and parses the command-line arguments
 */
public class PascalOptionParser extends PascalSettings {

	/** The option parser */
	OptionParser parser_;
	/** The options */
	OptionSet options;

	// ============================================================================
	// PUBLIC METHODS

	/** Parse the command-line arguments, initialize all settings */
	public PascalOptionParser() {

		// Defines the arguments
		defineArgs();
	}

	
	// ----------------------------------------------------------------------------

	/** Parses the command-line arguments, which were defined by defineArgs() */
	public void parse(String[] args) {

		// (1) Parse the options
		options = null;
		try {
			options = parser_.parse(args);
		} catch (Exception e) {
			displayHelp();
			throw new RuntimeException(e);
		}

		// Display help
		if (options.has("help") || options.has("h") || options.has("?")) {
			displayHelp();
			System.exit(0);
		}

		// (2-3) Set and load the settings file
		if (options.has("set"))
			loadSettings((String) options.valueOf("set"));

		// (4) Set command-line options (override settings in loaded settings file)

		if (options.has("outdir"))
			outputDirectory_ = getFileOption("outdir"); 
		if (options.has("outsuffix"))
			outputSuffix_ = (String) options.valueOf("outsuffix");		
		if (options.has("writesnpbed"))
			writeSnpBedFile_ = true;
		if (options.has("ucscdir"))
			ucscAnnotationFile_ = getFileOption("ucscdir");
		if (options.has("refpopdir"))
			refPopDirectory_ = getFileOption("refpopdir");
		if (options.has("genestobeloaded"))
			genesToBeLoadedFile_ = getFileOption("genestobeloaded");
		if (options.has("pval"))
			snpPvalFile_ = getFileOption("pval");
		if (options.has("pvalcol"))
			pvalCol_ = (Integer) options.valueOf("pvalcol");
		if (options.has("mafcutoff"))
			useMafCutoff_ = (Double) options.valueOf("mafcutoff");
		
		if (options.has("pop"))
			refPopFilePrefix_ = (String) options.valueOf("pop");
		if (options.has("popformat"))
		refPopFileExtension_ = (String) options.valueOf("popformat");
		if (options.has("chr"))
			chromosome_ = "chr" + options.valueOf("chr");
		
		if (options.has("maxsnp"))
			maxSnpsPerGene_ = (Integer) options.valueOf("maxsnp");
		
		if (options.has("usefakephenotype"))
			useFakePhenotype_ = true;
		if (options.has("usefakesignal"))
			useFakeSignal_ = true;
		if (options.has("chanceofsignal"))
			chanceOfSignal_ = (Double) options.valueOf("chanceofsignal");
		if (options.has("multipleofphenotype"))
			multipleOfPhenotype_ = (Integer) options.valueOf("multipleofphenotype");
		
		
		if (options.has("useimhof"))
			useImhof_=true;		
		if (options.has("usedavies"))
			useDavies_=true;
		if (options.has("usefarebrother"))
			useFarebrother_=true;
		
		if (options.has("genewisesnpdir"))
			writeGenewiseSnpFiles_= (Boolean) options.valueOf("genewisesnpdir");
		
		if (options.has("analyticvegas"))
			useAnalyticVegas_ = true;
		if (options.has("maxvegas"))
			useMaxVegas_ = true;
		if (options.has("maxeffvegas"))
		useMaxEffVegas_ = true;

		if (options.has("withzscore"))
			withZScore_ = true;
		if (options.has("prunelist"))
			snpFilterFile_ = getFileOption("filterlist");
		if (options.has("up"))
			geneWindowUpstream_ = (Integer) options.valueOf("up");
		if (options.has("down"))
			geneWindowDownstream_ = (Integer) options.valueOf("down");
		//Pathway options
		if (options.has("annotation"))
			genomeAnnotation_ = (String) options.valueOf("annotation");
		if (options.has("gencodefile"))
			gencodeAnnotationFile_ = getFileOption("gencodefile");
		
		if (options.has("finalize")) {
			concatenateChromosomeResultsDir_ = (String) options.valueOf("finalize");
			runConcatenateChromosomeResults_ = true;
		}
		if (options.has("delete"))
			deleteOriginals_ = true;
		
		//pathway options
		if (options.has("runpathway"))
			runPathwayAnalysis_= true;
		
		if (options.has("mergedistance"))
			mergeGenesDistance_ = (Double) options.valueOf("mergedistance");
		
		if (options.has("genesetfile"))
			geneSetFile_ = getFileOption("genesetfile");
		
		if (options.has("hypgeom"))
			useHypGeom_ = true;
		
		if (options.has("simulated"))
			useSimulation_ = true;
		
		if (options.has("loadscores"))
			loadScoresFromFiles_= true;
		
		if (options.has("scores"))
			geneScoreFile_ = getFileOption("scores");
		
		if (options.has("metascores"))
			metaGeneScoreFile_ = getFileOption("metascores");
		if (options.has("deflationrate"))
			deflationRate_ = (Double) options.valueOf("deflationrate");
		
		if (options.has("detailed"))
			writeDetailedOutput_ = true;
		if (options.has("rewritesettings"))
			writeUsedSettings_ = getFileOption("rewritesettings");

		checkOptions();
	}

	// ============================================================================
	// PRIVATE METHODS

	/** Display help on console */
	static private void displayHelp() {

		System.out.println("USAGE");
		System.out.println("   java -jar Ngsea.jar <file> -o <file> [OPTIONS]");
		System.out.println("OPTIONS");
		System.out.println("   --help          Display this usage information");
		System.out.println("   TBD");
	}

	
	// ----------------------------------------------------------------------------

	/** Defines the command-line arguments */
	private void defineArgs() {

		// The command-line parser
		parser_ = new OptionParser();

		// Display help
		parser_.accepts("help");
		// settingsFile_
		parser_.accepts("set").withRequiredArg();
		parser_.accepts("writesnpbed");
		parser_.accepts("genewisesnpdir").withRequiredArg();
		// outputDirectory_
		parser_.accepts("outdir").withRequiredArg();
		parser_.accepts("outsuffix").withRequiredArg();
				
		// snpPvalFile_
		parser_.accepts("pval").withRequiredArg();
		// snpPvalCol_
		parser_.accepts("pvalcol").withRequiredArg();
		// population_
		parser_.accepts("pop").withRequiredArg();
		parser_.accepts("refpopdir").withRequiredArg();
		// populationFormat_
		parser_.accepts("popformat").withRequiredArg();
		parser_.accepts("gencodefile").withRequiredArg();
		parser_.accepts("genestobeloaded").withRequiredArg();
		// chromosome_
		parser_.accepts("chr").withRequiredArg().ofType(Integer.class);
		
		// testStatisticNumSnps_
		parser_.accepts("nsnp").withRequiredArg().ofType(Integer.class);
		// maxSnpsPerGene_
		parser_.accepts("maxsnp").withRequiredArg().ofType(Integer.class);
		parser_.accepts("mafcutoff").withRequiredArg().ofType(Double.class);
		
		parser_.accepts("usedavies");
		parser_.accepts("usefarebrother");
		parser_.accepts("useimhof");
		parser_.accepts("ucscdir").withRequiredArg();			

		// useAnalyticVegas_				
		parser_.accepts("analyticvegas");
		parser_.accepts("maxvegas");
		parser_.accepts("maxeffvegas");
		parser_.accepts("analyticvegasweighted");
		parser_.accepts("analyticvegasweightedwrapped");
		parser_.accepts("orthostat");
		parser_.accepts("orthostatwrapped");
		parser_.accepts("orthosum");
		parser_.accepts("projectionvegas");
		parser_.accepts("usefakephenotype");
		parser_.accepts("usefakesignal");
		parser_.accepts("chanceofsignal").withRequiredArg().ofType(Double.class);
		parser_.accepts("multipleofphenotype").withRequiredArg().ofType(Integer.class);		
			
		parser_.accepts("randseed").withRequiredArg().ofType(Integer.class);		
		
		parser_.accepts("bedfilepath").withRequiredArg();
		parser_.accepts("bedweight").withRequiredArg().ofType(Double.class);
		parser_.accepts("filteronbed");
		parser_.accepts("bedbackgroundweight").withRequiredArg().ofType(Double.class);
		parser_.accepts("bedbackgroundextension").withRequiredArg().ofType(Integer.class);
		parser_.accepts("conditionfraction").withRequiredArg().ofType(Double.class);		
		parser_.accepts("conditionfractioncross").withRequiredArg().ofType(Double.class);		
		parser_.accepts("explainedfraction").withRequiredArg().ofType(Double.class);
		parser_.accepts("projectionmafcutoff").withRequiredArg().ofType(Double.class);
		parser_.accepts("useonlygwassnps");
		parser_.accepts("withzscore");
		parser_.accepts("prunelist").withRequiredArg();
		// geneWindowUpstream_
		parser_.accepts("up").withRequiredArg().ofType(Integer.class);
		parser_.accepts("down").withRequiredArg().ofType(Integer.class);		
		parser_.accepts("maxpruningcutoff").withRequiredArg().ofType(Integer.class);
		parser_.accepts("annotation").withRequiredArg();
		parser_.accepts("rewritesettings").withRequiredArg();
			// concatenateChromosomeResultsDir
			parser_.accepts("finalize").withRequiredArg();
		// deleteOriginals_
		parser_.accepts("delete");
		//pathway options
		parser_.accepts("runpathway");
		parser_.accepts("mergedistance").withRequiredArg().ofType(Double.class);
		//path to pathway
		parser_.accepts("hypgeom");
		parser_.accepts("hypgeomquantiles").withRequiredArg();
		parser_.accepts("simulated");
		parser_.accepts("loadscores");
		parser_.accepts("scores").withRequiredArg();
		parser_.accepts("metascores").withRequiredArg();
		parser_.accepts("deflationrate").withRequiredArg().ofType(Double.class);
		
		parser_.accepts("genesetfile").withRequiredArg();
		parser_.accepts("detailed");
		parser_.accepts("eqtlfile").withRequiredArg();
		parser_.accepts("runeqtl");		
		parser_.accepts("topeqtl");
		//prioritization
		parser_.accepts("netpath").withRequiredArg();
		
		// Example
		// parser_.accepts("cut").withRequiredArg().ofType(Integer.class).defaultsTo(-1);
	}
	
	
	// ----------------------------------------------------------------------------
    
    /** Get a file / directory saved as a string, throw exception if the name is empty */
	private File getFileOption(String param) {
		
		String filename = (String) options.valueOf(param);
    	if (filename.isEmpty() || filename.equals(" "))
    		throw new RuntimeException("Option '" + param + "': file/directory name is empty or has trailing whitespace");
    	
		return new File(filename.trim());
	}

}
