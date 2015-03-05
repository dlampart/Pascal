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
package ch.unil.genescore.main;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

/**
 * Defines and parses the command-line arguments
 */
public class GeneScoreOptionParser extends Settings {

	/** The option parser */
	OptionParser parser_ = null;

	// ============================================================================
	// PUBLIC METHODS

	/** Parse the command-line arguments, initialize all settings */
	public GeneScoreOptionParser() {

		// Defines the arguments
		defineArgs();
	}

	// ----------------------------------------------------------------------------

	/** Parses the command-line arguments, which were defined by defineArgs() */
	public void parse(String[] args) {

		// (1) Parse the options
		OptionSet options = parser_.parse(args);

		// Display help
		if (options.has("help")) {
			displayHelp();
			System.exit(0);
		}

		// (2) Set the settings file
		if (options.has("set"))
			settingsFile_ = (String) options.valueOf("set");

		// (3) Load the settings file
		loadSettings();

		// (4) Set command-line options (default is what has been loaded from
		// settings file)
		if (options.has("outdir"))
			outputDirectory_ = (String) options.valueOf("outdir");
		if (options.has("outsuffix"))
			outputSuffix_ = (String) options.valueOf("outsuffix");		
		if (options.has("writesnpbed"))
			writeSnpBedFile_ = true;
		if (options.has("ucscdir"))
			ucscAnnotationFile_ = (String) options.valueOf("ucscdir");
		if (options.has("refpopdir"))
			refPopDirectory_ = (String) options.valueOf("refpopdir");
		if (options.has("genestobeloaded"))
			genesToBeLoadedFile_ = (String) options.valueOf("genestobeloaded");
		if (options.has("pval"))
			snpPvalFile_ = (String) options.valueOf("pval");
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
		
		if (options.has("nsnp"))
			testStatisticNumSnps_ = (Integer) options.valueOf("nsnp");
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
		
		if (options.has("randseed")){
			randomSeed_ = (Integer) options.valueOf("randseed");
			initializeRandomNumberGenerators();
		}
		
		
		if (options.has("useimhof"))
			useImhof_=true;		
		if (options.has("usedavies"))
			useDavies_=true;
		if (options.has("usefarebrother"))
			useFarebrother_=true;
		
		if (options.has("genewisesnpdir"))
			writeGenewiseSnpFiles_= (String) options.valueOf("genewisesnpdir");
		
		if (options.has("analyticvegas"))
			useAnalyticVegas_ = true;
		if (options.has("maxvegas"))
			useMaxVegas_ = true;
		if (options.has("maxeffvegas"))
		useMaxEffVegas_ = true;

		if (options.has("withzscore"))
			withZScore_ = true;
		if (options.has("prunelist"))
			snpFilterFile_ = (String) options.valueOf("filterlist");
		if (options.has("up"))
			geneWindowUpstream_ = (Integer) options.valueOf("up");
		if (options.has("down"))
			geneWindowDownstream_ = (Integer) options.valueOf("down");
		//Pathway options
		if (options.has("annotation"))
			genomeAnnotation_ = (String) options.valueOf("annotation");
		if (options.has("gencodefile"))
			gencodeAnnotationFile_ = (String) options.valueOf("gencodefile");
		
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
			geneSetFile_ = (String) options.valueOf("genesetfile");
		
		if (options.has("hypgeom"))
			useHypGeom_ = true;
		
		if (options.has("simulated"))
			useSimulation_ = true;
		
		if (options.has("loadscores"))
			loadScoresFromFiles_= true;
		
		if (options.has("scores"))
			geneScoreFile_ = (String) options.valueOf("scores");
		
		if (options.has("metascores"))
			metaGeneScoreFile_ = (String) options.valueOf("metascores");
		if (options.has("deflationrate"))
			deflationRate_ = (Double) options.valueOf("deflationrate");
		
		//path to pathway
		if (options.has("detailed"))
			writeDetailedOutput_ = true;
		if (options.has("rewritesettings"))
		writeUsedSettings_ = (String) options.valueOf("scores");
		//path to pathway

	

		// TODO, write a method that checks consistency / if everything has been defined that we need
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
	
	

}
