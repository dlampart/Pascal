/*
Copyright (c) 2013 Daniel Marbach
 
We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper available at:
http://compbio.mit.edu/flynet

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
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Properties;
import java.util.Random;

import org.apache.commons.math3.random.Well19937c;

import cern.jet.random.tdouble.engine.MersenneTwister64;
import ch.unil.genescore.prioritization.ConnectionTreeGetter;
import ch.unil.genescore.prioritization.NetFromSets;
import ch.unil.genescore.prioritization.SparseNet;
import ch.unil.genescore.vegas.AnalyticVegas;
import ch.unil.genescore.vegas.GeneScoreEvaluator;
import ch.unil.genescore.vegas.MaxEffVegas;
import ch.unil.genescore.vegas.MaxSimulAndAnalyticVegas;


/** 
 * Offers global parameters (settings) and functions used by all classes of the
 * Main package.
 * 
 * Code adapted from GnwSettings.java by Thomas Schaffter and Daniel Marbach (gnw.sf.net)
 */
public class Settings {	
	
	/** The configuration file with the settings (leave empty for default settings) */
	static public String settingsFile_ = null;
	/** The properties (settings file) */
	static private Properties set_ = null;
	
	/** Colt Mersenne Twister random engine (should be used by all other random number generators) */
	static public MersenneTwister64 mersenneTwisterRng_ = null;
	/** Apache Commons random engine */
	static public Well19937c wellRng_ = null;
	/** Java random engine */
	static public Random jdkRng_ = null;

	// ----------------------------------------------------------------------------
	// VARIOUS
	
	/** Current version of genescore */
	static public String version_ = "0.1 Alpha";
	/** Seed for the random number generator. Set to -1 to use current time */
	static public int randomSeed_ = -1;
	/** Output directory to save stuff */
	static public String outputDirectory_ = "";
	/** A suffix/ending that is appended to all output files for this run (use to distinguish output files from multiple runs) */
	static public String outputSuffix_ = "";
	/** Set true to display more detailed information on console */
	static public boolean verbose_ = true;

	// INPUT FILES
	/** SNP p-values text file: ID in column 1; p-value in column pvalCol (see below) */ 
	static public String snpPvalFile_ = null;
	/** The column in the tab-separated text file with the p-value (usually column 2) */
	static public int pvalCol_ = -1;

	/** gwas-file has Zscores (rs_nr\tZScore\tPval) */
	static public boolean withZScore_ = false;
	/** fileName containing snp list to which analysis is restricted */
	static public String snpFilterFile_ = null;
	
	/** Directory with the reference population files (must be named 'refPopFilePrefix.chrXX.refPopFileExtension') */
	static public String refPopDirectory_ = null;
	/** Prefix of reference population files (e.g. EUR or some arbitrary name) */
	static public String refPopFilePrefix_ = null;			
	/** Extension of reference population files (defines format, must be 'ser.gz', 'tped.gz' or 'txt.gz') */
	static public String refPopFileExtension_ = null;
	
	/** Set of genes to be considered (leave empty to use all genes from the annotation file) */
	static public String genesToBeLoadedFile_ = null;

	/** The chromosome to be considered (chr1, ..., chr22, chrX, chrY), leave empty for all chromosomes */ 
	static public String chromosome_ = null;
	/** Ignore sex chromosomes */
	static public boolean ignoreAllosomes_ = true;
	/** Ignore mitochondrial genes */
	static public boolean ignoreChrM_ = true;

	/** The genes to be used for gene scoring ('gencode', 'ucsc' or [TODO] a bed file with a custom set of regions) */
	public static String genomeAnnotation_ = null;
	/** The file with the gencode annotation */
	static public String gencodeAnnotationFile_ = null;
	/** UCSC genome browser annotation (use for Entrez IDs) */ 
	static public String ucscAnnotationFile_ = null;
	/**  bed-annotation (use for gene symbols) */ 
	static public String bedAnnotationFile_ = null;
	/** File with Snp-weights for each gene */ 
	static public String geneWiseSnpWeightsFile_ = null;
	/** File format Snp-weights*/ 
	static public String weightFileFormat_ = null;
	//static public String identifyGenesByIdOrSymbol_ = null; 
	/** Set true to load only protein-coding genes */
	static public boolean loadOnlyProteinCodingGenes_ = false;

	/** Mapping file to convert Entrez IDs, ENSEMBL IDs and gene symbols */
	static public String geneIdMappingFile_ = null;

	// PARAMETERS
	/** if phased genotype is known, this option will artificially treat the genotype as unphased */
	static public boolean dePhase_ = false;
	/** Window size up-stream of genes */
	static public int geneWindowUpstream_ = -1;
	/** Window size down-stream of genes */
	static public int geneWindowDownstream_ = -1;
	/** Use the X most significant snps for the test statistic (1: only the most significant snp; -1: all snps) */
	static public int testStatisticNumSnps_ = -1;
		
	/** Max number of snps per gene (display warning and ignore genes with more snps) */
	static public int maxSnpsPerGene_ = -1;
	/** only use snps that have maf above this value in refererence population*/
	static public double useMafCutoff_ = 0;
	/**pruning of correl matrix in maxvegas approach*/
	static public double maxPruningCutoff_ = 0;
	
	/** Analytic solution to Vegas */
	static public boolean useAnalyticVegas_ = false;
	static public boolean useMaxVegas_ = false;
	static public boolean useMaxEffVegas_ = false;
	/** Original version of Vegas with Monte Carlo simulation */
	static public boolean useSimulationVegas_ = false;
	
	/**
	
//	/** all snps onto which we project have to have maf above this cutoff */

	static public double eigenValueFractionCut_ = 1e4; 
	/** Use only GWAS snps */
	static public boolean useOnlyGwasSnps_ = true;

	
	/** The requested absolute precision (used for Imhof and Farebrother methods) */
	static public double requestedAbsolutePrecision_ = -1;
	/** The requested relative precision (used for Imhof and Farebrother methods) */
	static public double requestedRelativePrecision_ = -1;
	/** Maximum number of subintervals in the partition of the given integration interval */
	static public int gslIntegrationLimit_ = -1;
	
	/** Maxiumum number of iterations for Farebrother method */
	static public int farebrotherMaxIterations_ = -1;
	/** Farebrother mode parameter */
	static public double farebrotherMode_ = 1.0;
	/** use Farebrother algorithm*/
	static public boolean useFarebrother_ = false;
	/** use Imhof algorithm*/
	static public boolean useImhof_ = false;
	/** use Davies algorithm*/
	static public boolean useDavies_ = false;
	/** Error bound for Davies method */
	static public double daviesErrorBound_ = -1;
	/** Number of integration terms for Davies method */
	static public int daviesIntegrationTerms_ = -1;

	/** Keep Imhof result if error is below this bound even if the requested precision could not be reached */  
	static public double toleratedAbsolutePrecision_ = -1;
	/** Keep Imhof result if error is below this bound even if the requested precision could not be reached */  
	static public double toleratedRelativePrecision_ = -1;

	/** The delta parameter for weighting SNPs (set to 0 for unweighted, can be a comma-separated list of values) */
	static public ArrayList<Double> snpWeightingDelta_ = null;

	/** The number of samples at each stage when computing empirical p-values */
	static public ArrayList<Integer> adaptiveNumSamples_ = null;
	/** Adaptive estimation of empirical p-values -- require at least this many samples beyond the observed test statistic */
	static public int numSamplesGreaterCutoff_ = -1;

	/** When looking at gene A, filter out coding SNPs of genes B */
	static public boolean removeCodingSnpsOfOtherGenes_ = false;
	/** File where coding snps are defined */
	static public String codingSnpsFile_ = null;

	// OUTPUT
	/** Write more detailed output to the console and the gene score file (estimated errors, Imhof and Farebrother status, ...) */
	static public boolean writeDetailedOutput_ = true; 
	/** Write more detailed output to the console and the gene score file (estimated errors, Imhof and Farebrother status, ...) */
	static public boolean writeDetailedErrorOutput_ = true; 
	/** Write a BED file with the coordinates for all considered SNPs (intersection of study SNPs and reference population SNPs) */
	static public boolean writeSnpBedFile_ = false;
	/** Write plink-Tped file for the genotype that was used */
	static public boolean writeTpedFile_ = false;
	/** Flag to output a separate file for each gene with snp positions, p-values and correlation matrix */
	static public boolean saveGeneReports_ = false;
	/** Path to were we would write all Correlation matrices (if null no writing) */
	static public String writeCorFiles_ = null;
	/** Path to were we would write all Correlation matrices (if null no writing) */
	static public String writeGenewiseSnpFiles_ = null;
	/** Path to were a settingsfile that was used should be written if empty no writing happens */
	static public String writeUsedSettings_ = null;
	
	// PATHWAY ANALYSIS
	/** Set true to run pathway analysis */
	static public boolean runPathwayAnalysis_ = false;
	/** The gene score file: ID in column 1; p-value in column pvalCol (see above) */ 
	//static public String geneScoreFile_ = null;
	/** The gene set library file (.gmt format used by GSEA) */
	static public String geneSetFile_ = null;
	/** Genes to be excluded from enrichment analysis (e.g., MHC region) */
	static public String excludedGenesFile_ = null;
	/** Create meta-genes by merging genes that are closer than the given distance (given in *megabases*; set to -1 to disable meta-genes) */
	static public double mergeGenesDistance_ = -1;
	/** Use empirical distribution of p-values by sampling over scores for gene set enrichment*/
	static public boolean useSimulation_ = false;
	/** Use empirical distribution of p-values by sampling over scores for gene set enrichment*/
	static public boolean useSimulationWeightedSampling_ = false;
	
	/** max nr of simulation runs.  at least 100'000 is required. */
	static public int maxNrOfSimulationsForEnrichment_ = 100000;
	/** Use chi-squared distribution to compute p-values for gene set enrichment */
	static public boolean useChi2_ = false; 
	/** Use rank-sum-test to compute p-values for gene set enrichment */
	static public boolean useRankSum_ = false; 
	/** Use hypergeometric distribution to compute p-values for gene set enrichment */
	static public boolean useHypGeom_ = false; 
	/** quantile to use for hypergeometric distribution (lower tail)*/
	static public double[] hypGeomQuantiles_ = {-1}; 
	
	static public boolean useGamma_ = false; 
	/** shape parameter to use for gamma distribution*/
	static public double[] gammaShapeParameters_ = {0.1,0.2,0.3,0.5,1,2,4}; 
	
	static public boolean useExpHyp_ = true; 
	/** shape parameter to use for gamma distribution (lower tail)*/
	static public double[] expHypParameters_ = {0.5,0.75,0.85,0.9};
	
	static public boolean loadScoresFromFiles_ = false; 
	static public String geneScoreFile_ = null; 
	static public String metaGeneScoreFile_ = null; 
	
	/** Write only gene sets that pass the given significance level */
	static public double writeSignificanceThreshold_ = -1;
	
	/** Set true to write gene and meta-gene scores to a file */ 
	static public boolean writeGeneScores_ = false;

	// EQTL-Analysis
	/** Set true to run eqtl analysis */
	static public boolean runEqtlAnalysis_ = false;
	static public String eqtlFile_ = null;
	public static boolean onlyTopOverlappedEqtl_ = false;
	public static boolean runEqtlProjection_ = true
			;
	
	
	//TopLdSnp
	static public boolean runTopLdSnp_ = false;
	static public String firstList_ = "resources/eqtl/BrainDownloads/snpList.txt";
	static public String secondList_ = "resources/gwas/dbgap/EUR.IBDGenetics.CD.txt";
	static public String topLdSnpOutFile_ = "resources/eqtl/BrainDownloads/snpListLdPartners.txt";
	
	
	//prioritizationAnalysis
	static public boolean runPrioritizationAnalysis_ = true;
	static public boolean useSparseNet_ = true;
	static public boolean useGeneSets_ = true;
	static public String netPath_ = "";
	
	// CONCATENATE CHROMOSOME RESULT FILES
	/** Set true to run concatenate chromosome results */
	static public boolean runConcatenateChromosomeResults_ = false;
	/** Concatenate individual chromosome result files in the given directory */
	static public String concatenateChromosomeResultsDir_ = null;	
	/** Delete original chromosome result files after concatenating them */
	static public boolean deleteOriginals_ = false;
	public static boolean useFakePhenotype_;
	public static boolean useFakeSignal_;
	public static double chanceOfSignal_=0.0;
	public static int multipleOfPhenotype_=0;
	

	private static Object randSeed_;
	public static double deflationRate_;
	public static int deflationDist_;
	public static Boolean onlyPathwayGenesAsBackground_;	
	
	// processed fields
	
	public static String gwasName_;	
	public static String chromFileExtension_;	
	
	private static  GeneScoreEvaluator GeneScoreEvaluatorInstance_ = null;
	private static  ConnectionTreeGetter  ConnectionGetterInstance_ = null;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Initialize settings */
	static public void initialize() {
		
		initializeRandomNumberGenerators();
	}
	
	// ----------------------------------------------------------------------------
	public void firstProcessingOfSettings(){
		gwasName_ = Utils.extractBasicFilename(Settings.snpPvalFile_, false);
		chromFileExtension_ = Settings.chromosome_.equals("") ? "" : "." + Settings.chromosome_;
		setGenescoreEvaluatorInstance();		
		setConnectionGetterInstance();
	}
	
	private void setConnectionGetterInstance(){
		if (Settings.useSparseNet_){
			SparseNet mySparseNet = new SparseNet();
			mySparseNet.setParser(netPath_);
			ConnectionGetterInstance_  = mySparseNet;
		}
			else if (Settings.useGeneSets_){
				netPath_ = geneSetFile_;
			ConnectionGetterInstance_ = new NetFromSets();
			
		}
	}
	
	private void setGenescoreEvaluatorInstance(){
		
		GeneScoreEvaluator evaluator = null;
		
		
		// Determine which gene scoring method to use
		if (Settings.useAnalyticVegas_) {
			
			evaluator = new AnalyticVegas();
			
		
		} else if (Settings.useMaxVegas_) {
		
			evaluator = new MaxSimulAndAnalyticVegas();
			
		} else if (Settings.useMaxEffVegas_) {
			
					evaluator = new MaxEffVegas();
		}
        else {
			throw new RuntimeException("No gene scoring method selected");
		}
		GeneScoreEvaluatorInstance_ = evaluator;
	}
	
	/** Load and initialize settings */
	static public void loadSettings() {
		
		try {
			InputStream in;
			if (settingsFile_ == null || settingsFile_.compareTo("") == 0) {
				Main.println("- No settings file specified");
				Main.println("- Loading default settings\n");
			//	in = Settings.class.getClassLoader().getResourceAsStream("ch/unil/genescore/main/settings.txt");				
				in = new FileInputStream("settings.txt");
			} else {
				Main.println("- Loading settings file: " + settingsFile_ + "\n");
				in = new FileInputStream(settingsFile_);
			}
			set_ = new Properties();
			set_.load(new InputStreamReader(in));
			
			setParameterValues();
			checkParamBounds();
			
		} catch (Exception e) {
			Main.warning(e.getMessage());
			String msg = "Failed to load settings file (a parameter may be missing or malformed)";
			if (settingsFile_ != null)
				msg += ": " + settingsFile_; 
			Main.error(msg);
		}
		
		// Reinitialize the random number generators
		initializeRandomNumberGenerators();		
		checkOptions();
	//	warnOptions();
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get folder of the default Main directory (in home) */
	public File getGseaDirectory() {
		
		File folder = new File(MainDirectoryPath());
		return folder;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get file associated to default settings file */
	public File getCustomMainSettings() {
		
		File file = new File(personalMainSettingsPath());
		return file;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get path to default Main directory (in home) */
	public String MainDirectoryPath() {
		
		return System.getProperty("user.home")
				+ System.getProperty("file.separator")
				+ "Main";
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get path to default settings file */
	public String personalMainSettingsPath() {
		
		return MainDirectoryPath()
				+ System.getProperty("file.separator")
				+ "settings.txt";
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Return true if a custom default Main settings file exists */
	public boolean personalMainSettingsExist() {
		
		return (new File(personalMainSettingsPath())).exists();
	}
	
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Set the user path. Could be with or without "/" terminal.
	 * @param absPath Absolute path
	 */
	public void setOutputDirectory(String absPath) {
		
		outputDirectory_ = absPath;
		String sep = System.getProperty("file.separator");
		if (outputDirectory_.charAt(outputDirectory_.length()-1) != sep.charAt(0))
			outputDirectory_ += sep;
	}
	
	
	// ============================================================================
	// PRIVATE METHODS

	/** Create new instances for the random number generators, initialize with randomSeed_ */
	protected static void initializeRandomNumberGenerators() {
		System.out.println("random seed:");
		System.out.println(Integer.toString(randomSeed_));
		if (randomSeed_ == -1) {
			mersenneTwisterRng_ = new MersenneTwister64(new java.util.Date());
			wellRng_ = new Well19937c();
			jdkRng_ = new Random();
		} else {
			mersenneTwisterRng_ = new MersenneTwister64(randomSeed_);
			wellRng_ = new Well19937c(randomSeed_);
			jdkRng_ = new Random(randomSeed_);
		}
		
		//uniformDistribution_ = new Uniform(mersenneTwister_);
		//normalDistribution_ = new Normal(0, 1, mersenneTwister_); // mean 0, stdev 1
	}
	
	
	// ----------------------------------------------------------------------------

	/** Set Main parameters based on the loaded properties */
	static private void setParameterValues() throws Exception {

		// VARIOUS
		randomSeed_ = getSettingInt("randomSeed");
		outputDirectory_ = getSetting("outputDirectory");
		if (outputDirectory_.equals("")) 
			outputDirectory_ = System.getProperty("user.dir");
		outputSuffix_ = getSetting("outputSuffix");
		verbose_ = getSettingBoolean("verbose");
		
		// ----------------------------------------------------------------------------
		// ENRICHMENT ANALYSIS
		
		// INPUT FILES
		snpPvalFile_ = getSetting("snpPvalFile");
		pvalCol_ = getSettingInt("pvalCol");
		
		refPopDirectory_ = getSetting("refPopDirectory");
		refPopFilePrefix_ = getSetting("refPopFilePrefix");
		refPopFileExtension_ =getSetting("refPopFileExtension");
		
		withZScore_ = getSettingBoolean("withZScore");
		snpFilterFile_ = getSetting("snpFilterFile");
		dePhase_ = getSettingBoolean("dePhase");
		
		genesToBeLoadedFile_ = getSetting("genesToBeLoadedFile");
		
		writeDetailedOutput_ = getSettingBoolean("writeDetailedOutput");
		writeDetailedErrorOutput_ = getSettingBoolean("writeDetailedErrorOutput");
		writeSnpBedFile_ = getSettingBoolean("writeSnpBedFile");
		writeTpedFile_ = getSettingBoolean("writeTpedFile");
		
		chromosome_ = getSetting("chromosome");
		ignoreAllosomes_ = getSettingBoolean("ignoreAllosomes");
		ignoreChrM_ = getSettingBoolean("ignoreChrM");
		
		geneWindowUpstream_ = getSettingInt("geneWindowUpstream");
		geneWindowDownstream_ = getSettingInt("geneWindowDownstream");
		testStatisticNumSnps_ = getSettingInt("testStatisticNumSnps");
				
		maxSnpsPerGene_ = getSettingInt("maxSnpsPerGene");
		useMafCutoff_ = getSettingDouble("useMafCutoff");
		
		useAnalyticVegas_ = getSettingBoolean("useAnalyticVegas");
		useMaxVegas_ = getSettingBoolean("useMaxVegas");
		useMaxEffVegas_ = getSettingBoolean("useMaxEffVegas");
		useSimulationVegas_ = getSettingBoolean("useSimultationVegas");
		
		maxPruningCutoff_=getSettingDouble("maxPruningCutoff");
		useFakePhenotype_=getSettingBoolean("useFakePhenotype");
		useFakeSignal_=getSettingBoolean("useFakeSignal");
		chanceOfSignal_=getSettingDouble("chanceOfSignal");
		multipleOfPhenotype_=getSettingInt("multipleOfPhenotype");
		
		

		eigenValueFractionCut_ = getSettingDouble("eigenValueFractionCut");

		requestedAbsolutePrecision_ = getSettingDouble("requestedAbsolutePrecision");
		requestedRelativePrecision_ = getSettingDouble("requestedRelativePrecision");
		gslIntegrationLimit_ = getSettingInt("gslIntegrationLimit");
		
		farebrotherMaxIterations_ = getSettingInt("farebrotherMaxIterations");
		farebrotherMode_ = getSettingDouble("farebrotherMode");
		
		useImhof_ = getSettingBoolean("useImhof");		
		useFarebrother_ = getSettingBoolean("useFarebrother");		
		useDavies_ = getSettingBoolean("useDavies");		
		daviesErrorBound_ = getSettingDouble("daviesErrorBound");
		daviesIntegrationTerms_ = getSettingInt("daviesIntegrationTerms");
		
		toleratedAbsolutePrecision_ = getSettingDouble("toleratedAbsolutePrecision");
		toleratedRelativePrecision_ = getSettingDouble("toleratedRelativePrecision");
		
		genomeAnnotation_ = getSetting("genomeAnnotation");
		gencodeAnnotationFile_ = getSetting("genecodeAnnotationFile");
		ucscAnnotationFile_ = getSetting("ucscAnnotationFile");
		bedAnnotationFile_ = getSetting("bedAnnotationFile");
		geneWiseSnpWeightsFile_ = getSetting("geneWiseSnpWeightsFile");
		weightFileFormat_= getSetting("weightFileFormat");
		//identifyGenesByIdOrSymbol_ = getSetting("identifyGenesByIdOrSymbol");
		loadOnlyProteinCodingGenes_ = getSettingBoolean("loadOnlyProteinCodingGenes");
		
		geneIdMappingFile_ = getSetting("geneIdMappingFile");
		
		snpWeightingDelta_ = getSettingDoubleArray("snpWeightingDelta", false);
		
		adaptiveNumSamples_ = getSettingIntArray("adaptiveNumSamples", true);
		numSamplesGreaterCutoff_ = getSettingInt("numSamplesGreaterCutoff");
		
		removeCodingSnpsOfOtherGenes_ = getSettingBoolean("removeCodingSnpsOfOtherGenes");
		codingSnpsFile_ = getSetting("codingSnpsFile");
		
		saveGeneReports_ = getSettingBoolean("saveGeneReports");
		writeCorFiles_ = getSetting("corFilesPath");
		writeGenewiseSnpFiles_=getSetting("genewiseSnpFilesPath");
		runPathwayAnalysis_ = getSettingBoolean("runPathwayAnalysis");
		deflationRate_ = getSettingDouble("deflationRate");
		deflationDist_ = getSettingInt("deflationDistance");
		//geneScoreFile_ = getSetting("geneScoreFile");
		geneSetFile_ = getSetting("geneSetFile");
		excludedGenesFile_ = getSetting("excludedGenesFile");
		mergeGenesDistance_ = getSettingDouble("mergeGenesDistance");
		useSimulation_ = getSettingBoolean("useSimulation");
		useSimulationWeightedSampling_ = getSettingBoolean("useSimulationWeightedSampling");
		maxNrOfSimulationsForEnrichment_ = getSettingInt("maxNrOfSimulationsForEnrichment");
		useChi2_ = getSettingBoolean("useChi2");
		useRankSum_ = getSettingBoolean("useRankSum");
		useHypGeom_ = getSettingBoolean("useHypGeom");
		hypGeomQuantiles_ = getSettingDoubleAr("hypGeomQuantiles");
		useGamma_ = getSettingBoolean("useGamma");
		gammaShapeParameters_ = getSettingDoubleAr("gammaShapeParameters");
		useExpHyp_ = getSettingBoolean("useExpHyp");
		expHypParameters_ = getSettingDoubleAr("expHypParameters");
		
		loadScoresFromFiles_ = getSettingBoolean("loadScoresFromFiles");
		geneScoreFile_ = getSetting("geneScoreFile");
		metaGeneScoreFile_ = getSetting("metaGeneScoreFile");
		onlyPathwayGenesAsBackground_ = getSettingBoolean("onlyPathwayGenesAsBackground");
		
		
		writeSignificanceThreshold_ = getSettingDouble("writeSignificanceThreshold");
		writeGeneScores_ = getSettingBoolean("writeGeneScores");
		writeUsedSettings_ = getSetting("writeUsedSettings");
		
		eqtlFile_ = getSetting("eqtlFile");
		runEqtlAnalysis_ = getSettingBoolean("runEqtlAnalysis");
		onlyTopOverlappedEqtl_ = getSettingBoolean("onlyTopOverlappedEqtl");
		
		//prioritizationAnalysis
		runPrioritizationAnalysis_ = getSettingBoolean("runPrioritizationAnalysis");
		useSparseNet_ = getSettingBoolean("useSparseNet");
		useGeneSets_ = getSettingBoolean("useGeneSets");
		netPath_ = getSetting("netPath");
		
		runConcatenateChromosomeResults_ = getSettingBoolean("runConcatenateChromosomeResults");
		concatenateChromosomeResultsDir_ = getSetting("concatenateChromosomeResultsDir");
		deleteOriginals_ = getSettingBoolean("deleteOriginals");
	}
	
	
	// ----------------------------------------------------------------------------

	/** Check if the parameters are valid */
	static private void checkParamBounds() {
				
		if (runPathwayAnalysis_) {
			// At least one enrichment evaluation method must be selected
			if (!useSimulation_ && !useChi2_ && !useRankSum_)
				throw new IllegalArgumentException("Set either 'useRankSum' and/or 'useChi2' and/or useSimulation to true");
		}
		checkBounds(writeSignificanceThreshold_, 0, 1, "writeSignificanceThreshold");

	}

	
	// ----------------------------------------------------------------------------

	/** Check if the parameters are valid */
	static private void checkBounds(double param, double min, double max, String paramName) {

		if (param < min || param > max)
			throw new IllegalArgumentException("Parameter " + paramName + "=" + param + " is outside of its valid range [" + min + ", " + max + "]");
	}

		
	// ----------------------------------------------------------------------------

	/** Get the string value of a parameter from the setting file */
	static private String getSetting(String param) {
		
		String value = set_.getProperty(param);
		if (value == null)
			Main.error("Parameter not found in setting file: " + param);
		
		return value; 
	}

	
	// ----------------------------------------------------------------------------

	/** Get the integer value of a parameter from the setting file */
	static private int getSettingInt(String param) {
		return Integer.valueOf(getSetting(param)); 
	}

	/** Get the double value of a parameter from the setting file */
	static private double getSettingDouble(String param) {
		return Double.valueOf(getSetting(param)); 
	}

	// ----------------------------------------------------------------------------

	/** Parse a boolean property */
	static private boolean getSettingBoolean(String name) {
		
		String value = getSetting(name);
		if (value.equals("1") || value.equalsIgnoreCase("true") || value.equalsIgnoreCase("t"))
			return true;
		else if (value.equals("0") || value.equalsIgnoreCase("false") || value.equalsIgnoreCase("f"))
			return false;
		else
			throw new IllegalArgumentException("Parse error for boolean parameter '" + name + "': expected '1' or '0', found '" + value + "'");
	}

	/** Parse a double array comma separated*/
	protected static double[] getSettingDoubleAr(String param) {
		
		String value = getSetting(param);
		String[] strAr= value.split(",");
		double[] myAr = new double[strAr.length];
		for (int i=0;i<strAr.length;i++){
			myAr[i]=Double.valueOf(strAr[i]);
		}
		return(myAr);
	}

	// ----------------------------------------------------------------------------

	/** Parse an int array property */
	static private ArrayList<Integer> getSettingIntArray(String name, boolean positiveSorted) {
		
		String[] propStr = getSetting(name).split(",");
		ArrayList<Integer> prop = new ArrayList<Integer>();
		
		if (propStr.length == 1 && propStr[0].compareTo("") == 0)
			return prop;

		for (int i=0; i<propStr.length; i++)
			prop.add(Integer.valueOf(propStr[i]));
			
		if (positiveSorted && !Utils.posIntIncreasing(prop))
			Main.error("Error parsing settings file, " + name + " has to be an ordered list of positive integers, given in increasing order");
		
		return prop;
	}


	// ----------------------------------------------------------------------------

	/** Parse an double array property */
	static private ArrayList<Double> getSettingDoubleArray(String name, boolean positiveSorted) {
		
		String[] propStr = getSetting(name).split(",");
		ArrayList<Double> prop = new ArrayList<Double>();
		
		if (propStr.length == 1 && propStr[0].compareTo("") == 0)
			return prop;

		for (int i=0; i<propStr.length; i++)
			prop.add(Double.valueOf(propStr[i]));
			
		if (positiveSorted && !Utils.posDoubleIncreasing(prop))
			Main.error("Error parsing settings file, " + name + " has to be an ordered list of positive numbers, given in increasing order");
		
		return prop;
	}

	
	// ----------------------------------------------------------------------------
	
	
	/** Parse a string array property */
	@SuppressWarnings("unused")
	static private String[] getStringArraySetting(Properties set, String name) {
		
		return set.getProperty(name).split(",");
	}
	public static void checkOptions(){

		int int1 = (useAnalyticVegas_) ? 1 : 0;
		int int2 = (useMaxVegas_) ? 1 : 0;
		int int3 = (useMaxEffVegas_) ? 1 : 0;
		int int4 = (useSimulationVegas_) ? 1 : 0;
		int tot = int1 + int2 + int3 + int4;
		if (tot>1){
			throw new RuntimeException("error: more than 1 vegas option active");
		}
		
		int1 = (useImhof_) ? 1 : 0;
		int2 = (useDavies_) ? 1 : 0;
		int3 = (useFarebrother_) ? 1 : 0;		
		tot = int1 + int2 + int3;
		if (tot>1){
			throw new RuntimeException("error: more than 1 algorithm option active");
		}
	}
public static void warnOptions(){
		
		if (useOnlyGwasSnps_){
			System.out.println("WARNING: <useonlygwassnps> option selected in file:");
			try {
				Thread.sleep(2400);
			} catch(InterruptedException ex) {
				Thread.currentThread().interrupt();
			}
		}
	}
public static GeneScoreEvaluator getGeneScoreEvaluator(){
	return GeneScoreEvaluatorInstance_;
}
public static ConnectionTreeGetter getConnectionGetter(){
	return ConnectionGetterInstance_;
}
}
