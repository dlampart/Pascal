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
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Properties;
import java.util.Random;

import org.apache.commons.math3.random.Well19937c;


import ch.unil.gpsutils.Settings;


/** 
 * Offers global parameters (settings) and functions used by all classes of the
 * Main package.
 */
public class PascalSettings extends Settings {	
	
	/** Current version of genescore */
	public final String version_ = "1.0 alpha";

	/** The annotation file (default file included in jar) */
	public final String annotationRsc = "edu/mit/magnum/gene/rsc/gene_coord.bed";
	
	/** Apache Commons random engine */
	public Well19937c wellRng_;
	/** Java random engine */
	public Random jdkRng_;


	// ----------------------------------------------------------------------------
	// VARIOUS
	
	/** Output directory to save stuff */
	public File outputDirectory_;
	/** A suffix/ending that is appended to all output files for this run (use to distinguish output files from multiple runs) */
	public String outputSuffix_;
	/** Set true to display more detailed information on console */
	public boolean verbose_;
	/** PRIVATE, NEEDS TO BE SET WITH setRandomSeed(), which initializes the random number generators. Set to -1 to use current time */
	private int randomSeed_;

	// ----------------------------------------------------------------------------
	// INPUT
	
	/** SNP p-values text file: ID in column 1; p-value in column pvalCol (see below) */ 
	public File snpPvalFile_;
	/** The column in the tab-separated text file with the p-value (default column 2) */
	public int pvalCol_;
	/** GWAS file has Zscores (rs_nr\tZScore\tPval) */
	public boolean withZScore_; // TODO are z-scores used? 
	
	/** Directory with the reference population files (must be named 'refPopFilePrefix.chrXX.refPopFileExtension') */
	public File refPopDirectory_;
	/** Prefix of reference population files (e.g. EUR or some arbitrary name) */
	public String refPopFilePrefix_; // TODO delete, extract from file names
	/** Extension of reference population files (defines format, must be 'ser.gz', 'tped.gz' or 'txt.gz') */
	public String refPopFileExtension_; // TODO delete?
	
	/** File containing snp list to which analysis is restricted */
	public File snpFilterFile_;
	/** Set of genes to be considered (leave empty to use all genes from the annotation file) */
	public File genesToBeLoadedFile_;

	/** The chromosome to be considered (chr1, ..., chr22, chrX, chrY), leave empty for all chromosomes */ 
	public String chromosome_;
	/** Ignore sex chromosomes */
	public boolean ignoreAllosomes_;
	/** Ignore mitochondrial genes */
	public boolean ignoreChrM_;

	/** The genes to be used for gene scoring ('gencode', 'ucsc' or [TODO] a bed file with a custom set of regions) */
	public String genomeAnnotation_;
	/** The file with the gencode annotation */
	public File gencodeAnnotationFile_;
	/** UCSC genome browser annotation (use for Entrez IDs) */ 
	public File ucscAnnotationFile_;
	/**  bed-annotation (use for gene symbols) */ 
	public File bedAnnotationFile_;
	/** Set true to load only protein-coding genes */
	public boolean loadOnlyProteinCodingGenes_;
	/** Mapping file to convert Entrez IDs, ENSEMBL IDs and gene symbols */
	public File geneIdMappingFile_;
	
	/** File with Snp-weights for each gene */ 
	public File geneWiseSnpWeightsFile_;
	/** File format Snp-weights*/ 
	public String weightFileFormat_;

	// ----------------------------------------------------------------------------
	// PARAMETERS

	/** Window size up-stream of genes */
	public int geneWindowUpstream_;
	/** Window size down-stream of genes */
	public int geneWindowDownstream_;
		
	/** Max number of snps per gene (display warning and ignore genes with more snps) */
	public int maxSnpsPerGene_;
	/** Only use snps that have maf above this value in refererence population */
	public double useMafCutoff_;
	/** if phased genotype is known, this option will artificially treat the genotype as unphased */
	public boolean dePhase_;
	/** Check converge of Max vegas externally */
	public boolean externalConvergenceCheck_;
	/** Matrix in maxvegas approach*/
	public double maxPruningCutoff_;
	
	/** Analytic solution to Vegas */
	public boolean useAnalyticVegas_;
	/** TODO describe */
	public boolean useMaxVegas_;
	/** TODO describe */
	public boolean useMaxEffVegas_;
	/** Original version of Vegas with Monte Carlo simulation */
	public boolean useSimulationVegas_;
	
	/** all snps onto which we project have to have maf above this cutoff */
	public double eigenValueFractionCut_; 
	/** Use only GWAS snps */
	public boolean useOnlyGwasSnps_;

	/** The requested absolute precision (used for Imhof and Farebrother methods) */
	public double requestedAbsolutePrecision_;
	/** The requested relative precision (used for Imhof and Farebrother methods) */
	public double requestedRelativePrecision_;
	/** Maximum number of subintervals in the partition of the given integration interval */
	public int gslIntegrationLimit_;
	
	/** Maxiumum number of iterations for Farebrother method */
	public int farebrotherMaxIterations_;
	/** Farebrother mode parameter */
	public double farebrotherMode_;
	/** use Farebrother algorithm*/
	public boolean useFarebrother_;
	/** use Imhof algorithm*/
	public boolean useImhof_;
	/** use Davies algorithm*/
	public boolean useDavies_;
	/** Error bound for Davies method */
	public double daviesErrorBound_;
	/** Number of integration terms for Davies method */
	public int daviesIntegrationTerms_;

	/** Keep Imhof result if error is below this bound even if the requested precision could not be reached */  
	public double toleratedAbsolutePrecision_;
	/** Keep Imhof result if error is below this bound even if the requested precision could not be reached */  
	public double toleratedRelativePrecision_;

	/** The delta parameter for weighting SNPs (set to 0 for unweighted, can be a comma-separated list of values) */
	public ArrayList<Double> snpWeightingDelta_;

	/** The number of samples at each stage when computing empirical p-values */
	public ArrayList<Integer> adaptiveNumSamples_;
	/** Adaptive estimation of empirical p-values -- require at least this many samples beyond the observed test statistic */
	public int numSamplesGreaterCutoff_;

	/** When looking at gene A, filter out coding SNPs of genes B */
	public boolean removeCodingSnpsOfOtherGenes_;
	/** File where coding snps are defined */
	public String codingSnpsFile_;

	// ----------------------------------------------------------------------------
	// OUTPUT

	/** Write more detailed output to the console and the gene score file (estimated errors, Imhof and Farebrother status, ...) */
	public boolean writeDetailedOutput_; 
	/** Write more detailed output to the console and the gene score file (estimated errors, Imhof and Farebrother status, ...) */
	public boolean writeDetailedErrorOutput_; 
	/** Write a BED file with the coordinates for all considered SNPs (intersection of study SNPs and reference population SNPs) */
	public boolean writeSnpBedFile_;
	/** Write plink-Tped file for the genotype that was used */
	public boolean writeTpedFile_;
	/** Flag to output a separate file for each gene with snp positions, p-values and correlation matrix */
	public boolean saveGeneReports_;
	/** Path to were we would write all Correlation matrices (if null no writing) */
	public boolean writeCorFiles_; // TODO delete if redundant
	/** Path to were we would write all Correlation matrices (if null no writing) */
	public boolean writeGenewiseSnpFiles_; // TODO delete if redundant with saveGeneReports
	/** Path to were a settingsfile that was used should be written if empty no writing happens */
	public File writeUsedSettings_; // TODO delete, check that it is always written
	
	// ----------------------------------------------------------------------------
	// PATHWAY ANALYSIS

	/** Set true to run pathway analysis */
	public boolean runPathwayAnalysis_;
	/** The gene set library file (.gmt format used by GSEA) */
	public File geneSetFile_;
	/** Genes to be excluded from enrichment analysis (e.g., MHC region) */
	public File excludedGenesFile_;
	/** Create meta-genes by merging genes that are closer than the given distance (given in *megabases*; set to -1 to disable meta-genes) */
	public double mergeGenesDistance_;
	
	/** Use empirical distribution of p-values by sampling over scores for gene set enrichment*/
	public boolean useSimulation_;
	/** Use empirical distribution of p-values by sampling over scores for gene set enrichment */
	public boolean useSimulationWeightedSampling_; // TODO difference from above
	/** max nr of simulation runs.  at least 100'000 is required. */
	public int maxNrOfSimulationsForEnrichment_;
	/** Use chi-squared distribution to compute p-values for gene set enrichment */
	public boolean useChi2_; 
	/** Use rank-sum-test to compute p-values for gene set enrichment */
	public boolean useRankSum_; 
	/** Use hypergeometric distribution to compute p-values for gene set enrichment */
	public boolean useHypGeom_; 
	/** quantile to use for hypergeometric distribution (lower tail) */
	public ArrayList<Double> hypGeomQuantiles_; 
	/** TODO describe */
	public boolean useGamma_; 
	/** shape parameter to use for gamma distribution */
	public ArrayList<Double> gammaShapeParameters_; 
	/** TODO describe */
	public boolean useExpHyp_; 
	/** shape parameter to use for gamma distribution (lower tail)*/
	public ArrayList<Double> expHypParameters_;
	
	/** Load pre-computed gene scores */
	public boolean loadScoresFromFiles_;
	/** File with pre-computed gene scores */
	public File geneScoreFile_; 
	/** TODO What is this? */
	public File metaGeneScoreFile_; 
	
	/** Write only gene sets that pass the given significance level */
	public double writeSignificanceThreshold_;
	/** Set true to write gene and meta-gene scores to a file */ 
	public boolean writeGeneScores_;

	
	// CONCATENATE CHROMOSOME RESULT FILES
	/** Set true to run concatenate chromosome results */
	public boolean runConcatenateChromosomeResults_;
	/** Concatenate individual chromosome result files in the given directory */
	public String concatenateChromosomeResultsDir_;	
	/** Delete original chromosome result files after concatenating them */
	public boolean deleteOriginals_;
	
	// TODO describe
	public boolean useFakePhenotype_;
	public boolean useFakeSignal_;
	public double chanceOfSignal_;
	public int multipleOfPhenotype_;

	public double deflationRate_;
	public int deflationDist_;
	public Boolean onlyPathwayGenesAsBackground_;	
	
	// processed fields // TODO move to respective classes?
	public String gwasName_;	
	public String chromFileExtension_;	
	
	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor */
	public PascalSettings() {
		resetToDefaults();
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Set default values for all settings */
	public void resetToDefaults() {

		// VARIOUS
		outputDirectory_ = new File(System.getProperty("user.dir"));
		outputSuffix_ = "";
		verbose_ = true;
		// Initializes the RNGs
		setRandomSeed(42);

		// INPUT FILES
		snpPvalFile_ = null;
		pvalCol_ = 2;
		withZScore_ = false;

		refPopDirectory_ = null;
		refPopFilePrefix_ = null;			
		refPopFileExtension_ = null;

		snpFilterFile_ = null;
		genesToBeLoadedFile_ = null;

		chromosome_ = null;
		ignoreAllosomes_ = true;
		ignoreChrM_ = true;

		genomeAnnotation_ = null;
		gencodeAnnotationFile_ = null;
		ucscAnnotationFile_ = null;
		bedAnnotationFile_ = null;
		loadOnlyProteinCodingGenes_ = false;
		geneIdMappingFile_ = null;

		geneWiseSnpWeightsFile_ = null;
		weightFileFormat_ = null;

		// PARAMETERS
		geneWindowUpstream_ = -1;
		geneWindowDownstream_ = -1;

		maxSnpsPerGene_ = -1;
		useMafCutoff_ = 0;
		dePhase_ = false;
		externalConvergenceCheck_ = false;
		maxPruningCutoff_ = 0;

		useAnalyticVegas_ = false;
		useMaxVegas_ = false;
		useMaxEffVegas_ = false;
		useSimulationVegas_ = false;

		eigenValueFractionCut_ = 1e4; 
		useOnlyGwasSnps_ = true;

		requestedAbsolutePrecision_ = -1;
		requestedRelativePrecision_ = -1;
		gslIntegrationLimit_ = -1;

		farebrotherMaxIterations_ = -1;
		farebrotherMode_ = 1.0;
		useFarebrother_ = false;
		useImhof_ = false;
		useDavies_ = false;
		daviesErrorBound_ = -1;
		daviesIntegrationTerms_ = -1;

		toleratedAbsolutePrecision_ = -1;
		toleratedRelativePrecision_ = -1;

		snpWeightingDelta_ = null;

		adaptiveNumSamples_ = null;
		numSamplesGreaterCutoff_ = -1;

		removeCodingSnpsOfOtherGenes_ = false;
		codingSnpsFile_ = null;

		// OUTPUT
		writeDetailedOutput_ = true; 
		writeDetailedErrorOutput_ = true; 
		writeSnpBedFile_ = false;
		writeTpedFile_ = false;
		saveGeneReports_ = false;
		writeCorFiles_ = false;
		writeGenewiseSnpFiles_ = false;
		writeUsedSettings_ = null;

		// PATHWAY ANALYSIS
		runPathwayAnalysis_ = false;
		geneSetFile_ = null;
		excludedGenesFile_ = null;
		mergeGenesDistance_ = -1;
		useSimulation_ = false;
		useSimulationWeightedSampling_ = false;

		maxNrOfSimulationsForEnrichment_ = 100000;
		useChi2_ = false; 
		useRankSum_ = false; 
		useHypGeom_ = false;
		hypGeomQuantiles_ = null; 

		useGamma_ = false; 
		gammaShapeParameters_ = new ArrayList<Double>();
		gammaShapeParameters_.add(0.1);
		gammaShapeParameters_.add(0.2);
		gammaShapeParameters_.add(0.3);
		gammaShapeParameters_.add(0.5);
		gammaShapeParameters_.add(1.0);
		gammaShapeParameters_.add(2.0);
		gammaShapeParameters_.add(4.0);

		useExpHyp_ = true; 
		expHypParameters_ = new ArrayList<Double>();
		expHypParameters_.add(0.5);
		expHypParameters_.add(0.75);
		expHypParameters_.add(0.85);
		expHypParameters_.add(0.9);
		
		loadScoresFromFiles_ = false; 
		geneScoreFile_ = null; 
		metaGeneScoreFile_ = null; 

		writeSignificanceThreshold_ = -1;
		writeGeneScores_ = false;

		// CONCATENATE CHROMOSOME RESULT FILES
		runConcatenateChromosomeResults_ = false;
		concatenateChromosomeResultsDir_ = null;	
		deleteOriginals_ = false;
		useFakePhenotype_ = false;
		useFakeSignal_ = false;
		chanceOfSignal_=0.0;
		multipleOfPhenotype_=0;

	}
	
	
	// ----------------------------------------------------------------------------

	/** Load and initialize settings */
	public void loadSettings(String settingsFile) {
		
		Pascal.println("SETTINGS FILE");
		Pascal.println("-------------\n");
				
		try {
			// Check that the specified settings file exists
			if (settingsFile == null || settingsFile.isEmpty())
				throw new RuntimeException("No settings file specified");
			else if (!new File(settingsFile).exists())
				throw new RuntimeException("Settings file not found: " + settingsFile);

			// Open file input stream
			Pascal.println("- Loading settings file: " + settingsFile + "\n");
			InputStream in = new FileInputStream(settingsFile);

			// Load the settings
			prop = new Properties();
			prop.load(new InputStreamReader(in));
			
			// Get the param values
			setParameterValues();
			checkOptions();

		} catch (Exception e) {
			Pascal.warning(e.getMessage());
			throw new RuntimeException("Failed to load settings file (a parameter may be missing or malformed): " + settingsFile);
		}		

	}
	
	
	// ----------------------------------------------------------------------------

	
	/** Set Main parameters based on the loaded properties */
	private void setParameterValues() throws Exception {

		// VARIOUS
		if (prop.containsKey("outputDirectory")) {
			outputDirectory_ = getFileSetting("outputDirectory");
			if (outputDirectory_.equals("")) 
				outputDirectory_ = new File(System.getProperty("user.dir"));
		}
		if (prop.containsKey("outputSuffix"))
			outputSuffix_ = getSetting("outputSuffix");
		if (prop.containsKey("verbose"))
			verbose_ = getSettingBoolean("verbose");
		if (prop.containsKey("randomSeed"))
			setRandomSeed(getSettingInt("randomSeed"));

		// ----------------------------------------------------------------------------
		// INPUT

		if (prop.containsKey("snpPvalFile"))
			snpPvalFile_ = getFileSetting("snpPvalFile");
		if (prop.containsKey("pvalCol"))
			pvalCol_ = getSettingInt("pvalCol");
		if (prop.containsKey("withZScore"))
			withZScore_ = getSettingBoolean("withZScore");
		
		if (prop.containsKey("refPopDirectory"))
			refPopDirectory_ = getFileSetting("refPopDirectory");
		if (prop.containsKey("refPopFilePrefix"))
			refPopFilePrefix_ = getSetting("refPopFilePrefix");
		if (prop.containsKey("refPopFileExtension"))
			refPopFileExtension_ = getSetting("refPopFileExtension");
		
		if (prop.containsKey("snpFilterFile"))
			snpFilterFile_ = getFileSetting("snpFilterFile");
		if (prop.containsKey("genesToBeLoadedFile"))
			genesToBeLoadedFile_ = getFileSetting("genesToBeLoadedFile");

		if (prop.containsKey("chromosome"))
			chromosome_ = getSetting("chromosome");
		if (prop.containsKey("ignoreAllosomes"))
			ignoreAllosomes_ = getSettingBoolean("ignoreAllosomes");
		if (prop.containsKey("ignoreChrM"))
			ignoreChrM_ = getSettingBoolean("ignoreChrM");

		if (prop.containsKey("genomeAnnotation"))
			genomeAnnotation_ = getSetting("genomeAnnotation");
		if (prop.containsKey("genecodeAnnotationFile"))
			gencodeAnnotationFile_ = getFileSetting("genecodeAnnotationFile");
		if (prop.containsKey("ucscAnnotationFile"))
			ucscAnnotationFile_ = getFileSetting("ucscAnnotationFile");
		if (prop.containsKey("bedAnnotationFile"))
			bedAnnotationFile_ = getFileSetting("bedAnnotationFile");
		if (prop.containsKey("loadOnlyProteinCodingGenes"))
			loadOnlyProteinCodingGenes_ = getSettingBoolean("loadOnlyProteinCodingGenes");
		if (prop.containsKey("geneIdMappingFile"))
			geneIdMappingFile_ = getFileSetting("geneIdMappingFile");

		if (prop.containsKey("geneWiseSnpWeightsFile"))
			geneWiseSnpWeightsFile_= getFileSetting("geneWiseSnpWeightsFile");
		if (prop.containsKey("weightFileFormat"))
			weightFileFormat_= getSetting("weightFileFormat");
		
		// ----------------------------------------------------------------------------
		// PARAMETERS
		
		if (prop.containsKey("geneWindowUpstream"))
			geneWindowUpstream_ = getSettingInt("geneWindowUpstream");
		if (prop.containsKey("geneWindowDownstream"))
			geneWindowDownstream_ = getSettingInt("geneWindowDownstream");

		if (prop.containsKey("maxSnpsPerGene"))
			maxSnpsPerGene_ = getSettingInt("maxSnpsPerGene");
		if (prop.containsKey("useMafCutoff"))
			useMafCutoff_ = getSettingDouble("useMafCutoff");
		if (prop.containsKey("dePhase"))
			dePhase_ = getSettingBoolean("dePhase");
		if (prop.containsKey("externalConvergenceCheck"))
			externalConvergenceCheck_=getSettingBoolean("externalConvergenceCheck");
		if (prop.containsKey("maxPruningCutoff"))
			maxPruningCutoff_=getSettingDouble("maxPruningCutoff");

		if (prop.containsKey("useAnalyticVegas"))
			useAnalyticVegas_ = getSettingBoolean("useAnalyticVegas");
		if (prop.containsKey("useMaxVegas"))
			useMaxVegas_ = getSettingBoolean("useMaxVegas");
		if (prop.containsKey("useMaxEffVegas"))
			useMaxEffVegas_ = getSettingBoolean("useMaxEffVegas");
		if (prop.containsKey("useSimultationVegas"))
			useSimulationVegas_ = getSettingBoolean("useSimultationVegas");

		if (prop.containsKey("eigenValueFractionCut"))
			eigenValueFractionCut_ = getSettingDouble("eigenValueFractionCut");
		if (prop.containsKey("useOnlyGwasSnps"))
			useOnlyGwasSnps_ = getSettingBoolean("useOnlyGwasSnps");

		if (prop.containsKey("requestedAbsolutePrecision"))
			requestedAbsolutePrecision_ = getSettingDouble("requestedAbsolutePrecision");
		if (prop.containsKey("requestedRelativePrecision"))
			requestedRelativePrecision_ = getSettingDouble("requestedRelativePrecision");
		if (prop.containsKey("gslIntegrationLimit"))
			gslIntegrationLimit_ = getSettingInt("gslIntegrationLimit");

		if (prop.containsKey("farebrotherMaxIterations"))
			farebrotherMaxIterations_ = getSettingInt("farebrotherMaxIterations");
		if (prop.containsKey("farebrotherMode"))
			farebrotherMode_ = getSettingDouble("farebrotherMode");
		if (prop.containsKey("useFarebrother"))
			useFarebrother_ = getSettingBoolean("useFarebrother");
		if (prop.containsKey("useImhof"))
			useImhof_ = getSettingBoolean("useImhof");
		if (prop.containsKey("useDavies"))
			useDavies_ = getSettingBoolean("useDavies");
		if (prop.containsKey("daviesErrorBound"))
			daviesErrorBound_ = getSettingDouble("daviesErrorBound");
		if (prop.containsKey("daviesIntegrationTerms"))
			daviesIntegrationTerms_ = getSettingInt("daviesIntegrationTerms");

		if (prop.containsKey("toleratedAbsolutePrecision"))
			toleratedAbsolutePrecision_ = getSettingDouble("toleratedAbsolutePrecision");
		if (prop.containsKey("toleratedRelativePrecision"))
			toleratedRelativePrecision_ = getSettingDouble("toleratedRelativePrecision");

		if (prop.containsKey("snpWeightingDelta"))
			snpWeightingDelta_ = getSettingDoubleArray("snpWeightingDelta", false, Pascal.log);

		if (prop.containsKey("adaptiveNumSamples"))
			adaptiveNumSamples_ = getSettingIntArray("adaptiveNumSamples", true, Pascal.log);
		if (prop.containsKey("numSamplesGreaterCutoff"))
			numSamplesGreaterCutoff_ = getSettingInt("numSamplesGreaterCutoff");

		if (prop.containsKey("removeCodingSnpsOfOtherGenes"))
			removeCodingSnpsOfOtherGenes_ = getSettingBoolean("removeCodingSnpsOfOtherGenes");
		if (prop.containsKey("codingSnpsFile"))
			codingSnpsFile_ = getSetting("codingSnpsFile");

		// ----------------------------------------------------------------------------
		// OUTPUT
		
		if (prop.containsKey("writeDetailedOutput"))
			writeDetailedOutput_ = getSettingBoolean("writeDetailedOutput");
		if (prop.containsKey("writeDetailedErrorOutput"))
			writeDetailedErrorOutput_ = getSettingBoolean("writeDetailedErrorOutput");
		if (prop.containsKey("writeSnpBedFile"))
			writeSnpBedFile_ = getSettingBoolean("writeSnpBedFile");
		if (prop.containsKey("writeTpedFile"))
			writeTpedFile_ = getSettingBoolean("writeTpedFile");
		if (prop.containsKey("saveGeneReports"))
			saveGeneReports_ = getSettingBoolean("saveGeneReports");
		if (prop.containsKey("writeCorFiles"))
			writeCorFiles_ = getSettingBoolean("writeCorFiles");
		if (prop.containsKey("writeGenewiseSnpFiles"))
			writeGenewiseSnpFiles_=getSettingBoolean("writeGenewiseSnpFiles");
		if (prop.containsKey("writeUsedSettings"))
			writeUsedSettings_ = getFileSetting("writeUsedSettings");

		// ----------------------------------------------------------------------------
		// PATHWAY ANALYSIS

		if (prop.containsKey("runPathwayAnalysis"))
			runPathwayAnalysis_ = getSettingBoolean("runPathwayAnalysis");
		if (prop.containsKey("geneSetFile"))
			geneSetFile_ = getFileSetting("geneSetFile");
		if (prop.containsKey("excludedGenesFile"))
			excludedGenesFile_ = getFileSetting("excludedGenesFile");
		if (prop.containsKey("mergeGenesDistance"))
			mergeGenesDistance_ = getSettingDouble("mergeGenesDistance");

		if (prop.containsKey("useSimulation"))
			useSimulation_ = getSettingBoolean("useSimulation");
		if (prop.containsKey("useSimulationWeightedSampling"))
			useSimulationWeightedSampling_ = getSettingBoolean("useSimulationWeightedSampling");
		if (prop.containsKey("maxNrOfSimulationsForEnrichment"))
			maxNrOfSimulationsForEnrichment_ = getSettingInt("maxNrOfSimulationsForEnrichment");
		if (prop.containsKey("useChi2"))
			useChi2_ = getSettingBoolean("useChi2");
		if (prop.containsKey("useRankSum"))
			useRankSum_ = getSettingBoolean("useRankSum");
		if (prop.containsKey("useHypGeom"))
			useHypGeom_ = getSettingBoolean("useHypGeom");
		if (prop.containsKey("hypGeomQuantiles"))
			hypGeomQuantiles_ = getSettingDoubleArray("hypGeomQuantiles", false, Pascal.log);
		if (prop.containsKey("useGamma"))
			useGamma_ = getSettingBoolean("useGamma");
		if (prop.containsKey("gammaShapeParameters"))
			gammaShapeParameters_ = getSettingDoubleArray("gammaShapeParameters", false, Pascal.log);
		if (prop.containsKey("useExpHyp"))
			useExpHyp_ = getSettingBoolean("useExpHyp");
		if (prop.containsKey("expHypParameters"))
			expHypParameters_ = getSettingDoubleArray("expHypParameters", false, Pascal.log);

		if (prop.containsKey("loadScoresFromFiles"))
			loadScoresFromFiles_ = getSettingBoolean("loadScoresFromFiles");
		if (prop.containsKey("geneScoreFile"))
			geneScoreFile_ = getFileSetting("geneScoreFile");
		if (prop.containsKey("metaGeneScoreFile"))
			metaGeneScoreFile_ = getFileSetting("metaGeneScoreFile");

		if (prop.containsKey("writeSignificanceThreshold"))
			writeSignificanceThreshold_ = getSettingDouble("writeSignificanceThreshold");
		if (prop.containsKey("writeGeneScores"))
			writeGeneScores_ = getSettingBoolean("writeGeneScores");

		// CONCATENATE CHROMOSOME RESULT FILES
		if (prop.containsKey("runConcatenateChromosomeResults"))
			runConcatenateChromosomeResults_ = getSettingBoolean("runConcatenateChromosomeResults");
		if (prop.containsKey("concatenateChromosomeResultsDir"))
			concatenateChromosomeResultsDir_ = getSetting("concatenateChromosomeResultsDir");
		if (prop.containsKey("deleteOriginals"))
			deleteOriginals_ = getSettingBoolean("deleteOriginals");

		
		if (prop.containsKey("useFakePhenotype"))
			useFakePhenotype_=getSettingBoolean("useFakePhenotype");
		if (prop.containsKey("useFakeSignal"))
			useFakeSignal_=getSettingBoolean("useFakeSignal");
		if (prop.containsKey("chanceOfSignal"))
			chanceOfSignal_=getSettingDouble("chanceOfSignal");
		if (prop.containsKey("multipleOfPhenotype"))
			multipleOfPhenotype_=getSettingInt("multipleOfPhenotype");
		
		if (prop.containsKey("deflationRate"))
			deflationRate_ = getSettingDouble("deflationRate");
		if (prop.containsKey("deflationDistance"))
			deflationDist_ = getSettingInt("deflationDistance");
		if (prop.containsKey("onlyPathwayGenesAsBackground"))
			onlyPathwayGenesAsBackground_ = getSettingBoolean("onlyPathwayGenesAsBackground");
		
		// PROCESSED FIELDS
		gwasName_ = Pascal.utils.extractBasicFilename(snpPvalFile_.getName(), false);
		chromFileExtension_ = chromosome_.equals("") ? "" : "." + chromosome_;
	}
	

	// ----------------------------------------------------------------------------

	/** Check if selected options / settings are valid */
	public void checkOptions(){

		if (runPathwayAnalysis_) {
			// At least one enrichment evaluation method must be selected
			if (!useSimulation_ && !useChi2_ && !useRankSum_)
				throw new IllegalArgumentException("Set either 'useRankSum' and/or 'useChi2' and/or useSimulation to true");
		}
		checkBounds(writeSignificanceThreshold_, 0, 1, "writeSignificanceThreshold");

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
	
	
	// ----------------------------------------------------------------------------

	/** Check if the parameters are valid */
	private void checkBounds(double param, double min, double max, String paramName) {

		if (param < min || param > max)
			throw new IllegalArgumentException("Parameter " + paramName + "=" + param + " is outside of its valid range [" + min + ", " + max + "]");
	}

	
	// ----------------------------------------------------------------------------

	/** Create new instances for the random number generators, initialize with randomSeed_ */
	public void setRandomSeed(int seed) {
		
		randomSeed_ = seed;
		if (randomSeed_ == -1) {
			//mersenneTwisterRng_ = new MersenneTwister();
			wellRng_ = new Well19937c();
			jdkRng_ = new Random();
		} else {
			//mersenneTwisterRng_ = new MersenneTwister(randomSeed_);
			wellRng_ = new Well19937c(randomSeed_);
			jdkRng_ = new Random(randomSeed_);
		}
	}
	
	public int getRandomSeed() { return randomSeed_; }


}
