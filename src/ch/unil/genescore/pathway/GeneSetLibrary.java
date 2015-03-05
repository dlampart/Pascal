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
package ch.unil.genescore.pathway;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import javastat.inference.nonparametric.RankSumTest;
import jsc.independentsamples.MannWhitneyTest;
import jsc.tests.H1;
import no.uib.cipr.matrix.DenseMatrix;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.DistributionMethods;
import ch.unil.genescore.vegas.WritingMethods;


/**
 * The gene sets
 */
public class GeneSetLibrary {

	
	/** The gene sets for which enrichment is computed */
	protected ArrayList<GeneSet> geneSets_ = null;
	//TODO: genes shouldn't contain meta-genes anymore
	/** The background set of genes, i.e., the union of all gene sets. Does NOT contain metaGenes. Does NOT contain genes without score. */
	protected HashSet<Gene> genes_ = null;
	/** The union of all meta genes, ordered by chromosome and start position of the first gene geneId is key. States are as follows:
	 * first, metaGenes are assembled via geneSet.createMetaGenes().
	 * second, metaGenes are calculated.
	 * third, metaGenes are removed again if they had no score.
	 *  */
	protected HashMap<String,MetaGene> metaGenes_ = null;
	
	protected DenseMatrix pathwayCorMat_ = null;
	private static Random rand = new Random();


	protected  ArrayList<Gene> genesForSimulation_ = null;
	/** Genes that were not loaded from the pathway / gene set file because they are not in the loaded genome annotation */
	private HashSet<String> notFoundInGenomeAnnot_ = null;
	
	/** The number of gene sets with either p-value below Settings.pathwaySignificanceThreshold */
	private int numSignificantSets_ = -1;
	/** key: gene_ids ; values: all pws that contain this gene (or metagene containing that gene.)  */
	private HashMap<String, HashSet<GeneSet>> samplingWeightHelper_ = null;
	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor */
	public GeneSetLibrary() {
	}
	
	/**input: genes: string is geneId within gene*/
	public GeneSetLibrary(String geneSetFile, HashMap<String, Gene> genes) {
		//TODO: does it really here cover all genes that do not have a score?, no
		load(geneSetFile, genes);
	}

	
	// ----------------------------------------------------------------------------

	/** 
	 * Create meta-genes genes and update genes_, which is the union of all gene sets
	 * (including meta-genes, excluding genes that do not appear individually [outside
	 * of meta-genes] anymore). genes_ is ordered by chromosome and position.
	 */
	public void createMetaGenes() {
		
		// The union of all meta-genes
		metaGenes_ = new HashMap<String, MetaGene>();		
		for (GeneSet set : geneSets_) {
			set.createMetaGenes(metaGenes_);
			for (MetaGene myMetaGene : set.getMetaGenes()){
			
			metaGenes_.put(myMetaGene.getId(), myMetaGene);
			}
		}

	}
	/** removes metaGenes that didn't have a score. has to be removed 
	 * from metaGenes_ and geneSets_.
	 * */
	public void removeMetaGenesWithoutScore() {
				
		HashSet<MetaGene> metaGenesToBeRemoved = new HashSet<MetaGene>();
		
		for (GeneSet set : geneSets_) {
			
			metaGenesToBeRemoved.addAll(set.removeMetaGenesWithoutScore());
		}
		for (MetaGene metaGeneToBeRemoved : metaGenesToBeRemoved)
			metaGenes_.remove(metaGeneToBeRemoved.id_);
		
	}
	
	
	// ----------------------------------------------------------------------------

	/** Count min, max and total number of merged genes per meta-gene */
	public int[] countMergedGenes() {
		
		// Initialize
		int[] minMaxTot = new int[3];
		minMaxTot[0] = Integer.MAX_VALUE;
		minMaxTot[1] = Integer.MIN_VALUE;
		minMaxTot[2] = 0;
		
		for (GeneSet set : geneSets_)
			set.countMergedGenes(minMaxTot);
		
		return minMaxTot;
	}

		
	// ----------------------------------------------------------------------------
	protected double[] getAllChisq(){
		
		double[] ret = new double[genesForSimulation_.size()];
		int i=0;
		for (Gene g: genesForSimulation_){
			ret[i]=g.getChi2Stat();
			i++;
		}
		return ret;		
	}
	
	protected double[] getAllSamplingWeights(){
		
		double[] ret = new double[genesForSimulation_.size()];
		int i=0;
		for (Gene g: genesForSimulation_){
			ret[i]=g.getSamplingWeight();
			i++;
		}
		return ret;		
	}
	
	
	protected double[] getSetChisq(GeneSet set){
	
		HashSet<Gene> gSet=set.getGenes();
		double[] ret= new double[gSet.size()];
					
		// Sum the chi2 stats
		int i = 0;
		for (Gene g : gSet){
			ret[i]=g.getChi2Stat();
			i++;
		}
		return ret;	
	}
	
	
	private void returnChiSqToPreviousVal(double[] oldVals){
		
		int i=0;
		for (Gene g: genesForSimulation_){			
			g.setChi2Stat(oldVals[i]);
			i++;
		}
	}
	
	/** Compute enrichment for each gene set */
	public void computeEnrichment() {
		
		calculateSamplingWeights();
		
		double[] allChiSqVals=null;
		double[] allChiSqValsCopy=null;
		double[] setChiSqVals=null;
		
		double[] allSamplingWeights = null;
		
		//if (Settings.useSimulation_){
			double p;
			double[] ps= null;
			for (GeneSet set : geneSets_){
				System.out.println(set.getId());
				fillUpGenesForSimulationArray(set);
				allChiSqVals=getAllChisq();
				allChiSqValsCopy=allChiSqVals.clone();
				if(Settings.deflationRate_!=0){
					
					Damper damper = new Damper(genesForSimulation_, Settings.deflationDist_, Settings.deflationRate_);
					damper.dampSet();
				}
				allChiSqVals=getAllChisq();
				setChiSqVals=getSetChisq(set);
				allSamplingWeights=getAllSamplingWeights();
				
				if (Settings.useRankSum_){
					computeRankSumPvalue(set);
				}
				
				if (Settings.useSimulation_){
					System.out.println("simulation");
					p = simulatePvalue(setChiSqVals,allChiSqVals);
					set.setSimulPvalue(p);
				}
				if (Settings.useSimulationWeightedSampling_){
					System.out.println("simulationWeightedSampling");
					p = simulatePvalueWeightedSampling(setChiSqVals, allChiSqVals, allSamplingWeights);
					set.setSimulPvalueWeighted(p);
				}
				
				if (Settings.useHypGeom_){
					computeHypGeomPvalue(set);
					ps=computeHypGeomPvalue2(setChiSqVals,allChiSqVals);
					
				}
				if (Settings.useChi2_ || Settings.useGamma_ || Settings.useExpHyp_) {
					rankGeneScoresAndMapToChi();				
					fillUpGenesForSimulationArray(set);
					if (Settings.useChi2_)
						computeChi2Pvalue(set);
				}		
	
				returnChiSqToPreviousVal(allChiSqValsCopy);
			}
		//}		
		
	}
	
	public void computeApproxChi2Values(){
		
		Percentile perc = new Percentile();
		EffTestCalculator effCalc = new EffTestCalculator();
		effCalc.addGeneIds(genes_);
		effCalc.addMetaGeneIds(metaGenes_);
		effCalc.setGeneSets(geneSets_);
		effCalc.initializeGeneVals();
		
		int nrOfLoops=10000;
		double[] mins = new double[nrOfLoops];
		for (int i=0; i<nrOfLoops;i++){			
			mins[i]=effCalc.calcMinVal();
		}
		perc.setData(mins);
		double p1 = perc.evaluate(1);
		double p5 = perc.evaluate(5);
		double p10 = perc.evaluate(10);
		double p15 = perc.evaluate(15);
		//double w1 = perc.evaluate(1.0/50);
		System.out.println(p1);
		System.out.println(p5);
		System.out.println(p10);
		System.out.println(p15);
		System.out.println("asdf");
		
	}

	// ----------------------------------------------------------------------------


	/** Write results of enrichment computation */
	public void writeOutput(String filename) {

		// Sort gene sets based on enrichment score
		sortGeneSets();

		// Open output file
		FileExport writer = new FileExport(filename);
		// Write header
		writer.println(GeneSet.getResultsAsStringHeader());

		for (GeneSet set : geneSets_) {
			if ((Settings.useChi2_ && set.getChi2Pvalue() <= Settings.writeSignificanceThreshold_) ||
					(Settings.useSimulation_ && set.getSimulPvalue() <= Settings.writeSignificanceThreshold_)
						|| (Settings.useRankSum_ && set.getRankSumPvalue() <= Settings.writeSignificanceThreshold_))
				writer.println(set.getResultsAsString());
		}
		writer.close();
	}
	public void writePathwayCorrelationMat(String fileName, String additonalDirectory){
		WritingMethods.writeMTJ(pathwayCorMat_,fileName, additonalDirectory);
	}
	
	// ----------------------------------------------------------------------------

	/** Sort gene sets by enrichment p-value */
	public void sortGeneSets() {
		
		class GeneSetComparator implements Comparator<GeneSet> {
			@Override
			public int compare(GeneSet set1, GeneSet set2) {
				if (Settings.useChi2_)
					return Double.compare(set1.getChi2Pvalue(), set2.getChi2Pvalue());
				else
					return Double.compare(set1.getSimulPvalue(), set2.getSimulPvalue());
			}
		}
		Collections.sort(geneSets_, new GeneSetComparator());
	}

	
	// ----------------------------------------------------------------------------

	/** Count number of sets with significant enrichment or depletion */
	public int countNumSignificantSets() {
		
		numSignificantSets_ = 0;
		for (GeneSet set : geneSets_)
			if ((Settings.useChi2_ && set.getChi2Pvalue() <= Settings.writeSignificanceThreshold_) ||
					(Settings.useSimulation_ && set.getSimulPvalue() <= Settings.writeSignificanceThreshold_)) {
				numSignificantSets_++;
		}
		return numSignificantSets_;
	}

	
	// ============================================================================
	// PRIVATE METHODS
	
	/** Transform gene scores to ranks and map them to the corresponding quantiles of 1-df chi2 distribution */
	private void rankGeneScoresAndMapToChi() {
		
		// Transform scores to empirical probabilities and map to 1-df chi2 distribution
		//ChiSquaredDistribution chi2 = new ChiSquaredDistribution(1);
		double prevScore = -1;
		int count = 1;

		// Sort genes by score
		GeneScoreList rankedGenes = new GeneScoreList(genesForSimulation_, true);
		
		for (Gene g : rankedGenes.getGenes()) {
			// Check that genes are ordered
			double score[] = g.getScore();
			if (score[0] < prevScore)
				throw new RuntimeException("Genes are not ordered by score");
			prevScore = score[0];
			
			double empiricalProb = count / ((double) genesForSimulation_.size()+1);
			g.setNormalizedScore(empiricalProb);
			g.setChi2Stat(DistributionMethods.chiSquared1dfInverseCumulativeProbabilityUpperTail(empiricalProb));
			count++;
		}
	}
	
	/** Transform gene scores to ranks and map them to the corresponding quantiles of 1-df chi2 distribution */
	private void rankGeneScoresAndMapToUniform() {
		
		// Transform scores to empirical probabilities and map to 1-df chi2 distribution
		//ChiSquaredDistribution chi2 = new ChiSquaredDistribution(1);
		double prevScore = -1;
		int count = 1;

		// Sort genes by score
		GeneScoreList rankedGenes = new GeneScoreList(genesForSimulation_, true);
		
		for (Gene g : rankedGenes.getGenes()) {
			// Check that genes are ordered
			double score[] = g.getScore();
			if (score[0] < prevScore)
				throw new RuntimeException("Genes are not ordered by score");
			prevScore = score[0];
			
			double empiricalProb = count / ((double) genesForSimulation_.size()+1);
			g.setNormalizedScore(empiricalProb);
			count++;
		}
	}
	
	public static List<Gene> getRandomSample(ArrayList<Gene> items, int m){
	    for(int i=0;i<m;i++){
	        int pos = i + rand.nextInt(items.size() - i);
	        Gene tmp = items.get(pos);
	        items.set(pos, items.get(i));
	        items.set(i, tmp);
	    }
	    return  items.subList(0, m);
	}
	
	public static double[] getRandomSampleDouble(double[] items, int m){
	    for(int i=0;i<m;i++){
	        int pos = i + rand.nextInt(items.length - i);
	        double tmp = items[pos];
	        items[pos]=items[i];
	        items[i]=tmp;
	    }
	    double[] subList = new double[m];
	    for (int i=0; i<m; i++)
	    	subList[i]=items[i];
	    return  subList;
	}
	
	/** Compute enrichment for the given set using hypergeometric distribution (magenta-option)
	 * 	quantile refers to the cutoff used to assign ins and outs.
	 * */
	protected double[] computeHypGeomPvalue2(double[] set, double[] totSet) {
		double[] quantiles = Settings.hypGeomQuantiles_;
		double[] pvals = new double[quantiles.length];
		if (set.length == 0){
			for (int i=0;i<quantiles.length;i++){
				pvals[i]=1;
			}						
			return pvals;
		}
			
		for (int j=0;j <quantiles.length ; j++){
			double[] valAr = totSet.clone();
			Arrays.sort(valAr);	
			int n=(int) Math.floor(valAr.length*quantiles[j]);
			double cutoffVal=valAr[n];
			int setSize= valAr.length-n-1;
					
		int q = 0;
		for (double g : set){
			if(g >cutoffVal)
				q += 1;
		}
		HypergeometricDistribution myHypgeom = new HypergeometricDistribution(valAr.length,setSize, set.length);
		double pval=1-myHypgeom.cumulativeProbability((q-1));		
		pvals[j]=pval;
		}
		return(pvals);			
	}
		
	
	
protected void computeHypGeomPvalue(GeneSet set) {
		double[] quantiles = Settings.hypGeomQuantiles_;
		double[] pvals = new double[quantiles.length];
		if (set.getGenes().size() == 0){
			for (int i=0;i<quantiles.length;i++){
				pvals[i]=1;
			}			
			set.setHypGeomPvalues(pvals);
			return;
		}
			
		for (int j=0;j <quantiles.length ; j++){
		double[] valAr = new double[genesForSimulation_.size()];
		for (int i=0;i<genesForSimulation_.size();i++){
			valAr[i]=genesForSimulation_.get(i).getChi2Stat();
		}
		Arrays.sort(valAr);
		int n=(int) Math.floor(valAr.length*quantiles[j]);
		double cutoffVal=valAr[n];
		int setSize= valAr.length-n-1;
		HashSet<Gene> pathwayGenes = set.getGenes();
		if (pathwayGenes.size() == 0)
			return;
				

		int q = 0;
		for (Gene g : pathwayGenes){
			if(g.getChi2Stat()>cutoffVal)
				q += 1;
		}
		
		
		HypergeometricDistribution myHypgeom = new HypergeometricDistribution(valAr.length,setSize, pathwayGenes.size());
		double pval=1-myHypgeom.cumulativeProbability((q-1));
		pvals[j]=pval;
		}
		set.setHypGeomPvalues(pvals);
	}
/** get a weighted random sample; implementing Efraimidis et al. 2006 */
double[] getWeightedRandomSample(int setLength, double[] totSet,double[] totWeights){
	//NaturalRanking ranker = new NaturalRanking();
	double[] draws = new double[totWeights.length];
	int[] ranks = new int[totWeights.length];
	double[] out = new double[setLength];
	TreeMap<Double,Integer> myRankTree = new TreeMap<Double,Integer>();
	int treesize=0;
	for (int i=0 ; i<totWeights.length ; i++){
		draws[i]=Math.log(rand.nextDouble())/totWeights[i];
		if (treesize<setLength){
			treesize++;
			myRankTree.put(draws[i], i);
		}
		else if(treesize==setLength)
			if(myRankTree.firstKey() < draws[i]){
				myRankTree.pollFirstEntry();
				myRankTree.put(draws[i], i);
			}	
	}
	Iterator<Entry<Double, Integer>> it = myRankTree.entrySet().iterator();
	Entry<Double, Integer> ent;
	int count=0;
	while (it.hasNext()){		
		ent = it.next();
		ranks[count]=ent.getValue();
		count++;
		//int rank = it.;
	}
	for (int i=0 ; i<setLength ; i++){
		
		out[i]=totSet[ranks[i]];
		
	}	
	return out;
	
}
/** Compute enrichment for the given set using chi2 distribution using weighted sampling implementing Efraimidis et al. 2006 */
protected double simulatePvalueWeightedSampling(double[] set, double[] totSet,double[] totWeights) {

	
	if (set.length == 0)
		return 1;
			
	// Sum the chi2 stats
	double q = 0;
	for (double g : set)
		q += g;
	
	double nrOfHigherSimuls=0;
	boolean enoughSimulated=false;
	int nrOfSimulRuns=1000;
			if (nrOfSimulRuns*100>Settings.maxNrOfSimulationsForEnrichment_){
				throw new RuntimeException("minimum nr of simulation runs for enrichment is 100000. at least 1000'000 is recommended.");
			}
	while (!enoughSimulated && nrOfSimulRuns <= Settings.maxNrOfSimulationsForEnrichment_){
		System.out.println(nrOfSimulRuns);
		nrOfSimulRuns=nrOfSimulRuns*10;
		nrOfHigherSimuls=0;
		
		int aboveZero = 0;
		for (int i=0; i<totWeights.length; i++)
			if (totWeights[i]>0)
				aboveZero++;
				
			
		double[] totSetPositiveWeight = new double[aboveZero];
		double[] totWeightsPositiveWeight = new double[aboveZero];
		int count=0;
		for (int i=0; i<totWeights.length; i++)
			if (totWeights[i]>0){
				totSetPositiveWeight[count]=totSet[i];
				totWeightsPositiveWeight[count]=totWeights[i];
				count++;
			}
		
		
		double[] rndSample = null;
		
		for (int i=0 ; i < nrOfSimulRuns ; i++){
			
		//	System.out.println(i);
			rndSample=getWeightedRandomSample(set.length, totSetPositiveWeight,totWeightsPositiveWeight);

			double qsimul=0;
			for (double g : rndSample)
				qsimul += g;
			if (qsimul >= q)
				nrOfHigherSimuls += 1;		
			
			}
		if (nrOfHigherSimuls>50){
			enoughSimulated=true;
		}
	}	
	double pval = nrOfHigherSimuls/nrOfSimulRuns;
		
	//set.setSimulPvalue(pval);
	return pval;
	
}
	/** Compute enrichment for the given set using chi2 distribution */
protected double simulatePvalue(double[] set, double[] totSet) {	
		
		if (set.length == 0)
			return 1;
				
		// Sum the chi2 stats
		double q = 0;
		for (double g : set)
			q += g;
		
		double nrOfHigherSimuls=0;
		boolean enoughSimulated=false;
		int nrOfSimulRuns=10000;
				if (nrOfSimulRuns*100>Settings.maxNrOfSimulationsForEnrichment_){
					throw new RuntimeException("minimum nr of simulation runs for enrichment is 100000. at least 1000'000 is recommended.");
				}
		while (!enoughSimulated && nrOfSimulRuns <= Settings.maxNrOfSimulationsForEnrichment_){
			System.out.println(nrOfSimulRuns);
			nrOfSimulRuns=nrOfSimulRuns*10;
			nrOfHigherSimuls=0;
			double[] rndSample = null;
			
			for (int i=0 ; i < nrOfSimulRuns ; i++){
				
				rndSample=getRandomSampleDouble(totSet, set.length);

				double qsimul=0;
				for (double g : rndSample)
					qsimul += g;
				if (qsimul >= q)
					nrOfHigherSimuls += 1;		
				
				}
			if (nrOfHigherSimuls>50){
				enoughSimulated=true;
			}
		}	
		double pval = nrOfHigherSimuls/nrOfSimulRuns;
			
		//set.setSimulPvalue(pval);
		return pval;
	}
	
	// ----------------------------------------------------------------------------
/** Compute enrichment for the given set using rank-sum-Test */
private void computeRankSumPvalue(GeneSet set) {

	Set<Gene> complementSet = new HashSet<Gene>(genesForSimulation_);
	complementSet.removeAll(set.genes_);
	double[] complementVals = new double[complementSet.size()];
	double[] setVals = new double[set.genes_.size()];
	int count=0;
	for (Gene g : set.genes_){
		setVals[count]=g.getChi2Stat();
		count++;	
	}
	count=0;
	for (Gene g : complementSet){
		complementVals[count]=g.getChi2Stat();
		count++;	
	}	
//	
	double pval=1;
	double pval6=1;
	
	
	MannWhitneyUTest test = new MannWhitneyUTest();
	try{
		pval = test.mannWhitneyUTest(setVals,complementVals);
		RankSumTest test2 = new RankSumTest();
		double pval2=test2.pValue("greater",setVals,complementVals);
		double pval3=test2.pValue("equal",setVals,complementVals);
		MannWhitneyTest test3 = new	MannWhitneyTest(setVals,complementVals);
		double pval5=test3.approxSP();			
		MannWhitneyTest test4 = new	MannWhitneyTest(setVals,complementVals, H1.GREATER_THAN);
		 pval6=test4.approxSP();		
	}catch(Exception e){
		
		System.out.println("Couldn't compute ManWhitneyU:"+ e.getMessage());
	}
	
	
	set.setRankSumPvalue(pval6);
}


	/** Compute enrichment for the given set using chi2 distribution */
	private void computeChi2Pvalue(GeneSet set) {

		HashSet<Gene> pathwayGenes = set.getGenes();
		if (pathwayGenes.size() == 0)
			return;
				
		// Sum the chi2 stats
		double q = 0;
		for (Gene g : pathwayGenes)
			q += g.getChi2Stat();
			//q += g.getScore(0);
		
		ChiSquaredDistribution chi2 = new ChiSquaredDistribution(pathwayGenes.size());
		// P(X > q) -- enrichment of genes with low p-values
		double pval = 1 - chi2.cumulativeProbability(q);
//		boolean depletion = false;
		
//		if (pval > 0.5) {
			//pval = 1 - pval;
//			depletion = true;
//		}
		set.setChi2Pvalue(pval);
//		set.setDepletion(depletion);
	}
	
	// ----------------------------------------------------------------------------

	/** Load gene set library, only genes that have scores are included */
	private void load(String geneSetFile, HashMap<String, Gene> genes) {

		geneSets_ = new ArrayList<GeneSet>();
		
		genes_ = new HashSet<Gene>();
		notFoundInGenomeAnnot_ = new HashSet<String>();	
		if(!Settings.onlyPathwayGenesAsBackground_){
			for (Map.Entry<String, Gene> entry : genes.entrySet()){
				if(entry.getValue().getChi2Stat()!=-1)
					genes_.add(entry.getValue());
			}
		}
		// The column where the gene list starts
		int firstGeneCol = 1; // standard format: id gene1 gene2 ...
		if (geneSetFile.endsWith(".gmt"))
			firstGeneCol = 2; // gmt format (msigdb): id url gene1 gene2 ...
		
		// Open the file
		FileParser parser = new FileParser(geneSetFile);

		// Check that gene set ids are unique
		HashSet<String> geneSetIds = new HashSet<String>();
		
		while (true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;

			// The next gene set
			String geneSetId = nextLine[0];
			// Check that ids are unique
			if (geneSetIds.contains(geneSetId))
				throw new RuntimeException("Gene set " + geneSetId + " is listed twice");
			
			// Create the gene set
			GeneSet geneSet = new GeneSet(geneSetId);
			geneSets_.add(geneSet);
			geneSetIds.add(geneSetId);
			
			// Add the genes in pathways to gene sets
			for (int i=firstGeneCol; i<nextLine.length; i++) {
				String id = nextLine[i];
				Gene gene = genes.get(id);
				
				if (gene != null && gene.getChi2Stat()!=-1) {					
					geneSet.add(gene);
					if(Settings.onlyPathwayGenesAsBackground_){
						genes_.add(gene);
					}
				
				} else {
					// The gene has no score, skip
					notFoundInGenomeAnnot_.add(id);
				}
			}
			
		}
		
	
		parser.close();
	}

	/**genesForSimulation */
	private void fillUpGenesForSimulationArray(GeneSet geneSet){
		//if("KEGG_CITRATE_CYCLE_TCA_CYCLE".equals(geneSet.getId())){
		//	System.out.println("asfasf");
		//}
		//get all genes contained in meta genes.
		HashSet<Gene> genesInMetaGenes= new HashSet<Gene>();
		if (metaGenes_!=null){
			for (MetaGene metaGene : geneSet.getMetaGenes()){			
				TreeSet<Gene> allGenesInMetaGene = metaGene.getGenes();
				for (Gene gene : allGenesInMetaGene){
					genesInMetaGenes.add(gene);
				}
			}
		}
		//add genes not belonging to meta-gene.
		genesForSimulation_ = new ArrayList<Gene>();
		for (Gene gene : genes_){
			if (!genesInMetaGenes.contains(gene)){
				genesForSimulation_.add(gene);
			}
		}
		//add meta genes.
		if (geneSet.getMetaGenes()!=null){
			for (MetaGene metaGene : geneSet.getMetaGenes()){
			genesForSimulation_.add(metaGene);
			}		
		}
	}	
	// ----------------------------------------------------------------------------
	/** calculates the weights  for sampling the genes:
	 * 
	 * for each gene x:
	 * 1) we use samplingWeightHelper_ to find all gene sets where the gene is present (as gene or as part of metagene)
	 * 2) for each of these genesets : we update gene x's samplingWeight by adding 1/(sizeofGeneSet^2)
	 * 
	 * for each metaGene : replace 1) above with :find all gene sets where any uf the subgene is present.
	 *  */
	protected void calculateSamplingWeights(){
		
		calculateSamplingWeightHelper();
		HashSet<GeneSet> gSetSet = null;
		double size;
		double weight;
		for (Gene g: genes_){				
			if(getSamplingWeightHelper().containsKey(g.id_)){
				gSetSet = getSamplingWeightHelper().get(g.id_);
				for (GeneSet gSet : gSetSet){
					size = gSet.genes_.size();
					weight = 1.0/(size*size);
					g.updateSamplingWeight(weight);
				}
			}	
		}		
		if(metaGenes_!=null){
			for (MetaGene mg: metaGenes_.values()){
				gSetSet = new HashSet<GeneSet>();
				for (Gene g: mg.getGenes()){
				
					gSetSet.addAll(getSamplingWeightHelper().get(g.id_));						
				}
				for (GeneSet gSet : gSetSet){
					size = gSet.genes_.size();
					weight = 1.0/(size*size);
					mg.updateSamplingWeight(weight);
				}
			}
		}
		
			
	}
	/**
	 * fills up  samplingWeightHelper_
	 * 
	 * key: gene_ids ; values: all pws that contain this gene (or metagene containing that gene.)  */
	private void calculateSamplingWeightHelper(){
		HashMap<String, HashSet<GeneSet>> metaSet = new HashMap<String, HashSet<GeneSet>>();
		//Gene geneToProcess =null;		
		TreeSet<Gene> subGenes = null;
		for (GeneSet set : geneSets_){					
			for (Gene g : set.genes_){
				if (g instanceof MetaGene){
					MetaGene mg = (MetaGene) g;
					subGenes  = mg.getGenes();
				}
				else {
					subGenes = new TreeSet<Gene>();
					subGenes.add(g);
				}
						
					//for (Gene subg : (MetaGene) g.getMetaGenes()))
				for (Gene geneToProcess : subGenes){
				
					if (!metaSet.containsKey(geneToProcess.id_)){
						HashSet<GeneSet> setForG  = new HashSet<GeneSet>();
						setForG.add(set);
						metaSet.put(geneToProcess.id_, setForG);
					}
					else{
						metaSet.get(geneToProcess.id_).add(set);				
					}
				}		
			}
		}
		setSamplingWeightHelper(metaSet);		
	}
	// ----------------------------------------------------------------------------

	/** Create a vector with the scores of these genes */
		private void printSimulationArray(GeneSet set) {

			if(set.getMetaGenes()!=null || set.getMetaGenes().isEmpty()){
			FileExport writer = new FileExport("geneSimulFiles/SimulationAr_" + set.getId() + ".txt");
			// Write header
			for (Gene el:genesForSimulation_){
				
				writer.println(el.id_+"\t" + el.getChi2Stat());
				
			}
			writer.close();
			}
			//double[] scores = new double[genes.size()];
			//int i = 0;
			//for (Gene g : genes)
			//	scores[i++] = g.getScore(0);
			
			//return scores;
		}
	// ----------------------------------------------------------------------------

	/** Create a vector with the scores of these genes */
	private double[] asScoreVector(Collection<Gene> genes) {

		double[] scores = new double[genes.size()];
		int i = 0;
		for (Gene g : genes)
			scores[i++] = g.getScore(0);
		
		return scores;
	}

	public void computeApproxPathwayCorrelation() {
		
		DenseMatrix corMat = new DenseMatrix(geneSets_.size(),geneSets_.size());
		for (int i=0; i<geneSets_.size();i++){
			GeneSet leftSet = geneSets_.get(i);
			double leftSize = leftSet.genes_.size();
			for (int j=0; j<geneSets_.size();j++){
				GeneSet rightSet = geneSets_.get(j);
				double rightSize = rightSet.genes_.size();
				HashSet<Gene> unpackedMetaGenes = new HashSet<Gene>();
				HashSet<Gene> allRightGenes = new HashSet<Gene>();
				if (null!=rightSet.getMetaGenes())
				for (MetaGene mg : rightSet.getMetaGenes()){
					unpackedMetaGenes.addAll(mg.getGenes());
				}
				
				allRightGenes.addAll(unpackedMetaGenes);
				allRightGenes.addAll(rightSet.genes_);
				allRightGenes.removeAll(rightSet.getMetaGenes());
				
				HashSet<Gene> copiedLeftGenes = new HashSet<Gene>(leftSet.genes_);
				copiedLeftGenes.retainAll(allRightGenes);
				double count=copiedLeftGenes.size();
				if (null!=leftSet.getMetaGenes())
				for (MetaGene mg : leftSet.getMetaGenes()){
					TreeSet<Gene> mgSetCopy =new TreeSet<Gene>(mg.getGenes());
					mgSetCopy.retainAll(allRightGenes);
					if (!mgSetCopy.isEmpty()){
						count++;
					}
				}				
				double corr = count/Math.sqrt(leftSize*rightSize);
				corMat.set(i, j, corr);
				//corMat.set(j, i, corr);
			}			
		}
		pathwayCorMat_ = corMat;
	}
	
public void computeApproxPathwayCorrelationLowerBound() {
		
		DenseMatrix corMat = new DenseMatrix(geneSets_.size(),geneSets_.size());
		for (int i=0; i<geneSets_.size();i++){
			HashSet<Gene> leftSet = geneSets_.get(i).genes_;
			double leftSize = leftSet.size();
			for (int j=0; j<geneSets_.size();j++){
				HashSet<Gene> rightSet = geneSets_.get(j).genes_;
				double rightSize = rightSet.size();
				HashSet<Gene> copyLeftSet= new HashSet<Gene>(leftSet);
				copyLeftSet.retainAll(rightSet);
				double count = copyLeftSet.size(); 
				double corr = count/Math.sqrt(leftSize*rightSize);
				//double corr = count;
					
				corMat.set(i, j, corr);
			}
		}
		pathwayCorMat_ = corMat;
	}
	// ============================================================================
	// GETTERS AND SETTERS
	
	public ArrayList<GeneSet> getGeneSets() { return geneSets_; }
	public HashSet<Gene> getGenes() { return genes_; }
	public Collection<MetaGene> getMetaGenes() { return metaGenes_.values(); }
	
	public HashSet<String> getNotFoundInGenomeAnnot() { return notFoundInGenomeAnnot_; }
	public int getNumSignificantSets() { return numSignificantSets_; }
	public void setGenesForSimulation(ArrayList<Gene> genesForSimulation){
	genesForSimulation_=genesForSimulation;
	}

	public HashMap<String, HashSet<GeneSet>> getSamplingWeightHelper() {
		return samplingWeightHelper_;
	}

	public void setSamplingWeightHelper(HashMap<String, HashSet<GeneSet>> samplingWeightHelper_) {
		this.samplingWeightHelper_ = samplingWeightHelper_;
	}

}