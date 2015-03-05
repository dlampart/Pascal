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

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Settings;
import no.uib.cipr.matrix.DenseVector;

public class GeneDataFakePhenotype extends GeneData {
	
	static int nrOfMultiples_;
	static double[] fakePheno_=null;
	static double[][] fakePhenoAr_=null;
	static double fakePhenoSd_;
	static double fakePhenoMean_;
	//static Random r = new Random();
	static String snpStr= "snpValsSeed"+ Settings.randomSeed_ + ".txt";
	static String fakePhenoStr= "fakePhenoSeed"+ Settings.randomSeed_ + ".txt";
	static FileExport writerSnpScores = new FileExport(snpStr);
	static FileExport writerFakePheno = new FileExport(fakePhenoStr);
	private static Random rand = new Random(Settings.randomSeed_);
	static TDistribution myT= null;
	static double[] fakeSignal_ = null;
	public GeneDataFakePhenotype(ArrayList<Snp> snpList) {
		super(snpList);				
		nrOfMultiples_=Settings.multipleOfPhenotype_;		
	}
	
	public GeneDataFakePhenotype(ArrayList<Snp> snpList, double[] fakeSignal) {
		super(snpList);	
		fakeSignal_=fakeSignal;				
		nrOfMultiples_=Settings.multipleOfPhenotype_;
		
	}	
	
	private void recalcPvalsAndZscoresMultipleTimes(Snp mySnp){
		
		double zscoreSum = 0;
		for (int i=0; i<nrOfMultiples_;i++){
			mySnp.recalcPvalsAndZscoresDirect(fakePhenoAr_[i]);
			zscoreSum = zscoreSum + mySnp.getZscore();
		}
		double zscore = zscoreSum/Math.sqrt(nrOfMultiples_);
		double chi2Stat=zscore*zscore;
		//double pval=1-chiSquared1df_.cumulativeProbability(chi2Stat);
		double pval=DistributionMethods.chiSquared1dfCumulativeProbabilityUpperTail(chi2Stat);				
		mySnp.setChi2Stat(chi2Stat);
		mySnp.setPval(pval);
		mySnp.setZscore(zscore);
		
	}
	
	public void processData(){
		
		if (fakePhenoAr_==null){
			//constructFakePhenotype();
			//if (Settings.multipleOfPhenotype_>1)
				constructFakePhenotypeAr();
			
			//writeFakePhenotype();
		}			
			for (Snp mySnp : snpList_){
				//mySnp.recalcPvalsAndZscores(fakePheno_, fakePhenoMean_,fakePhenoSd_, myT);
				//mySnp.recalcPvalsAndZscoresDirect(fakePheno_);
				recalcPvalsAndZscoresMultipleTimes(mySnp);
				String writingString = Double.toString(mySnp.getZscore());
				writerSnpScores.println(writingString);
				
			}
							
		 cov_ = LinkageDisequilibrium.computeCorrelationMatrixMTJ(snpList_);
		setSnpScores();
	 }
	
	public void constructFakePhenotype(){		
		int nrOfSamples = snpList_.get(0).getGenotypes().length;
		fakePheno_ = new double[nrOfSamples];
		simulateNormals(fakePheno_);
		if (fakeSignal_!=null)
			addSignal(fakePheno_);
		double len = fakePheno_.length-2;
		myT = new TDistribution(len);
	}
	public void constructFakePhenotypeAr(){			
		int nrOfSamples = snpList_.get(0).getGenotypes().length;
		fakePhenoAr_ = new double[nrOfMultiples_][nrOfSamples];
		for (int i=0; i<nrOfMultiples_;i++)
			simulateNormals(fakePhenoAr_[i]);
		if (fakeSignal_!=null){
			for (int i=0; i<nrOfMultiples_;i++)
				addSignal(fakePhenoAr_[i]);
		}
		double len = fakePhenoAr_[0].length-2;
		myT = new TDistribution(len);
	}
	private void calcFakePhenoMeanAndSd(){
		int lmin1=fakePheno_.length-1;
		double adj=((double) fakePheno_.length)/((double) lmin1);
		double curSum=0;
		double curSumSq=0;
		for (double d : fakePheno_){
			curSum=curSum+d;
			curSumSq=curSumSq+d*d;
		}
		fakePhenoMean_=curSum/((double) fakePheno_.length);		
		double curMeanSq=curSumSq/((double) fakePheno_.length);
		fakePhenoSd_=curMeanSq-(fakePhenoMean_*fakePhenoMean_);		
		fakePhenoSd_=Math.sqrt(fakePhenoSd_*adj);
	}
	private void simulateNormalsRnd(double[] fakePheno) {
		for (int i=0; i<fakePheno.length;i++){
			fakePheno[i]=Settings.jdkRng_.nextGaussian();
			
			calcFakePhenoMeanAndSd();
		}
	}
	
	public static void randomSample(double[] items){
	    for(int i=0;i<items.length;i++){
	        int pos = i + rand.nextInt(items.length - i);
	        double tmp = items[pos];
	        items[pos]=items[i];
	        items[i]= tmp;
	    }	    
	}
	private void addSignal(double[] fakePheno){
				
		for (int i=0; i<fakePheno.length;i++){
			fakePheno[i]=fakePheno[i]+fakeSignal_[i];				
		}		
		NaturalRanking ranking = new NaturalRanking(NaNStrategy.MINIMAL,TiesStrategy.MAXIMUM);
		double[] ranks = ranking.rank(fakePheno);
		for (int i=0; i<fakePheno.length;i++){
			//fakePheno[i]=myNormal.inverseCumulativeProbability(ranks[i]/(ranks.length+1));			
			fakePheno[i]=DistributionMethods.normalInverseCumulativeProbability(ranks[i]/(ranks.length+1));
					
		}			
	}
	
	private void simulateNormals(double[] fakePheno) {
		double p = 0.5/fakePheno.length;
		for (int i=0; i<fakePheno.length;i++){
			fakePheno[i]=DistributionMethods.normalInverseCumulativeProbability(p);			
			p=p+1.0/fakePheno.length;
		}
		randomSample(fakePheno);
	}
	
	private void writeFakePhenotype(){
		for (int i=0; i<fakePheno_.length;i++){
			writerFakePheno.println(String.valueOf(fakePheno_[i]));			
		}
		writerFakePheno.close();
	}
	
	private void runGenescoreMultipleTimesForSameGene(){
		
		
	}
	
	void writeSnpGenotypes(String filename){
		FileExport writerSnpGenotypes = new FileExport(filename);
		int nrOfSnps=snpList_.size();
		int nrOfGenotypes=snpList_.get(0).getGenotypeLength();
		//for (int i=0;i < nrOfGenotypes ; i++){
		String myStr = null;
		for (int i=0;i < nrOfSnps ; i++){
			Snp mySnps=snpList_.get(i);
			byte[] genotypes = mySnps.getGenotypes();
			System.out.println(genotypes.length);
			myStr = "";
			myStr=myStr.concat(mySnps.id_);
			
			for (int j=0;j < nrOfGenotypes ; j++){				
				int currentGenotype = genotypes[j];
				myStr=myStr.concat("\t");
				myStr=myStr.concat(String.valueOf(currentGenotype));				
			}	
			writerSnpGenotypes.println(myStr);
		}
		
		 writerSnpGenotypes.close();
	}
}
