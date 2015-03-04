package ch.unil.genescore.vegas;

import java.util.ArrayList;
import java.util.LinkedList;

import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.main.Settings;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

public class GeneDataProjection extends GeneData {

	
	
	private LinkedList<Snp> currentSnpsWithGenotypes_;
	private SnpWeightPairs SnpsAndWeights_;
	private DenseMatrix weightedCrossCorr_;
	protected double[] snpWiseVar_ = null;
	private ProjectionVegas projectionVegas_;
	protected LdMat weightedCrossCorrMat_;
	private ArrayList<Double> transformedScores_;
	private DenseMatrix matToDecompose_;
	private double snpWiseVarCutoff_;
	private double[] weights_;
	

	public GeneDataProjection(ArrayList<Snp> snpList, LinkedList<Snp> currentSnpsWithGenotypes,SnpWeightPairs SnpsAndWeights) {
		super(snpList);		
		setSnpWiseVarCutoff(0.0);
		currentSnpsWithGenotypes_ = currentSnpsWithGenotypes;
		SnpsAndWeights_ = SnpsAndWeights;
		
        
	}
	/** calculates imputation quality in terms of the variance that is explained away by imputation. (have to reverse weighting the matToDecompose to achieve this) */
	protected void calculateSnpwiseVar(){
		
		double[] weightedVar=MTJConvenienceMethods.getDiagDar(getMatToDecompose());
		double[] rowWeights=weightedCrossCorrMat_.getRowWeights();
		snpWiseVar_ = new double[rowWeights.length];
		for (int i=0;i<rowWeights.length ;i++){
			snpWiseVar_[i]=weightedVar[i]/(rowWeights[i]*rowWeights[i]);
		}
	}
	protected DenseMatrix getMatToDecompose(){
		 return matToDecompose_;
	}
	protected void calculateMatToDecompose(){		
				 DenseMatrix matToDecompose=projectionVegas_.getMatToDecomposeMTJ();
				 matToDecompose_= matToDecompose;
	}
	
	/**outputs filtered version of matToDecompose  (throw out snps with imputation quality at or below snpWiseVarCutoff_)*/
    public DenseMatrix getCorr(){
    	ArrayList<Integer> aboveVarCutoffIndicesAr = new ArrayList<Integer>();
    	
    	for (int i=0; i<snpWiseVar_.length;i++){
    		if(snpWiseVar_[i]>snpWiseVarCutoff_){
    			aboveVarCutoffIndicesAr.add(i);
    		}
    	}    	
    	int[] aboveVarCutoffIndices = new int[aboveVarCutoffIndicesAr.size()];    	
    	for (int i =0 ; i < aboveVarCutoffIndicesAr.size(); i++)
    		aboveVarCutoffIndices[i]=aboveVarCutoffIndicesAr.get(i);
        return MTJConvenienceMethods.getSubDenseMatrixMTJ(getMatToDecompose(),aboveVarCutoffIndices,aboveVarCutoffIndices);
    }
    /**outputs filtered version of transformedScores (throw out snps with imputation quality at or below snpWiseVarCutoff_) */
    public void calculateScores(){
    	DenseVector transformedScoresVect = projectionVegas_.getTransformedZscores();
    	ArrayList<Double> transformedScores = new ArrayList<Double>();
    	for (int i=0 ; i < transformedScoresVect.size() ; i++){
    		if (snpWiseVar_[i]>snpWiseVarCutoff_)
    		transformedScores.add(transformedScoresVect.get(i));
    				
    	}
    	transformedScores_= transformedScores;
    }
    
    public void calculateWeights(){
    	double[] weights = weightedCrossCorrMat_.getRowWeights();
    	ArrayList<Double> weightAr = new ArrayList<Double>();
    	for (int i=0 ; i < weights.length ; i++){
    		if (snpWiseVar_[i]>snpWiseVarCutoff_)
    		weightAr.add(weights[i]);    				
    	}
    	weights_ = ConvenienceMethods.arListToDoubleAr(weightAr);    	
    }
   
    
    public void writeGeneToFile(String fileName){
    	
    	ArrayList<Snp> allSnps =  SnpsAndWeights_.getSnps();
    	ArrayList<Double> allWeights = SnpsAndWeights_.getWeights();
    	ArrayList<Double> varExplainedPerSnp = MTJConvenienceMethods.getDiag(getMatToDecompose());
       	ArrayList<String> valStrings = new ArrayList<String>(); 
   
       	for (int i=0 ; i < allSnps.size() ; i++){
    		Snp currentSnp=allSnps.get(i);
    		int snpInGwas;
    		if (snpList_.contains(currentSnp)){
    			snpInGwas=1;
    		}
    		else{
    			snpInGwas=0;
    		}
    		   		   		
    		String currentString = "";
    		currentString += String.valueOf(transformedScores_.get(i)) + "\t"; 
    		currentString += String.valueOf(snpInGwas) + "\t";
    		currentString += String.valueOf(allWeights.get(i)) + "\t"; 
    		currentString += String.valueOf(varExplainedPerSnp.get(i)); 
    		
    		valStrings.add(currentString);
    		
    	}
    	String header="scores\tinGwas\tweighting\tvarExplained";
    	
    	WritingMethods.writeSnpPosWithValToFile(allSnps, valStrings,header, fileName,Settings.writeGenewiseSnpFiles_);
    	String corfileName = "corMat_" + fileName;
    	
    	//	WritingMethods.writeMTJ(projectionVegas_.getMatToDecomposeMTJ(), corfileName,Settings.writeGenewiseSnpFiles_);
    	
        
    }	
    public void processData(){
    	super.processData();
    	ArrayList<Snp> SnpList =  SnpsAndWeights_.getSnps();
    	insertWeightedCrossCorrMat();        
        weightedCrossCorr_ = weightedCrossCorrMat_.getWeighedCrossCorrelationMatrix();
        // in case of normalized genotypes this is proportianal to the prior variance of the vegas statistic  
    	// because the trace of a matrix is equal to the sum of eigenvalues.        	 
        double priorVar = weightedCrossCorrMat_.getSumOfWeightsSquared(true);
        ArrayList<Double> snpScores = super.getScores();
        insertProjectionVegas(snpScores, priorVar);      
        calculateSnpwiseVar();
        calculateScores();
        calculateWeights();
        
    }
    /*
    public DenseMatrix pruneMat(UpperSymmDenseMatrix correl){
		double pruningCutoff=0.99;		
		int n=correl.numColumns();		
		boolean[] indices = new boolean[n];		
		for (int i=0;i<n ; i++){
			for (int j=(i+1);j<n; j++){
				if(correl.get(i,j) <= pruningCutoff){
					indices[j]=true;
					double currentMax=Math.max(upper_[i], upper_[j]);
					upper_[i]=currentMax;
					upper_[j]=currentMax;
				}
			}			
		}
		int counter=0;
		for (boolean index : indices){
			if (index)
				counter++;
		}
		int[] keepIndices = new int[counter];		
		int counter2=0;
		for (int i=0; i< n; i++){
			
			if (indices[i]){				
				keepIndices[counter2]=i;
				counter2++;
			}
		}
		double[] upperPruned = new double[counter];
		for (int i=0; i < counter; i++){
				upperPruned[i]=upper_[keepIndices[i]];
		}
		DenseMatrix correlOut=MTJConvenienceMethods.getSubDenseMatrixMTJ(MTJConvenienceMethods.DenseFromSymmetric(correl), keepIndices, keepIndices);
		return correlOut;
	}

    */
    
    protected void insertWeightedCrossCorrMat(){
    	LdMat myLdMat = new LdMat(SnpsAndWeights_,snpList_);
    	myLdMat.processData();
    	weightedCrossCorrMat_ = myLdMat;
    }
    protected void insertProjectionVegas(ArrayList<Double> snpScores, double priorVar){
    	projectionVegas_ = new ProjectionVegas(snpScores, cov_, weightedCrossCorr_, priorVar);    	
    	projectionVegas_.processData();
    	calculateMatToDecompose();
    }
    public double[] returnWeights(boolean DownWeighWithImputationQuality){
    	double[] weights = weights_;
    	if (DownWeighWithImputationQuality){
    		for (int i=0; i<weights_.length ;i++){
    			weights[i]=weights[i]*Math.sqrt(snpWiseVar_[i]);
    		}
    	}
       	return weights;
     }    
    public ArrayList<Double> getScores(){
    	return transformedScores_;}
    public void setSnpWiseVarCutoff(double cutoff){snpWiseVarCutoff_=cutoff;}    
    public double[] getSnpWiseVar(){return snpWiseVar_;}
    
    
    
 
    		
}

	

