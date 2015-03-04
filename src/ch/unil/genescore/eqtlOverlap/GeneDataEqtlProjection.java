package ch.unil.genescore.eqtlOverlap;

import java.util.ArrayList;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import ch.unil.genescore.main.ConvenienceMethods;
import ch.unil.genescore.vegas.GeneData;
import ch.unil.genescore.vegas.GeneDataInterface;
import ch.unil.genescore.vegas.LinkageDisequilibrium;
import ch.unil.genescore.vegas.MTJConvenienceMethods;
import ch.unil.genescore.vegas.Snp;

public class GeneDataEqtlProjection implements GeneDataInterface{

	//public GeneDataEqtlProjection(ArrayList<Snp> snpList) {
		//super(snpList);
		// TODO Auto-generated constructor stub
	//}


	private Set<Snp> geneSnps_ = null;
	private Set<Snp> eqtlSnps_ =null;
	private Set<Snp> eqtlNeighbourSnps_ = null;
	private DenseMatrix cov_ =null;
	private ArrayList snpScores_ =null;
	

	public  GeneDataEqtlProjection(Set<Snp> geneSnps,  Set<Snp> eqtlSnps, Set<Snp> eqtlNeighbourSnps){
		
		geneSnps_=geneSnps;
		eqtlSnps_=eqtlSnps;
		eqtlNeighbourSnps_=eqtlNeighbourSnps;
		
	}
	private ArrayList<Snp> hashToAr(Set<Snp> H){
		ArrayList<Snp> newList = new ArrayList<Snp>();
		for (Snp h : H){
			newList.add(h);
		}
		return(newList);
	}	
	
	
	//private void makeMatrices(){
	@Override
	public void processData(){
		eqtlNeighbourSnps_.removeAll(eqtlSnps_);
		geneSnps_.addAll(eqtlSnps_);
		geneSnps_.removeAll(eqtlNeighbourSnps_);
		
		ArrayList<Snp> gSnps = hashToAr(geneSnps_);
		ArrayList<Snp> eSnps = hashToAr(eqtlSnps_);
		//TODO: remove again
		//ArrayList<Snp> eSnps = new ArrayList<Snp>();
		ArrayList<Snp> enSnps = hashToAr(eqtlNeighbourSnps_);
		
		int numG= gSnps.size();
		int numE= eSnps.size();
		int numEN= enSnps.size();
		
		DenseMatrix Z_G = LinkageDisequilibrium.getZscores(gSnps);
		DenseMatrix Z_E = LinkageDisequilibrium.getZscores(eSnps);
		DenseMatrix Z_EN = LinkageDisequilibrium.getZscores(enSnps);
		
		DenseMatrix Sig_EN_E = LinkageDisequilibrium.computeCrossCorrelationMatrixMTJ(enSnps, eSnps);		
		DenseMatrix Sig_G_E = LinkageDisequilibrium.computeCrossCorrelationMatrixMTJ(gSnps, eSnps);		
		DenseMatrix Sig_G_EN = LinkageDisequilibrium.computeCrossCorrelationMatrixMTJ(gSnps, enSnps);	
		
		DenseMatrix Sig_E = LinkageDisequilibrium.computeCorrelationMatrixMTJ(eSnps);
		DenseMatrix Sig_EN = LinkageDisequilibrium.computeCorrelationMatrixMTJ(enSnps);
		DenseMatrix Sig_G = LinkageDisequilibrium.computeCorrelationMatrixMTJ(gSnps);
		
		if(numEN==0){
			System.out.println("no eqtls for other genes in neighbourhood");
			cov_=Sig_G;
			snpScores_ = ConvenienceMethods.doubleArToArList(Z_G.getData());	
			return;
		}
		
		//TODO: regularize abit:
		double eps=0.05;
		DenseMatrix InvSig_E = MTJConvenienceMethods.regularizeScaleInvertMat(Sig_E,eps);

		//calculate Sig_ENdivE
		DenseMatrix  interim = null;
		DenseMatrix  s_interim = null;
		DenseMatrix  z_interim = null;
		
		interim = MTJConvenienceMethods.doGammaLambdaGammaTMTJ(Sig_EN_E,InvSig_E);				
		DenseMatrix Sig_ENdivE=Sig_EN.copy();
		Sig_ENdivE.add(-1,interim);		
		DenseMatrix InvSig_ENdivE=MTJConvenienceMethods.regularizeScaleInvertMat(Sig_ENdivE,eps);
		
	
		//calculate Sig_G_ENdivE;
		interim= new DenseMatrix(numG,numE);
		Sig_G_E.mult(InvSig_E, interim);
		DenseMatrix Sig_G_ENdivE = Sig_G_EN.copy();	
		s_interim = new DenseMatrix(numG,numEN);
		interim.transBmult(Sig_EN_E, s_interim);
		Sig_G_ENdivE.add(-1, s_interim);		
		
		//calculate Z_ENdivE
		interim=new DenseMatrix(numEN, numE);
		Sig_EN_E.mult(InvSig_E, interim);
		z_interim= new DenseMatrix(numEN,1);
		interim.mult(Z_E,z_interim);		
		DenseMatrix Z_ENdivE = Z_EN.copy();
		Z_ENdivE.add(-1,z_interim);
		
		
		DenseMatrix Z_GdivENdivE = Z_G.copy();
		z_interim = new DenseMatrix(numG, 1);		
		interim=new DenseMatrix(numG,numEN);
		Sig_G_ENdivE.mult(InvSig_ENdivE, interim);		
		interim.mult(Z_ENdivE, z_interim);
		Z_GdivENdivE.add(-1,z_interim);
		
		DenseMatrix Sig_GdivENdivE = Sig_G.copy();
		
		s_interim = MTJConvenienceMethods.doGammaLambdaGammaTMTJ(Sig_G_ENdivE,InvSig_ENdivE); 
		Sig_GdivENdivE.add(-1,s_interim);
		cov_ = Sig_GdivENdivE; 
		snpScores_ = ConvenienceMethods.doubleArToArList(Z_GdivENdivE.getData());
	}
	
	@Override
	public DenseMatrix getCorr(){return cov_;}
	
	@Override	
	public ArrayList<Double> getScores() {		
		return snpScores_;
	}
	@Override
	public double[] getWeights() {
		double[] d = new double[snpScores_.size()];
		for(int i=0; i<snpScores_.size();i++){
			d[i]=1;
		}
		return d;
	}	 
	

}
