package ch.unil.genescore.vegas;

import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

public class MTJConvenienceMethods {
	
	static public DenseMatrix hadamardProduct(Matrix A, Matrix B){
		assert(A.numColumns()==B.numColumns());
		assert(A.numRows()==B.numRows());
		DenseMatrix matOut = new DenseMatrix(A.numRows(),A.numColumns());
		double val;
		for (int i = 0; i < matOut.numRows();i++){			
			for (int j = 0; j < matOut.numColumns();j++){
				val=A.get(i,j)*B.get(i,j);
				matOut.set(i,j, val);
			}
		}
		return matOut;
	}
	
	static public DenseMatrix repMTJ(double[] vect, int nrToRep, boolean repOverRows){
		double val;
		DenseMatrix matOut = new DenseMatrix(nrToRep,vect.length);
		for (int i = 0; i < matOut.numRows();i++){
			for (int j = 0; j < matOut.numColumns();j++){
				val=vect[j];
				matOut.set(i,j, val);
			}
		}
		if (repOverRows==false){			
			DenseMatrix matOutT = new DenseMatrix(matOut.numColumns(),matOut.numRows());
			matOut.transpose(matOutT);
			return matOutT;
		}
		return matOut;
	}

	static public DenseMatrix diagMTJ(double[] diag){
		
		DenseMatrix MatOut = new DenseMatrix(diag.length, diag.length);
		for (int i = 0;i < diag.length; ++i)			
			MatOut.set(i,i, diag[i]);
		return MatOut;
	}
	
	static public double traceMTJ(UpperSymmDenseMatrix matToDecomposeMTJ_){
		if (matToDecomposeMTJ_.numColumns()!=matToDecomposeMTJ_.numRows()){
			throw new RuntimeException("matrix not square, cannot calculate trace");
		}
		double out=0;			
		for (int i = 0;i < matToDecomposeMTJ_.numRows(); ++i)			
			out += matToDecomposeMTJ_.get(i,i);
		return out;
	}
	static public double traceMTJ(DenseMatrix matToDecomposeMTJ_){
		if (matToDecomposeMTJ_.numColumns()!=matToDecomposeMTJ_.numRows()){
			throw new RuntimeException("matrix not square, cannot calculate trace");
		}
		double out=0;			
		for (int i = 0;i < matToDecomposeMTJ_.numRows(); ++i)			
			out += matToDecomposeMTJ_.get(i,i);
		return out;
	}

	static public DenseMatrix regularizeMat(DenseMatrix Mat, double eps){
		DenseMatrix matOut= new DenseMatrix(Mat);		
		assert(matOut.numColumns()==matOut.numRows());
		double[] ar = new double[matOut.numRows()];
		for (int i=0 ; i<matOut.numRows(); i++)
			ar[i]=eps;
		Matrix myDiag = diagMTJ(ar);
		matOut.add(myDiag);
		return matOut;
	}
static public UpperSymmDenseMatrix regularizeMat(UpperSymmDenseMatrix Mat, double eps){
		
		assert(Mat.numColumns()==Mat.numRows());
		double[] ar = new double[Mat.numRows()];
		for (int i=0 ; i<Mat.numRows(); i++)
			ar[i]=eps;
		Matrix myDiag = diagMTJ(ar);
		Mat.add(myDiag);
		return Mat;
	}
	public static DenseMatrix getSubDenseMatrixMTJ(DenseMatrix DenseMat,int[] rowIndices,int[] colIndices){
		DenseMatrix DenseMatOut = new DenseMatrix(rowIndices.length,colIndices.length);  
		int col=0;
		int row=0;
		for (int i=0;i<rowIndices.length;i++){
			row=rowIndices[i];
			for (int j=0;j<colIndices.length;j++){
				col=colIndices[j];			
				DenseMatOut.set(i,j,DenseMat.get(row,col));
			}			
		}
		return DenseMatOut;
	}
	
	public static UpperSymmDenseMatrix getSubDenseMatrixSymmMTJ(UpperSymmDenseMatrix SymmDenseMat,int[] indices){
		UpperSymmDenseMatrix SymmDenseMatOut = new UpperSymmDenseMatrix(indices.length);  
		int col=0;
		int row=0;
		for (int i=0;i<indices.length;i++){
			row=indices[i];
			for (int j=0;j<indices.length;j++){
				col=indices[j];			
				SymmDenseMatOut.set(i,j,SymmDenseMat.get(row,col));
			}			
		}
		return SymmDenseMatOut;
	}			
	public static DenseMatrix denseFromSymmetric(UpperSymmDenseMatrix mat){		
		DenseMatrix outMat = new DenseMatrix(mat.numColumns(),mat.numColumns());
		outMat.set(mat);
		return outMat;
	}
	public static UpperSymmDenseMatrix symmetricFromDense(DenseMatrix mat){		
		UpperSymmDenseMatrix outMat = new UpperSymmDenseMatrix(mat.numColumns());
		outMat.set(mat);
		return outMat;
	}
	/** create subMatrix from a given matrix. keeps all rows with index greater or equal to rowMinIndex. dito for colMinIndex.  */
	public static DenseMatrix getSubDenseMatrixMTJ(DenseMatrix DenseMat,int rowMinIndex,int colMinIndex){
		int numRows = DenseMat.numRows();
		int numColumns = DenseMat.numColumns();		
		int[] colsSelected = new int[numColumns-colMinIndex];
		int[] rowsSelected = new int[numRows-rowMinIndex];
		for (int i=0;i < rowsSelected.length ; i++)
			rowsSelected[i]=rowMinIndex+i;
		for (int i=0;i < colsSelected.length ; i++)
			colsSelected[i]=colMinIndex+i;
		DenseMatrix DenseMatOut = getSubDenseMatrixMTJ(DenseMat, rowsSelected, colsSelected);
		return DenseMatOut;
	}
	
	/**calculates \Gamma \Lambda \Gamma^T */
	public static DenseMatrix doGammaLambdaGammaTMTJ(DenseMatrix Gamma, DenseMatrix Lambda){
	// MTJ typical notation. method without argument manipulates matrix in place
	// method with extra matrix argument copies manipulation result into last argument.
		//see e.g. Gamma.transpose(GammaT);
		DenseMatrix interim = new DenseMatrix(Gamma.numRows(), Lambda.numColumns());
		DenseMatrix outMat = new DenseMatrix(Gamma.numRows(),Gamma.numRows());
		DenseMatrix  GammaT = new DenseMatrix(Gamma.numColumns(),Gamma.numRows());
		Gamma.mult(Lambda, interim);
		Gamma.transpose(GammaT); 
		interim.mult(GammaT, outMat);
		return outMat;
	}
	/**calculates \Gamma^T \Lambda \Gamma */
	public static DenseMatrix doGammaTLambdaGammaMTJ(DenseMatrix Gamma, DenseMatrix Lambda){
	// MTJ typical notation. method without argument manipulates matrix in place
	// method with extra matrix argument copies manipulation result into last argument.
		//see e.g. Gamma.transpose(GammaT);
		DenseMatrix interim = new DenseMatrix(Gamma.numRows(), Lambda.numColumns());
		DenseMatrix outMat = new DenseMatrix(Gamma.numRows(),Gamma.numRows());
		DenseMatrix  GammaT = new DenseMatrix(Gamma.numColumns(),Gamma.numRows());
		Gamma.transpose(GammaT); 
		GammaT.mult(Lambda, interim);		
		interim.mult(Gamma, outMat);
		return outMat;
	}
	public static DenseMatrix demean(DenseMatrix X){
		DenseMatrix matOut = X.copy();  				
		DenseMatrix myOnes = ones(1, X.numRows());
		myOnes.scale(1.0/myOnes.numColumns());
		DenseMatrix myScales = new DenseMatrix(1,X.numColumns());
		myOnes.mult(matOut, myScales);
		DenseMatrix meanMat= repMTJ(myScales.getData(),X.numRows(), true);		
		matOut.add(-1, meanMat);
		return matOut; 
	}
	public static DenseMatrix scaleToDiagOne(DenseMatrix X){
		DenseMatrix matOut = null;
		
		double[] myScaler=getDiagDar(X);
		for (int i=0; i< myScaler.length;i++)
			myScaler[i]=1.0/Math.sqrt(myScaler[i]);
		DenseMatrix scalerMat = hadamardProduct(repMTJ(myScaler, myScaler.length, false), repMTJ(myScaler, myScaler.length, true));		
		matOut= hadamardProduct(scalerMat,X);
		return matOut;
	}
	public static UpperSymmDenseMatrix scaleToDiagOne(UpperSymmDenseMatrix X){
		UpperSymmDenseMatrix matOut = null;
		
		double[] myScaler=getDiagDar(X);
		for (int i=0; i< myScaler.length;i++)
			myScaler[i]=1.0/Math.sqrt(myScaler[i]);
		DenseMatrix scalerMat = hadamardProduct(repMTJ(myScaler, myScaler.length, false), repMTJ(myScaler, myScaler.length, true));		
		matOut = new UpperSymmDenseMatrix(hadamardProduct(scalerMat,X));
		return matOut;
	}
	
	
	/**descale i.e. divide row and columnwise with myScaler entries (such that diagonal-entry (i,i) is divided by my myScaler[i]^2) */
	public static UpperSymmDenseMatrix deScale(UpperSymmDenseMatrix X, double[] myScaler){
		double[] myScalerInv = new double[myScaler.length];
		UpperSymmDenseMatrix matOut = null;
		for (int i=0; i< myScaler.length;i++)
			myScalerInv[i]=1.0/(myScaler[i]);
		DenseMatrix scalerMat = hadamardProduct(repMTJ(myScalerInv, myScalerInv.length, false), repMTJ(myScalerInv, myScalerInv.length, true));		
		matOut = new UpperSymmDenseMatrix(hadamardProduct(scalerMat,X));
		return matOut;
	}
	
	

	
	public static DenseMatrix ones(int numRows, int numColumns){	
		DenseMatrix ones = new DenseMatrix(numRows,numColumns);
		for (int  i=0; i < numRows; i++){
			for (int  j=0; j < numColumns; j++)
				ones.set(i, j, 1);
		}	
		return ones;					
	}
	/**calculates covariance matrix */
	public static DenseMatrix calculateCovarianceMat(DenseMatrix X){
	// MTJ typical notation. method without argument manipulates matrix in place
	// method with extra matrix argument copies manipulation result into last argument.
		//see e.g. Gamma.transpose(GammaT);
		DenseMatrix Xdemeaned=demean(X);
		DenseMatrix  covMat = new 	DenseMatrix (Xdemeaned.numColumns(),Xdemeaned.numColumns());
		double invNumberOfRows = 1.0/((double) X.numRows());
		Xdemeaned.transAmult(invNumberOfRows, Xdemeaned, covMat);
		return covMat;
		
	}
	/**returns diag as list*/
	public static ArrayList<Double> getDiag(DenseMatrix X){
		
		assert(X.numColumns()==X.numRows());
		ArrayList<Double> ar = new ArrayList<Double>();
		for (int i=0 ; i<X.numRows(); i++)
			ar.add(X.get(i, i));
		return ar;
		
	}
	/**returns diag as list*/
	public static double[] getDiagDar(Matrix X){
		
		assert(X.numColumns()==X.numRows());
		double[] ar = new double[X.numColumns()];
		for (int i=0 ; i<X.numRows(); i++)
			ar[i]=X.get(i, i);
		return ar;
		
	}
	public static DenseMatrix arToMat(double[] input){
		DenseMatrix mat = new DenseMatrix(input.length, 1);
		for (int i=0 ; i< input.length; i++)
			mat.set(i, 0,input[i]);
		return mat;
	}
	
	/**calculates cross-covariance matrix */
	public static DenseMatrix calculateCrossCovarianceMat(DenseMatrix X,DenseMatrix Y){
	// MTJ typical notation. method without argument manipulates matrix in place
	// method with extra matrix argument copies manipulation result into last argument.
		//see e.g. Gamma.transpose(GammaT);
			if (X.numRows()!=Y.numRows()){
				throw new RuntimeException("error: not same number of observations in both matrices.");
			}
			DenseMatrix Xdemeaned=demean(X);
			DenseMatrix Ydemeaned=demean(Y);
					
			DenseMatrix  crossCovMat = new DenseMatrix(Xdemeaned.numColumns(),Ydemeaned.numColumns());
			double invNumberOfRows = 1.0/((double) Xdemeaned.numRows());
			Xdemeaned.transAmult(invNumberOfRows, Ydemeaned, crossCovMat);
			return crossCovMat;
		
	}
	
	public static DenseMatrix invertMatrix(DenseMatrix A){
		
		DenseMatrix I = Matrices.identity(A.numColumns());
		DenseMatrix AI = I.copy();
		A.solve(I, AI);	
		return(AI);
	}
	
public static DenseMatrix regularizeScaleInvertMat(DenseMatrix A, double eps){
		DenseMatrix AI = invertMatrix(scaleToDiagOne(regularizeMat(A,eps)));		
		return(AI);
	}		
}
