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
package ch.unil.genescore.vegas.test;

import static org.junit.Assert.*;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

import org.junit.*;

import ch.unil.genescore.main.*;
import ch.unil.genescore.vegas.MTJConvenienceMethods;


/**
 * Unit tests for AnalyticVegas.
 * TODO @David As inspiration... :)
 */
public class MTJConvenienceMethodsTest {
	static DenseMatrix matA = null;
	static DenseMatrix matB = null;
	static DenseMatrix matC = null;
	static DenseMatrix matZero = null;
	static double delta;
	
	// ============================================================================
		// SETUP
		
		@BeforeClass
		public  static void testSetup() {
			Settings.loadSettings();
			delta = 1E-14;
			matA = new DenseMatrix(4,4);
			for (int i=0; i < 4 ; i++){
				for (int j=0; j < 4 ; j++){
					matA.set(i, j, (double) (i*j));		
				}
			}
			matB = new DenseMatrix(4,4);
			for (int i=0; i < 4 ; i++){
				for (int j=0; j < 4 ; j++){
					matB.set(i, j, (double) (i-j));		
				}
			}	
			matC = new DenseMatrix(4,6);
			for (int i=0; i < 4 ; i++){
				for (int j=0; j < 6 ; j++){
					matC.set(i, j, (double) (i-j));		
				}
			}
			matZero = new DenseMatrix(4,4);
			for (int i=0; i < 4 ; i++){
				for (int j=0; j < 4 ; j++){
					matZero.set(i, j, 0);		
				}
			}
		}
		

		@AfterClass
		public static void testCleanup() { }
		
		
		// ============================================================================
		// TESTS
		
		/** Test HadamardProduct */
		@Test
		public void testHadamardProduct() {
			DenseMatrix matTarget = new DenseMatrix(4,4);
			matTarget.set(1,2,-2);matTarget.set(1,3,-6);
			matTarget.set(2,1,2);matTarget.set(2,3,-6);
			matTarget.set(3,1,6);matTarget.set(3,2,6);			
			DenseMatrix matResult = null;
			matResult = MTJConvenienceMethods.hadamardProduct(matA,matB);
			double[] dataTarget = matTarget.getData();
			double[] dataResult = matResult.getData();
			
			for (int i=0; i < dataTarget.length;i++)
				assertEquals(dataTarget[i],dataResult[i], delta);			
		}
		/**Test regularizeMat*/
		@Test
		public void testRegularizeMat(){
			DenseMatrix matTarget = new DenseMatrix(4,4);
			for (int i=0; i < 4 ; i++){
				for (int j=0; j < 4 ; j++){
					if (i==j)
						matTarget.set(i,j,0.5);
					else 
						matTarget.set(i,j,0);							
				}
			}
			DenseMatrix matResult = MTJConvenienceMethods.regularizeMat(matZero, 0.5);
			UpperSymmDenseMatrix matResultUpperSymm = MTJConvenienceMethods.regularizeMat(new UpperSymmDenseMatrix(matZero), 0.5);
			double[] dataTarget = matTarget.getData();
			double[] dataResult = matResult.getData();
			double[] dataResultUpperSymm = matResultUpperSymm.getData();
			
			for (int i=0; i < dataTarget.length;i++){
				assertEquals(dataTarget[i],dataResult[i],delta);
				assertEquals(dataTarget[i],dataResultUpperSymm[i], delta);
			
			}
		
		}
		/**Test getSubDenseMat*/
		@Test
		public void testGetSubDenseMatrixMTJ(){
			//test general version
			DenseMatrix matTarget = new DenseMatrix(2,2);
			matTarget.set(0,0,1);matTarget.set(0,1,2);
			matTarget.set(1,0,3);matTarget.set(1,1,6);
			int[] rowIndices = {1,3};
			int[] colIndices = {1,2};
			DenseMatrix matResult = MTJConvenienceMethods.getSubDenseMatrixMTJ(matA, rowIndices,colIndices);
			double[] dataTarget = matTarget.getData();
			double[] dataResult = matResult.getData();
			for (int i=0; i < dataTarget.length;i++)
				assertEquals(dataTarget[i],dataResult[i], delta);
			//test version with only one index/
			matTarget.set(0,0,4);matTarget.set(0,1,6);
			matTarget.set(1,0,6);matTarget.set(1,1,9);
			matResult = MTJConvenienceMethods.getSubDenseMatrixMTJ(matA, 2,2);
			dataTarget = matTarget.getData();
			dataResult = matResult.getData();
			for (int i=0; i < dataTarget.length;i++)
				assertEquals(dataTarget[i],dataResult[i],delta);			
		}
		/**Test doGammaLambdaGammaTMTJ*/
		@Test
		 public void testDoGammaLambdaGammaTMTJ(){
			 // exploit special structure of the defined mats.
			 DenseMatrix matResult= MTJConvenienceMethods.doGammaLambdaGammaTMTJ(matA, matB);			
				double[] dataResult = matResult.getData();
				for (int i=0; i < dataResult.length;i++)
					assertEquals(dataResult[i],0,delta);
		 }
		/**Test repMTJ*/
		@Test
		 public void testRepMTJ(){			 
			//test row direction
			DenseMatrix matTarget = new DenseMatrix(2,2);
			matTarget.set(0,0,3);matTarget.set(0,1,6);
			matTarget.set(1,0,3);matTarget.set(1,1,6);			
			double[] vecToRep = {3,6};
			 DenseMatrix matResult= MTJConvenienceMethods.repMTJ(vecToRep,2,true);			
			 double[] dataResult = matResult.getData();
			 double[] dataTarget = matTarget.getData();
			for (int i=0; i < dataResult.length;i++)
				assertEquals(dataResult[i],dataTarget[i],delta);
		 	
			//test row direction		
		matTarget.set(0,0,3);matTarget.set(0,1,3);
		matTarget.set(1,0,6);matTarget.set(1,1,6);			
		matResult= MTJConvenienceMethods.repMTJ(vecToRep,2,false);			
		 dataResult = matResult.getData();
		 dataTarget = matTarget.getData();
		for (int i=0; i < dataResult.length;i++)
			assertEquals(dataResult[i],dataTarget[i],delta);
	 	}
		/**Test diagMTJ*/
		@Test
		public void testDiagMTJ(){
		//		double[] diag)
			DenseMatrix matTarget = new DenseMatrix(3,3);
			double[] myDiag = {0.4, 1, 2};
			for (int i=0 ; i<3 ; i++)
				matTarget.set(i,i,myDiag[i]);
			DenseMatrix matResult= MTJConvenienceMethods.diagMTJ(myDiag);
			double[] dataResult = matResult.getData();
			double[] dataTarget = matTarget.getData();
			for (int i=0; i < dataResult.length;i++)
				assertEquals(dataResult[i],dataTarget[i],delta);
			
		}
		@Test
		public void testScaleToDiagOne(){
			DenseMatrix ld= new DenseMatrix(3,3);
			DenseMatrix ldTarget= new DenseMatrix(3,3);
			int n=3;
			ld.set(0,0,1);ld.set(1,1,1);ld.set(2,2,1);
			ld.set(0,1,0.9);ld.set(1,2,0.9);
			ld.set(1,0,0.9);ld.set(2,1,0.9);
			ld.set(0,2,0.81);
			ld.set(2,0,0.81);
			
			ldTarget.set(0,0,1);ldTarget.set(1,1,0.81);ldTarget.set(2,2,1);
			ldTarget.set(0,1,0.81);ldTarget.set(1,2,0.81);
			ldTarget.set(1,0,0.81);ldTarget.set(2,1,0.81);
			ldTarget.set(0,2,0.81);
			ldTarget.set(2,0,0.81);
			DenseMatrix myLd=MTJConvenienceMethods.scaleToDiagOne(ldTarget);
			for (int i=0; i<n; i++){
				for (int j=0; j<n; j++){
					assertEquals(myLd.get(j, i),ld.get(j,i), delta);
					
				}
			}
			UpperSymmDenseMatrix myLd2=MTJConvenienceMethods.scaleToDiagOne(new UpperSymmDenseMatrix(ldTarget));
			for (int i=0; i<n; i++){
				for (int j=0; j<n; j++){
					assertEquals(myLd2.get(j, i),ld.get(j,i), delta);
					
				}
			}
			
		}
		
		@Test
		public void testDeScale(){
			DenseMatrix ld= new DenseMatrix(3,3);
			DenseMatrix ldTarget= new DenseMatrix(3,3);
			int n=3;
			ld.set(0,0,1);ld.set(1,1,1);ld.set(2,2,1);
			ld.set(0,1,0.9);ld.set(1,2,0.9);
			ld.set(1,0,0.9);ld.set(2,1,0.9);
			ld.set(0,2,0.81);
			ld.set(2,0,0.81);
			
			ldTarget.set(0,0,1);ldTarget.set(1,1,0.81);ldTarget.set(2,2,1);
			ldTarget.set(0,1,0.81);ldTarget.set(1,2,0.81);
			ldTarget.set(1,0,0.81);ldTarget.set(2,1,0.81);
			ldTarget.set(0,2,0.81);
			ldTarget.set(2,0,0.81);
			
			double[] weights={1,0.9,1};
			UpperSymmDenseMatrix myLd2=MTJConvenienceMethods.deScale(new UpperSymmDenseMatrix(ldTarget), weights);
			for (int i=0; i<n; i++){
				for (int j=0; j<n; j++){
					assertEquals(myLd2.get(j, i),ld.get(j,i), delta);
					
				}
			}
			
		}
		
		/**Test calculateCovarianceMat*/
		@Test
		public void testCalculateCovarianceMat(){
			
			DenseMatrix dat = new DenseMatrix(8,2);					
			for (int i=0; i < 8 ; i++){
				for (int j=0; j < 2 ; j++){
					dat.set(i, j, ((i-j)*i));		
				}
			}
			
			DenseMatrix testMat = MTJConvenienceMethods.calculateCovarianceMat(dat);
			assertEquals(testMat.get(0, 0),278.25,delta);			
			assertEquals(testMat.get(1, 0),241.5,delta);
			assertEquals(testMat.get(0, 1),241.5,delta);
			assertEquals(testMat.get(1, 1),210.0,delta);				
			
		}
		
		/**Test calculateCrossCovarianceMat*/
		@Test
		public void testCalculateCrossCovarianceMat(){
		//		double[] diag)
			DenseMatrix datX = new DenseMatrix(8,2);					
			for (int i=0; i < 8 ; i++){
				for (int j=0; j < 2 ; j++){
					datX.set(i, j, ((i-j)*i));		
				}
			}
			DenseMatrix datY = new DenseMatrix(8,3);					
			for (int i=0; i < 8 ; i++){
				for (int j=0; j < 3 ; j++){
					datY.set(i, j,((i-j)*(i*i)));		
				}
			}
			DenseMatrix testMat = MTJConvenienceMethods.calculateCrossCovarianceMat(datX,datY);
			assertEquals(testMat.get(0,0),1911.0,delta);
			assertEquals(testMat.get(1,0),1669.5,delta);
			assertEquals(testMat.get(0,1),1632.75,delta);
			assertEquals(testMat.get(1,1),1428,delta);
			assertEquals(testMat.get(0,2),1354.5,delta);
			assertEquals(testMat.get(1,2),1186.5,delta);			
		}
		
		
}
