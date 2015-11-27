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
package ch.unil.genescore.vegas;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Pascal;
import ch.unil.genescore.main.Settings;

public class WritingMethods {

	public static void writeMTJ(Matrix mat, String fileName, String additonalDirectory) {
		
		String filename = Settings.outputDirectory_ + "/" + additonalDirectory + "/" + fileName;
		System.out.println(filename); 
		File file = new File(filename);
	      try {
			file.createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
		FileExport writer = new FileExport(filename);
		
		for (int i = 0; i < mat.numRows(); i++){
			String line = "";
			for (int j = 0; j < mat.numColumns(); j++){
				if (j == 0){
					line += mat.get(i,j);
				}
				else {
					line += "\t" + mat.get(i,j);
				}				
			}
			writer.println(line);
		}
		writer.close();
		Pascal.println("");
	}

	public static void writeGenotypeToFile(ArrayList<Snp> snpList, String fileName, String additonalDirectory) {
	
		int numRows = snpList.get(0).getGenotypes().length;
		int numColumns = snpList.size();
		DenseMatrix printMat = new DenseMatrix(numRows, numColumns);		
		double val=0;
		for (int j=0; j<snpList.size();j++){			
			byte[] genotype = snpList.get(j).getGenotypes();
			for (int i=0 ; i<genotype.length ; i++){
				val=genotype[i];	
				printMat.set(i,j,val);				
			}		
		}
		writeMTJ( printMat, fileName, additonalDirectory);
	}
	public static void writeZscoreToFile(ArrayList<Snp> snpList, String fileName, String additonalDirectory) {
		
		
		int numSnps = snpList.size();
		DenseMatrix printMat = new DenseMatrix(numSnps,1);		
		double val=0;					
		for (int i=0 ; i< numSnps ; i++){
			val=snpList.get(i).getZscore();	
			printMat.set(i,0,val);				
		}				
		writeMTJ( printMat, fileName, additonalDirectory);
	}
	
	public static void writeSnpPosWithValToFile(ArrayList<Snp> snpList, ArrayList<String> valString,String header, String fileName, String additonalDirectory) {
		
		int numSnps= snpList.size();
		if (numSnps!=valString.size()){
			throw new RuntimeException("can't write; not same length");
		}
		String filename = Settings.outputDirectory_ + "/" + additonalDirectory + "/" + fileName;
		System.out.println(filename); 
		File file = new File(filename);
	      try {
			file.createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
	    
		FileExport writer = new FileExport(filename);
		String headerline= "chr\tstart\tend\tid\t" + header;
		writer.println(headerline);
		for (int i=0 ; i < numSnps ; i++){
			String line = "";
			Snp curSnp = snpList.get(i);
			String curVal = valString.get(i);
			
			line += curSnp.chr_;
			line += "\t" + curSnp.start_;
			line += "\t" + curSnp.end_;
			line += "\t" + curSnp.id_;
		
			line += "\t" + curVal;
			writer.println(line);
		}
		writer.close();
	}
	
	// ----------------------------------------------------------------------------

	/** Write correlation matrix of a gene to file */
	public static void writeLD(double[][] cor, String geneName, String additonalDirectory) {
						
		String filename = Settings.outputDirectory_ + "/" + additonalDirectory + "/" + "run_" + geneName + ".ld";
		System.out.println(filename); 
		File file = new File(filename);
	      try {
			file.createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
		FileExport writer = new FileExport(filename);
		
		for (int i = 0; i < cor.length; i++){
			String line = "";
			for (int j = 0; j < cor[i].length; j++){
				if (j == 0){
					line += cor[i][j];
				}
				else {
					line += "\t" + cor[i][j];
				}				
			}
			writer.println(line);
		}
		writer.close();
		Pascal.println("");
	}

	public static void writeLdMTJ(DenseMatrix cor, String geneName, String additonalDirectory) {
		
		String filename = Settings.outputDirectory_ + "/" + additonalDirectory + "/" + "run_" + geneName + ".ld";
		System.out.println(filename); 
		File file = new File(filename);
	      try {
			file.createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
		FileExport writer = new FileExport(filename);
		
		for (int i = 0; i < cor.numRows(); i++){
			String line = "";
			for (int j = 0; j < cor.numColumns(); j++){
				if (j == 0){
					line += cor.get(i,j);
				}
				else {
					line += "\t" + cor.get(i,j);
				}				
			}
			writer.println(line);
		}
		writer.close();
		Pascal.println("");
	}

}
