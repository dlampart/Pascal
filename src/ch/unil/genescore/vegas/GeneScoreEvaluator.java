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

/** 
 * Interface for classes that can evaluate gene scores
 */
public abstract class GeneScoreEvaluator {
	
	
	/** Compute the gene score / p-value, return true if success */
	public abstract boolean computeScore();
	/** Get the score(s) */
	public abstract double[] getScore();

	/** Get a string representation of the results */
	public abstract String getResultsAsString();
	/** Get a header line corresponding to getResultsAsString() */
	public abstract String getResultsAsStringHeader();
	public abstract void setDataFromGeneData(GeneData Dat);		
	
	
	/** Get output to be printed on console after computing score for a gene */
	public String getConsoleOutput() { return ""; }
	
	/** Get output to be printed for genes where no score could be computed (error) */
	public String getNoScoreOutput() { return ""; }
	/** Get a header line corresponding to getNoScoreOutput() */
	public String getNoScoreOutputHeader() { return ""; }

}
