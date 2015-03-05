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

import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;


/**
 * Summarize p-values of multiple snps at a given locus
 */
public class EmpiricalPvalue {

	/** The observed test statistic */
	private double testStatisticReal_ = -1;
	/** The number of samples that are greater or equal than the test statistic */
	private int numSamplesGreater_ = -1;
	/** The total number of samples done (how many times addSample() was called since the last reset()) */
	private int numSamples_ = -1;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public EmpiricalPvalue(double testStatisticObserved) {

		testStatisticReal_ = testStatisticObserved;
	}

	
	// ----------------------------------------------------------------------------

	/** Add the next sample */
	public void addSample(double x) {
		
		// Doesn't make a difference if we use > or >= because we're comparing random doubles
		if (x > testStatisticReal_)
			numSamplesGreater_++;
	}

	
	// ----------------------------------------------------------------------------

	/** Reset the sample counts */
	public void start() {
		
		numSamplesGreater_ = 0;
		numSamples_ = 0;
	}

	
	// ----------------------------------------------------------------------------

	/** Give the number of samples that were done since the last reset(), returns true if stopping criteria are met */
	public boolean stop(int numSamples) {
		
		numSamples_ = numSamples;
		return numSamplesGreater_ > Settings.numSamplesGreaterCutoff_;
		//return false;
	}

	
	// ----------------------------------------------------------------------------

	/** Get the current p-value (assumes that addSample() has been called numSample times) */
	public double getPval() {
		
		// Check that there are more samples greater than the observed statistic than the cutoff, or the max number of samples was reached
		assert numSamplesGreater_ > Settings.numSamplesGreaterCutoff_ || numSamples_ == Settings.adaptiveNumSamples_.get(Settings.adaptiveNumSamples_.size()-1);

		// Note, the correct (unbiased) estimate of a p-value from monte carlo sampling is not r/n, but (r+1)/(n+1)
		if (numSamples_ <= 0)
			Main.error("Number of samples <= 0, stop() may not have been called");
		
		return (numSamplesGreater_ + 1.0) / (numSamples_ + 1.0);
	}

	
	// ============================================================================
	// PRIVATE METHODS
	

	// ============================================================================
	// GETTERS AND SETTERS

	public double getTestStatisticReal() { return testStatisticReal_; }
	public int getNumSamplesGreater() { return numSamplesGreater_; }
	public int getNumSamples() { return numSamples_; }
	
}
