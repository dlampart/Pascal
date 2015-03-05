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

package ch.unil.genescore.main;


/**
 * Small class to output progress when running a loop
 */
public class ProgressMonitor {
	
	/** Total iterations */
	private int totalIterations_ = 0;
	/** Output frequency */
	private int freq_ = -1;
	/** Starting time */
	private long t0_ = 0;
	/** The string written at each iteration */
	private String asterisks_ = "*";

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public ProgressMonitor(int totalIterations) {

		totalIterations_ = totalIterations;
		
		freq_ = totalIterations / 40;
		if (freq_ == 0) {
			freq_ = 1;
			asterisks_ = "";
			for (int i=0; i<40/totalIterations_; i++)
				asterisks_ += "*";
		}
		
		t0_ = System.currentTimeMillis();
		
		Main.println("|-------------- Progress --------------|");
	}
	
	
	/** Constructor */
	public ProgressMonitor(int totalIterations, int freq) {
		
		this(totalIterations);
		freq_ = freq;
	}

	
	// ----------------------------------------------------------------------------

	/** Print progress */
	public void iteration(int i) {
		
		if (i % freq_ == 0)
			System.out.print(asterisks_);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Estimated total runtime */
	public void estimatedTotalRuntime(int i) {
		if (i % freq_ == 0) {
			long t1 = System.currentTimeMillis();
			Main.println(i + "\tERT: \t" + Utils.chronometer(totalIterations_*(t1-t0_)/i));
		}
	
	}
	
	
	// ----------------------------------------------------------------------------

	/** Print progress */
	public void done() {
		Main.print("\n\n");
	}

}
