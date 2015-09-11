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
