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

import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

//import cern.colt.matrix.tdouble.DoubleMatrix2D;
//import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;


/**
 * Some static utility functions
 */
public class Utils {

	/** Scientific format with plenty of digits (no loss in precision) */
	private static DecimalFormat scientific_ = new DecimalFormat("0.###############E0#####");
	/** Scientific format with around 10 digits (good for writing doubles to file) */
	private static DecimalFormat scientific10_ = new DecimalFormat("0.#######E0#");
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Take a time in [ms] and convert it into [h min s ms].
	 * @param dt [ms]
	 * @return Time [h m s ms]
	 */
	public static String chronometer(long dt) {
		
		int numHours = 0;
		int numMinutes = 0;
		int numSeconds = 0;
		
		//System.out.println(dt);
		
		numHours = (int)Math.floor(dt / 3600000.0);
		dt -= numHours * 3600000.0;
		
		numMinutes = (int)Math.floor(dt / 60000.0);
		dt -= numMinutes * 60000.0;
		
		numSeconds = (int)Math.floor(dt / 1000.0);
		dt -= numSeconds * 1000.0;
		
		String time = Integer.toString(numHours) + "h ";
		time += Integer.toString(numMinutes) + "min ";
		time += Integer.toString(numSeconds) + "s ";
		time += Integer.toString((int)dt) + "ms";
		
		return time;
	}
		
	// ----------------------------------------------------------------------------
	
	/** String representation of array, separate */
	public static String array2string(ArrayList<Double> v, String sep) {	
		
		if (v == null || v.size() == 0)
			return "";
		
		// The first element
		String str = Double.toString(v.get(1));
		if (v.size() == 1)
			return str;
		
		for (int i=1; i<v.size(); i++)
			str += sep + v.get(i);
		
		return str;
	}

	
	// ----------------------------------------------------------------------------

	/** Return true if the given array is an ordered list of positive integers, given in increasing order */
	public static boolean posIntIncreasing(ArrayList<Integer> x) {
		
		if (x == null || x.size() == 0)
			return true;
		
		for (int i=0; i<x.size(); i++) {
			if (x.get(i) < 1)
				return false;
			if (i > 0 && x.get(i) <= x.get(i-1))
				return false;
		}
		return true;
	}

	
	// ----------------------------------------------------------------------------

	/** Return true if the given array is an ordered list of positive doubles, given in increasing order */
	public static boolean posDoubleIncreasing(ArrayList<Double> x) {
		
		if (x == null || x.size() == 0)
			return true;
		
		for (int i=0; i<x.size(); i++) {
			if (x.get(i) <= 0)
				return false;
			if (i > 0 && x.get(i) <= x.get(i-1))
				return false;
		}
		return true;
	}

	
	// ----------------------------------------------------------------------------

	/** Extract the basic file name without path and without extension */
	static public String extractBasicFilename(String filename, boolean includePath) {
		 
		// The beginning of the filename (without the path) 
		int start = filename.lastIndexOf("/") + 1;
		if (start == -1)
			start = filename.lastIndexOf("\\") + 1; // windows
		if (start == -1)
			start = 0;
		
		// Remove .gz
		if (filename.endsWith(".gz"))
			filename = filename.substring(0, filename.length()-3);
		
		// The end of the filename (without file extension)
		int end = filename.lastIndexOf(".");
		if (end == -1 || end <= start) // not found or part of the path
			end = filename.length();

		// The output filename
		String basicFilename = filename.substring((includePath ? 0 : start), end);
		
		// Add custom suffix
		
		
	//	if (Settings.outputSuffix_ != null && Settings.outputSuffix_.compareTo("") != 0)
	//		basicFilename += Settings.outputSuffix_;
		
		return basicFilename;
	}

	
	// ----------------------------------------------------------------------------

	/** Scientific format */
	static public String toStringScientific(double x) {
		Double xd = new Double(x);
		if(xd.isNaN()){
			return "NA";
		}
		return scientific_.format(x);
	}

		
	// ----------------------------------------------------------------------------

	/** Scientific format with limited precision (around 10 digits, good for writing doubles to file) */
	static public String toStringScientific10(double x) {
		Double xd = new Double(x);
		if(xd.isNaN()){
			return "NA";
		}
		return scientific10_.format(x);
	}

	
	// ----------------------------------------------------------------------------

	/** Fixed length string (pad with white space if shorter) */
	public static String padRight(String s, int n) {
	
		return String.format("%1$-" + n + "s", s);  
	}
	
	
	// ----------------------------------------------------------------------------

	/** Convert dense Colt to Apache Commons matrix 
	static public RealMatrix colt2apache(DoubleMatrix2D colt) {
		
		RealMatrix apache = MatrixUtils.createRealMatrix(colt.rows(), colt.columns());
		for (int i=0; i<colt.rows(); i++)
			for (int j=0; j<colt.columns(); j++)
				apache.setEntry(i, j, colt.get(i, j));
		
		return apache;
	}

	
	// ----------------------------------------------------------------------------

	/** Convert dense Apache Commons to Colt matrix 
	static public DoubleMatrix2D apache2colt(RealMatrix apache) {
		
		// Convert apache commons matrix to colt
		DoubleMatrix2D colt = new DenseDoubleMatrix2D(apache.getRowDimension(), apache.getColumnDimension());
		for (int i=0; i<apache.getRowDimension(); i++)
			for (int j=0; j<apache.getColumnDimension(); j++)
				colt.set(i, j, apache.getEntry(i, j));
		
		return colt;
	}
*/
	
//	// ----------------------------------------------------------------------------
//
//	/** Get a URL from a relative path */
//	public URL getURL(String path)
//	{	
//		try
//		{
//			return new File(path).toURI().toURL();
//		} 
//		catch (MalformedURLException e)
//		{
//			log_.log(Level.WARNING, "NgseaSettings::getURL(): Bad URL", e);
//		}
//		return null;
//	}


}
