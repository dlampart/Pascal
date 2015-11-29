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

import java.util.ArrayList;

public class ConvenienceMethods {
	/**converting between ArrayList and double[] or int[]*/
	public static ArrayList<Double> doubleArToArList(double[] myList){
		ArrayList<Double> arOut = new ArrayList<Double>(myList.length);
		for (double el : myList){
			arOut.add(el);
		}
		return arOut;
	}
	public static double[] arListToDoubleAr(ArrayList<Double> myList){
		double[] arOut = new double[myList.size()];
		for (int i=0 ; i < myList.size(); i++){
			arOut[i]=myList.get(i);
		}
		return arOut;
	}
	public static ArrayList<Integer> intArToArList(int[] myList){
		ArrayList<Integer> arOut = new ArrayList<Integer>(myList.length);
		for (int el : myList){
			arOut.add(el);
		}
		return arOut;
	}
	public static int[] arListToIntAr(ArrayList<Integer> myList){
		int[] arOut = new int[myList.size()];
		for (int i=0 ; i < myList.size(); i++){
			arOut[i]=myList.get(i);
		}
		return arOut;
	}
	/**subindexing double[] with other int[]*/
	public static double[] subIndexing(double[] ar,int[] indices){
		double[] arOut = new double[indices.length];
		for (int i=0; i<indices.length; i++){
			arOut[i]=ar[indices[i]];			
		}
		return arOut;
	}
	public static int[] subIndexing(int[] ar,int[] indices){
		int[] arOut = new int[indices.length];
		for (int i=0; i<indices.length; i++){
			arOut[i]=ar[indices[i]];			
		}
		return arOut;
	}
	
	// ----------------------------------------------------------------------------

	/** Fixed length string (pad with white space if shorter) */
	public static String padRight(String s, int n) {

		return String.format("%1$-" + n + "s", s);  
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
}
