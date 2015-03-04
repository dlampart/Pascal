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
}
