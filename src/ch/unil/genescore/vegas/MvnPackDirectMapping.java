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

import com.sun.jna.Native;
//import com.sun.jna.NativeLibrary;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

public class MvnPackDirectMapping  {//implements MvnPack {

	static {
		
        //Native.register("/Users/dlampart/Documents/workspace/genescore/fortranlibs/libmvtpack.dlyb");
		//System.loadLibrary(libname);
	//	NativeLibrary.addSearchPath("mvtpack", "fortranlibs/libmvtpack.dlyb");
	//	String str = NativeLibrary.getNativeLibraryResourcePrefix();
		
		Native.register("mvtpack");
    }

	public static native void mvtdst_(IntByReference n,IntByReference df, double[] lower, double[] upper, int[] infin, double[] correl, double[] delta,
			IntByReference maxpts, DoubleByReference abseps, DoubleByReference releps, 
			DoubleByReference error, DoubleByReference value, IntByReference inform);

	

}
