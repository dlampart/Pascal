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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;

public class StreamMethods {
	
	
	public static DataOutputStream openDataOutputStream(String filename , String binaryFileVersionID)  {

		Main.println("Writing file: " + filename);

		try {
			FileOutputStream outfile = new FileOutputStream(filename);
			GZIPOutputStream gzip = new GZIPOutputStream(outfile);
			BufferedOutputStream buf = new BufferedOutputStream(gzip);
			DataOutputStream outStream = new DataOutputStream(buf);
			
			// Write the version, used as a check when reading files
			outStream.writeUTF(binaryFileVersionID);
			
			return outStream;
			
		} catch (Exception e) {
			throw new RuntimeException("Could not open binary output file: " + filename);
		}
	}
	
	/** Open a data input stream */
	public static DataInputStream openDataInputStream(String filename,  String binaryFileVersionID)  {
		
		if (Settings.verbose_)
			Main.println("Reading file: " + filename);

		try {
			FileInputStream infile = new FileInputStream(filename);
			GZIPInputStream gzip = new GZIPInputStream(infile);
			BufferedInputStream buf = new BufferedInputStream(gzip);
			DataInputStream inStream = new DataInputStream(buf);
			
			// Check that the file starts with the right versionId
			String versionId = inStream.readUTF();
			if (!versionId.equals(binaryFileVersionID))
				throw new RuntimeException("Incompatible version ID of binary file, delete the file to create a new one");
			
			return inStream;
			
		} catch (Exception e) {
			throw new RuntimeException("Could not open binary input file: " + filename);
		}
	}

}
