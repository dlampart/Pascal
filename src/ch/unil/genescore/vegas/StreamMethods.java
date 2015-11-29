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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ch.unil.genescore.main.Pascal;

public class StreamMethods {
	
	
	public static DataOutputStream openDataOutputStream(String filename , String binaryFileVersionID)  {

		Pascal.println("Writing file: " + filename);

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
		
		if (Pascal.set.verbose_)
			Pascal.println("Reading file: " + filename);

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
