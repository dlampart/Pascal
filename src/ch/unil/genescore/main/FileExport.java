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

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;


/**
 * Write a text file
 */
public class FileExport {

	/** The buffered file writer */
	BufferedWriter writer_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	    
	
	public FileExport(){};
	public FileExport(String filename, boolean gzip) {

		setWriter(filename, gzip);
	}
	
	/** Constructor for uncompressed file */
	public FileExport(String filename) {
		
		this(filename, false);
	}

	public void setWriter(String filename){
		setWriter(filename, false);
	}
	
	public void setWriter(String filename, boolean gzip){
		
		try {
			if (gzip)
				filename += ".gz";
			System.out.println("Writing file: " + filename);
		
			if (gzip) {
				FileOutputStream output = new FileOutputStream(filename);
				writer_ = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8"));
			} else {				
				FileWriter fstream = new FileWriter(filename);
				writer_ = new BufferedWriter(fstream);
			}
			
		} catch (Exception e) {
			Main.error(e);
		}
	}
	
	
    // ----------------------------------------------------------------------------

	/** Write a line to the file */
	public void println(String str) {
		
		print(str + "\n");
	}


	// ----------------------------------------------------------------------------

	/** Write the given string to the file */
	public void print(String str) {
		
		try {
			writer_.write(str);
		} catch (IOException e) {
			Main.error(e);
		}
	}
	
	
    // ----------------------------------------------------------------------------

	/** Be polite and close the file writer */
	public void flush() {
		
		try {
			writer_.flush();
		} catch (IOException e) {
			Main.error(e);
		}
	}

	
    // ----------------------------------------------------------------------------

	/** Be polite and close the file writer */
	public void close() {
		
		try {
			writer_.close();
		} catch (IOException e) {
			Main.error(e);
		}
	}
	  
}
