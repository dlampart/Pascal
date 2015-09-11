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
package ch.unil.genescore.gene;

import java.util.HashMap;
import java.util.HashSet;

import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;


/**
 * 
 */
public class GeneIdMapping {

	/** The unique instance of the mapping (Singleton design pattern) */
	static private GeneIdMapping instance_ = null;
	
	/** Mapping ensembl to entrez */
	private HashMap<String, HashSet<String>> ensembl2entrez_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Get the unique instance */
	public static GeneIdMapping getInstance() {
	
		if (instance_ == null)
			instance_ = new GeneIdMapping();
		
		return instance_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Map ensembl to entrez ids */
	public HashSet<String> ensembl2entrez(String ensemblId) {
		
		return ensembl2entrez_.get(ensemblId);
	}

		
	// ============================================================================
	// PRIVATE METHODS

	/** Private constructor, loads mapping */
	private GeneIdMapping() {
		
		load(Settings.geneIdMappingFile_);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Load the mapping */
	private void load(String filename) {
		
		ensembl2entrez_ = new HashMap<String, HashSet<String>>();
		FileParser parser = new FileParser(filename);
		
		while(true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			
			// Check number of columns
			if (nextLine.length != 3)
				Main.error("Expected three columns (ensembl id, entrez id, gene symbol)");
			
			// Parse ensembl id
			String ensg = nextLine[0];
			if (!(ensg.length() > 4 && ensg.substring(0, 4).equals("ENSG")))
				Main.error("Invalid ENSEMBL gene ID (expected 'ENSG...'): " + ensg);
			ensg = GeneAnnotationGencode.removeEnsemblVersion(ensg);

			// Parse entrez id
			String entrez = nextLine[1];
			if (entrez.length() > 0) {
				try {
					Integer.valueOf(entrez);
				} catch (NumberFormatException e) {
					Main.error("Invalid Entrez gene ID (expected an integer number): " + entrez);
				}
			}
			
			// Parse gene symbol
			//String symbol = nextLine[2];
			
			HashSet<String> entrezSet = ensembl2entrez_.get(ensg);
			if (entrezSet == null) {
				entrezSet = new HashSet<String>(1);
				ensembl2entrez_.put(ensg, entrezSet);
			}
			entrezSet.add(entrez);
		}
	}


	
	// ============================================================================
	// GETTERS AND SETTERS
		
	
}
