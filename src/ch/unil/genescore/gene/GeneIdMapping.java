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
