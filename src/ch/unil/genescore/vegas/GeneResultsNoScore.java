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

import java.util.ArrayList;

import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;

public class GeneResultsNoScore {
	
	private ArrayList<String> noScore_ = null;
	private FileExport exporter_ = null;
	public GeneResultsNoScore(){
		noScore_ = new ArrayList<String>(); 
		exporter_ = new FileExport();
	}

	public void add(String str){
		noScore_.add(str);
	}
	
	public void setExporter(FileExport exp) {
		exporter_=exp;
	}
	
	public void setExporter(String additionalOutputFileSuffix) {
		String filename = Settings.outputDirectory_ + "/" + Settings.gwasName_+ additionalOutputFileSuffix + ".scoreComputeError" + Settings.chromFileExtension_ + ".txt";		
		exporter_.setWriter(filename);		
	}
	
	
public void writeResultsToFile(String additionalOutputFileSuffix){
		
		if (noScore_.size() == 0)
			return;
		
		Main.warning("Gene score computation did not converge at specified precision for some genes");
		setExporter(additionalOutputFileSuffix);
		String header = "chromosome\tstart\tend\tstrand\tgene_id\tsymbol\tScore\tStatus";		
		exporter_.println(header);		
		for (String line : noScore_)
			exporter_.println(line);
		exporter_.close();
		Main.println("");
	}
	
}
