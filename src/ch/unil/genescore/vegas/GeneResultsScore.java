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

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Settings;

public class GeneResultsScore {
	
	private boolean headerWritten_ = false;
	protected FileExport exporter_ = null;
			//new FileExport();
	
	public GeneResultsScore(){	
		exporter_ = new FileExport();
	}
	

	protected String addDot(String withoutDot){
		String outStr=withoutDot;
		if (withoutDot!="")
			outStr = "." + withoutDot;
		return outStr;
	}
	
	public void setExporter(String additionalOutputFileSuffix){
		
		String filename = Settings.outputDirectory_ + "/" + Settings.gwasName_ + addDot(additionalOutputFileSuffix) + ".genescores" + Settings.chromFileExtension_  + ".txt";
		exporter_.setWriter(filename);			
	}
	public void writeLine(GeneScoreEvaluator evaluator, Gene gene){

		if (!getHeaderWritten()) {
			getExporter().println("chromosome\tstart\tend\tstrand\tgene_id\tgene_symbol" + evaluator.getResultsAsStringHeader());
			setHeaderWritten(true);
		}

	// Write results for this gene to file
		String nextLine = gene.toString();
		nextLine += evaluator.getResultsAsString();
		getExporter().println(nextLine);
		getExporter().flush();
	}
	
	/**write GeneScoreEvaluator output directly without gene information.*/
	public void writeLine(GeneScoreEvaluator evaluator){

		if (!getHeaderWritten()) {
			getExporter().println(evaluator.getResultsAsStringHeader());
			setHeaderWritten(true);
		}

		String nextLine = evaluator.getResultsAsString();
		getExporter().println(nextLine);
		getExporter().flush();
	}
	
	
	
	
	public FileExport getExporter() {
		return exporter_;
	}
	
	public boolean getHeaderWritten(){
		return headerWritten_;
	}
	public void setHeaderWritten(boolean h){
		headerWritten_=h;
	}

	public void setExporter(FileExport exporter) {
		this.exporter_ = exporter;
	}
	

}
