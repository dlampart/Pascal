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

import java.io.File;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileExport;

public class GeneResultsScore {
	
	private boolean headerWritten_ = false;
	protected FileExport exporter_ = null;
	
	public GeneResultsScore(){	
		//exporter_ = new FileExport();
	}
	

	protected String addDot(String withoutDot){
		String outStr=withoutDot;
		if (withoutDot!="")
			outStr = "." + withoutDot;
		return outStr;
	}
	
	public void setExporter(String additionalOutputFileSuffix){
		File file = new File(Pascal.set.outputDirectory_, 
				Pascal.set.gwasName_ + addDot(additionalOutputFileSuffix) + ".genescores" + Pascal.set.chromFileExtension_  + ".txt");
		//exporter_.setWriter(filename);
		exporter_ = new FileExport(Pascal.log, file);
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
