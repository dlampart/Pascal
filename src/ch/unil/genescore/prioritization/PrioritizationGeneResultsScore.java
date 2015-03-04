package ch.unil.genescore.prioritization;

import org.apache.commons.io.FilenameUtils;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.GeneResultsScore;
import ch.unil.genescore.vegas.GeneScoreEvaluator;

public class PrioritizationGeneResultsScore extends GeneResultsScore {
	
	@Override
	public void setExporter(String additionalOutputFileSuffix){
		String filename = getOutputString(additionalOutputFileSuffix);
		exporter_.setWriter(filename);			
	}
	
	
	
	private String getOutputString(String additionalOutputFileSuffix){
		
		
		if(additionalOutputFileSuffix==null){
			throw new RuntimeException("null not accepted as argument.");
		}
		
		String netpath = FilenameUtils.getBaseName(Settings.netPath_);
		String inputData = null;
		if (Settings.loadScoresFromFiles_)
			inputData = FilenameUtils.getBaseName(Settings.geneScoreFile_);
		else
			inputData = Settings.gwasName_;
		String outStr = Settings.outputDirectory_ + "/" + inputData + "." + Settings.outputSuffix_ + addDot(netpath) + addDot(additionalOutputFileSuffix) +".prioritization" + Settings.chromFileExtension_  + ".txt";
		return outStr;
	}	
}
