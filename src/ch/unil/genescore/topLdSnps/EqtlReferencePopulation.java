package ch.unil.genescore.topLdSnps;

import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.ReferencePopulation;

public class EqtlReferencePopulation extends ReferencePopulation{
	
	
	public void loadEqtlFile(){
		FileParser parser = new FileParser(Settings.eqtlFile_);
	}	
}
