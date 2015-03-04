package ch.unil.genescore.vegas;

import java.util.LinkedList;
import java.util.List;

public interface AccessGenewiseSnpWeights {

	SnpWeightPairs getOverlappedSnpsWithWeights(String geneId, List<Snp> SnpList);
		

}
