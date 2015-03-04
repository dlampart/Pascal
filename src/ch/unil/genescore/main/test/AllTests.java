package ch.unil.genescore.main.test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import ch.unil.genescore.eqtlOverlap.test.EqtlEvaluatorTest;
import ch.unil.genescore.eqtlOverlap.test.EqtlGenePairTest;
import ch.unil.genescore.eqtlOverlap.test.EqtlResultsTest;
import ch.unil.genescore.eqtlOverlap.test.GeneDataEqtlProjectionTest;
import ch.unil.genescore.gene.test.GenomicElementTest;
import ch.unil.genescore.pathway.test.DamperTest;
import ch.unil.genescore.pathway.test.PathwayMainTest;
import ch.unil.genescore.pathway.test.RankSumTestTest;
import ch.unil.genescore.pathway.test.geneSetLibraryTest;
import ch.unil.genescore.prioritization.test.PrunedConnectionTreGetterTest;
import ch.unil.genescore.prioritization.test.SparseNetTest;
import ch.unil.genescore.projection.SnpWeightCreatorTest;
import ch.unil.genescore.projection.test.ProjectionAcceptanceTest;
import ch.unil.genescore.projection.test.ProjectionAcceptanceWithoutImhofTest;
import ch.unil.genescore.vegas.GeneDataFakePhenotypeIntegrationTest;
import ch.unil.genescore.vegas.test.AllOverlappedTest;
import ch.unil.genescore.vegas.test.AnalyticVegasTest;
import ch.unil.genescore.vegas.test.BedFileStreamTest;
import ch.unil.genescore.vegas.test.DistributionMethodsTest;
import ch.unil.genescore.vegas.test.FarebrotherTest;
import ch.unil.genescore.vegas.test.FileParserTest;
import ch.unil.genescore.vegas.test.GeneDataProjectionTest;
import ch.unil.genescore.vegas.test.GeneResultsNoScoreTest;
import ch.unil.genescore.vegas.test.GeneResultsSnpsOutOfBoundsTest;
import ch.unil.genescore.vegas.test.GenewiseSnpWeightsTest;
import ch.unil.genescore.vegas.test.GenomeWideScoringTest;
import ch.unil.genescore.vegas.test.LdMatTest;
import ch.unil.genescore.vegas.test.LinkageDisequilibriumTest;
import ch.unil.genescore.vegas.test.MTJConvenienceMethodsTest;
import ch.unil.genescore.vegas.test.MaxVegasTest;
import ch.unil.genescore.vegas.test.OverlappedCollectionStreamTest;
import ch.unil.genescore.vegas.test.OverlappedGenomicElementTest;
import ch.unil.genescore.vegas.test.ProjectionVegasPipelineTest;
import ch.unil.genescore.vegas.test.ProjectionVegasTest;
import ch.unil.genescore.vegas.test.ReferencePopulationTest;
import ch.unil.genescore.vegas.test.SnpTest;
import ch.unil.genescore.vegas.test.SnpWeightMapTest;
import ch.unil.genescore.vegas.test.VegasSimulationNewTest;
import ch.unil.genescore.vegas.test.tmpTest;

@RunWith(Suite.class)
//@SuiteClasses({ GenomicElementTest.class, SnpTest.class, AnalyticVegasTest.class, PathwayMainTest.class,MTJConvenienceMethodsTest.class})
@SuiteClasses({
	//TODO: classes not properly covered:
	//Main, PathwayMain,SnpWeightCreator, 
	
	//prioritizationTests
	PrunedConnectionTreGetterTest.class,
	SparseNetTest.class,
	
	
	//eqtlOverlap/test/
	EqtlEvaluatorTest.class,
	EqtlGenePairTest.class,
	EqtlResultsTest.class,
	GeneDataEqtlProjectionTest.class,

	//gene/test
	GenomicElementTest.class,
	
	//pathway/test/
	DamperTest.class,
	geneSetLibraryTest.class,
	//PathwayMainTest.class, isn't proper unit-test
	RankSumTestTest.class,
	
	//projection
	//SnpWeightCreatorTest.class, isn't proper unit-test
	// ProjectionAcceptanceTest.class, isn't proper unit-test
	//ProjectionAcceptanceWithoutImhofTest.class, isn't proper unit-test
	
	//vegas
	//GeneDataFakePhenotypeIntegrationTest.class, isn't proper unit-test
	AllOverlappedTest.class,
	AnalyticVegasTest.class,
	BedFileStreamTest.class,
	DistributionMethodsTest.class,
	FarebrotherTest.class,
	GeneDataProjectionTest.class,
	GenewiseSnpWeightsTest.class,
	GenomeWideScoringTest.class,
	LdMatTest.class,
	LinkageDisequilibriumTest.class,
	MaxVegasTest.class,
	MTJConvenienceMethodsTest.class,
	OverlappedCollectionStreamTest.class,
	OverlappedGenomicElementTest.class,
	ProjectionVegasTest.class,
	//ReferencePopulationTest.class, isn't proper unit test
	SnpTest.class,
	SnpWeightMapTest.class,
	tmpTest.class,
	VegasSimulationNewTest.class,
	GeneResultsSnpsOutOfBoundsTest.class,
	GeneResultsNoScoreTest.class,
	FileParserTest.class,

	
	
	
})
public class AllTests {
	
}
