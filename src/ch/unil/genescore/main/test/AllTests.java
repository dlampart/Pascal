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
package ch.unil.genescore.main.test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import ch.unil.genescore.gene.test.GenomicElementTest;
import ch.unil.genescore.pathway.test.RankSumTestTest;
import ch.unil.genescore.pathway.test.geneSetLibraryTest;
import ch.unil.genescore.vegas.test.AnalyticVegasTest;
import ch.unil.genescore.vegas.test.DistributionMethodsTest;
import ch.unil.genescore.vegas.test.FarebrotherTest;
import ch.unil.genescore.vegas.test.FileParserTest;
import ch.unil.genescore.vegas.test.GeneResultsNoScoreTest;
import ch.unil.genescore.vegas.test.GeneResultsSnpsOutOfBoundsTest;
import ch.unil.genescore.vegas.test.MTJConvenienceMethodsTest;
import ch.unil.genescore.vegas.test.MaxVegasTest;
import ch.unil.genescore.vegas.test.OverlappedCollectionStreamTest;
import ch.unil.genescore.vegas.test.OverlappedGenomicElementTest;
import ch.unil.genescore.vegas.test.SnpTest;
@RunWith(Suite.class)
//@SuiteClasses({ GenomicElementTest.class, SnpTest.class, AnalyticVegasTest.class, PathwayMainTest.class,MTJConvenienceMethodsTest.class})
@SuiteClasses({
	
	GenomicElementTest.class,
	
	
	RankSumTestTest.class,
	geneSetLibraryTest.class,
	FarebrotherTest.class,
	AnalyticVegasTest.class,
	DistributionMethodsTest.class,
	
	MaxVegasTest.class,
	MTJConvenienceMethodsTest.class,
	OverlappedCollectionStreamTest.class,
	OverlappedGenomicElementTest.class,
	
	SnpTest.class,
	GeneResultsSnpsOutOfBoundsTest.class,
	GeneResultsNoScoreTest.class,
	FileParserTest.class,

	
	
	
})
public class AllTests {
	
}
