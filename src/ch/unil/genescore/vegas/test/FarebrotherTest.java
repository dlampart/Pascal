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
package ch.unil.genescore.vegas.test;
import static org.junit.Assert.assertEquals;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.Farebrother;
import ch.unil.genescore.vegas.FarebrotherExperiment;


public class FarebrotherTest {

	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
	}

	@AfterClass
	public static void testCleanup() { }
	
	
	// ============================================================================
	// TESTS

	/** Test test statistic */
	@Test
	public void FarebrotherTest() {
		Settings.requestedAbsolutePrecision_=1e-25;
		double w=1e-17;
		System.out.println(w);
		double q=1-w;
		double t=1-q;
		 t=w/10.0;
		System.out.println(t);
		double[] lambda={4,3,2,1};		
		
		Farebrother myfarebrother = new Farebrother(lambda);
		double res = myfarebrother.probQsupx(300);
		System.out.println(res);
		FarebrotherExperiment myfarebrotherExp = new FarebrotherExperiment(lambda);
		double res2 = myfarebrotherExp.probQsupx(300);			
		System.out.println(res2);
	}
	
}
