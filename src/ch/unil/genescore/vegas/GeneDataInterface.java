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

import no.uib.cipr.matrix.DenseMatrix;

public interface GeneDataInterface {

	public DenseMatrix getCorr();
	public ArrayList<Double> getScores();
	public double[] getWeights();
	public void processData();	
}
