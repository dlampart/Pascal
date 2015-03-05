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
package ch.unil.genescore.main.test;

import java.util.ArrayList;

import ch.unil.genescore.main.FileParser;

public class ParserStub extends FileParser {

	static ArrayList<String> TestData_ = null;
	int counter = 0;

	public ParserStub() {
	}

	public ParserStub(String filename) {
		super();
	}

	public static void setTestData(ArrayList<String> TestData) {
		TestData_ = TestData;
	}

	@Override
	public String[] readLine() {

		if (counter < TestData_.size()) {
			String[] out = TestData_.get(counter).split("\t", -1);
			counter++;
			return out;
		} else
			return null;
	}

	@Override
	public void close() {
	}
}
