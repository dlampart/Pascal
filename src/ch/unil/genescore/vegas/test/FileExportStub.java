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

import java.util.ArrayList;

import ch.unil.genescore.main.FileExport;

class FileExportStub extends FileExport {
	ArrayList<String> myStrs = new ArrayList<String>();
	boolean isClosed=false;
	@Override
	public void setWriter(String filename){};
	@Override
	public void println(String str){
		myStrs.add(str);			
	};
	@Override
	public void close(){isClosed=true;};
	public boolean isclosed(){return isClosed;}
	public ArrayList<String> getStrings(){return(myStrs);}
}