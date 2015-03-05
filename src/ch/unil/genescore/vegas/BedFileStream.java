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

import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.FileParser;


public class BedFileStream 
extends OverlappedGenomicElementFileStream
implements OverlappedGenomicElementStream {
	static String fakeIds = ".";
	FileParser parser_ = null;
	String filePath_ = null;
	//OverlappedGenomicElement currentElement_ = null;
	
	
	public BedFileStream(String filePath){
		filePath_= filePath;
    	parser_ = new FileParser(filePath,"\t");
    	streamOpen_=true;
    	           		
	}	
	public void loadElement(){
		if (!streamOpen_)
			throw new RuntimeException("never use on closed stream");
		String[] nextLine = parser_.readLine();
		if (nextLine==null){
			currentElement_ = null;
			parser_.close();
			streamOpen_=false;	
			return;
		}
		GenomicElement currentElement = null;
		if (nextLine.length>3){
			currentElement =  new GenomicElement(nextLine[3]);
		}
		else {
			currentElement =  new GenomicElement(fakeIds);
		}
		currentElement.setPosition(nextLine[0], Integer.parseInt(nextLine[1]), Integer.parseInt(nextLine[2]), true);
		currentElement_= new OverlappedGenomicElement(currentElement);			
	}
	public boolean streamOpen(){ return streamOpen_;}
	@Override
	public void reOpenStream() {
		
		if(streamOpen()){throw new RuntimeException("never call reopen on open stream");}
    	parser_ = new FileParser(filePath_,"\t");
    	streamOpen_=true;
		// TODO Auto-generated method stub
		
	}	
}
