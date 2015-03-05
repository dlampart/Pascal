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

import java.io.DataInputStream;
import java.io.IOException;

import ch.unil.genescore.main.Main;

public class SnpPositionStream
	extends OverlappedGenomicElementFileStream
	implements OverlappedGenomicElementStream {
	
	DataInputStream inStream_ = null;
	String binaryFileVersionID_= null;
	String filename_ = null;
	
	
	
	public SnpPositionStream(String filename,String binaryFileVersionID){
		System.out.println("reading snp positions from file:" + filename);
		DataInputStream inStream = null;
		try {						
			inStream = StreamMethods.openDataInputStream(filename, binaryFileVersionID);
			}
			catch(Exception e){
				Main.error(e, "Error loading SNP position (try deleting the binary reference population files)");
			}	
		inStream_ = inStream;
		streamOpen_= true;
		filename_ = filename;
		binaryFileVersionID_ = binaryFileVersionID;
	} 
	protected void loadElement()  {
		
		Snp currentSnp = null;
		if (!streamOpen_)
			throw new RuntimeException("never use on closed stream");
	
		String snpId = null;
		
		try {
			snpId = inStream_.readUTF();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// End of file is marked with the verion ID
		if (snpId.equals(binaryFileVersionID_)){
			streamOpen_=false;
			currentElement_ = null;		
			try {
				inStream_.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return;
		}		
		currentSnp = new Snp(snpId);		
		currentElement_ = new OverlappedGenomicElement(currentSnp);
		
		try {
			currentSnp.readPosAndAllele(inStream_);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		catch (DataInconsistencyException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}
	
	public Snp loadSnp(){
		
		loadElement();		
		return ((Snp) currentElement_.getMainElement());		
	}
	
	public boolean streamOpen(){ return streamOpen_;}
	public void reOpenStream(){
		if(streamOpen()){throw new RuntimeException("never call reopen on open stream");}
		DataInputStream inStream = null;
		try {						
			inStream = StreamMethods.openDataInputStream(filename_, binaryFileVersionID_);
			}
			catch(Exception e){
				Main.error(e, "Error loading SNP position (try deleting the binary reference population files)");
			}	
		inStream_ = inStream;
		streamOpen_= true;
		streamInitialized_ = false;

		
	}
}
