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

public class SnpSerializedPositionStream {

	private int lineCounter_ = 0;
	/** Next line */
	private Snp nextLine_ = null;
	/** current line */
	private Snp currentLine_ = null;
	
	DataInputStream inStream_ = null;
	String binaryFileVersionID_= null;
	String filename_ = null;
	boolean streamOpen_ = true;
	
	
	public SnpSerializedPositionStream(String filename,String binaryFileVersionID){
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
		readLine();
	}
		public boolean lineAvailable(){
			return(currentLine_!=null);
		}
		
		private Snp createSnpFromStream(){
			
			if (!streamOpen_)
				throw new RuntimeException("Don't call on closed stream");
		
			String snpId = null;
			
			try {
				snpId = inStream_.readUTF();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(-1);
			}
			// End of file is marked with the verion ID
			if (snpId.equals(binaryFileVersionID_)){
				streamOpen_=false;
				
				try {
					inStream_.close();
				} catch (IOException e) {
					Main.error(e);
				}
				return null;
			}		
			Snp currentSnp = new Snp(snpId);		
			
			try {
				currentSnp.readPosAndMinorAllele(inStream_);
			} catch (IOException e) {
				Main.error(e, "Error loading SNP position (try deleting the binary reference population files)");				
			} 
			catch (DataInconsistencyException e) {
				Main.error(e);
			 }
			return(currentSnp);
		}
					
		public Snp readLine() {
			Snp outLine = currentLine_;
			if (lineCounter_==0){
				currentLine_ = createSnpFromStream();				
				outLine=currentLine_;
			}
			else{
				currentLine_=nextLine_;				
			}
			lineCounter_++;
			if (streamOpen_)
				nextLine_ = createSnpFromStream();				
			return outLine;	
		} 
}
