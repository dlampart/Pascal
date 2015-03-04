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