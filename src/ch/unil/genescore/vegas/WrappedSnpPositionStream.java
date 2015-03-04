package ch.unil.genescore.vegas;

import java.io.DataInputStream;
import java.util.ArrayList;

import ch.unil.genescore.main.Main;

public class WrappedSnpPositionStream 
//extends OverlappedGenomicElementFileStream
implements OverlappedGenomicElementStream {
	
	private SnpPositionStream posStream_ = null;
	private int currentIndex_ = 0;	
	private ArrayList<String> filenames_ = null;
	private String binaryFileVersionID_ = null;
	
	WrappedSnpPositionStream(ArrayList<String> filenames,String binaryFileVersionID){
		filenames_ = filenames;
		binaryFileVersionID_=binaryFileVersionID;
		currentIndex_=0;
		posStream_ = new SnpPositionStream(filenames_.get(currentIndex_),binaryFileVersionID_);
	}
	@Override
	public OverlappedGenomicElement getNextAsOverlappedGenomicElement() {
		boolean currentStreamOpen = posStream_.streamOpen();
		if (!currentStreamOpen) {
			if ((currentIndex_+1)==filenames_.size())
				throw new RuntimeException("never call on closed stream.");			
			currentIndex_++;
			posStream_ = new SnpPositionStream(filenames_.get(currentIndex_),binaryFileVersionID_);			
		}
		return posStream_.getNextAsOverlappedGenomicElement();			
	}

	@Override
	public boolean streamOpen() {
		// TODO Auto-generated method stub
		return (posStream_.streamOpen() || (currentIndex_+1)!=filenames_.size());
	}
	@Override
	public void reOpenStream() {
		if(streamOpen()){throw new RuntimeException("never call reopen on open stream");}
		currentIndex_=0;
		posStream_ = new SnpPositionStream(filenames_.get(currentIndex_),binaryFileVersionID_);
	}

	
}


