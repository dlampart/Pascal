package ch.unil.genescore.vegas;
//TODO: orderCheck throw error from getNextAsOverlappedGenomicElement();	
public interface OverlappedGenomicElementStream {
	
	public OverlappedGenomicElement getNextAsOverlappedGenomicElement();
	public boolean streamOpen();
	public void reOpenStream();

}
