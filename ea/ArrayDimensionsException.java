package ea;
public class ArrayDimensionsException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public ArrayDimensionsException() {
		super("Matrix must be nxn i.e. equal number of columns and rows.");
	}
	
	public ArrayDimensionsException(String message) {
		super(message);
	}
	
	
}