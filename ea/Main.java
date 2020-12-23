package ea;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

public class Main {

	public static void main(String[] args) throws FileNotFoundException, ArrayDimensionsException {
		//constructor: population size, recombination probability, mutation probability
		Tsp tsp = new Tsp("cities90_11.csv", 100,1,0.7);
		//parameters: tsp instance, termination condition (number of generations)
		tsp.run(tsp,100);
	}
}
