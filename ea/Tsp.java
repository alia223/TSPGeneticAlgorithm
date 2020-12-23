package ea;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;

public class Tsp {
	//2D array storing all routes through the map of cities
	private double[][] graph;
	//String representing name of file that contains all cities
	private String fileName;
	//Random used to generate random permutations of cities (routes through the map)
	private Random random;
	//population
	private List<List<Integer>> population = new ArrayList<>();
	//population size
	private int populationSize;
	//generated offspring
	private List<List<Integer>> offspringPopulation = new ArrayList<>();
	//cost of all routes in population
	private double[] fitness;
	//global best cost/fitness
	private double lowestCost = Integer.MAX_VALUE;
	//global best route
	private List<Integer> shortestRoute = new ArrayList<>(); 
	//probability of applying crossover to produce variation
	private double recombinationProbability;
	//probability of applying swap mutation to produce variation
	private double mutationProbability;

	/**
	 * Constructs Travelling Salesperson problem
	 * @param fileName String representing name of file that contains all cities
	 * @param population Size size of population
	 * @param recombinationProbability Probability of applying crossover to produce variation
	 * @param mutationProbability Probability of applying swap mutation to produce variation
	 * @throws ArrayDimensionsException
	 * @throws FileNotFoundException
	 */
	public Tsp(String fileName, int populationSize, double recombinationProbability, double mutationProbability) throws ArrayDimensionsException, FileNotFoundException{
		this.fileName = fileName;
		//Initialise map of cities (2D array)
		this.graph = createGraph();
		//Check that the map of cities being given is n x n
		for(int i = 0;i < graph.length;i++) {
			if(graph.length != graph[0].length) {
				throw new ArrayDimensionsException("n x n Matrix required i.e. equal number of columns and rows.");
			}
		}
		//Initialise random variable for later use when generating random routes.
		random = new Random();
		this.populationSize = populationSize;
		this.fitness = new double[populationSize];
		this.recombinationProbability = recombinationProbability;
		this.mutationProbability = mutationProbability;
	}
	
	/**
	 * Gets file name that contains distances between all cities
	 * @return fileName String representing name of file that contains distances between all cities
	 */
	private String getFileName() {
		return this.fileName;
	}
	
	/**
	 * Gets population
	 * @return population
	 */
	public List<List<Integer>> getPopulation() {
		return this.population;
	}
	
	/**
	 * Gets population size
	 * @return populationSize
	 */
	public int getPopulationSize() {
		return this.populationSize;
	}
	
	/**
	 * Gets offspring
	 * @return offspringPopulation Population of offspring generated via crossover and muation of parents
	 */
	public List<List<Integer>> getOffspringPopulation() {
		return this.offspringPopulation;
	}
	
	/**
	 * Gets shortest route
	 * @return shortestRotue Global best route
	 */
	public List<Integer> getShortestRoute() {
		return this.shortestRoute;
	}
	
	/**
	 * Gets cost of shortest route
	 * @return recordCost Global best cost
	 */
	public double getLowestCost() {
		return this.lowestCost;
	}
	
	/**
	 * Calculate cost of each route in population
	 */
	public void calculateFitness() {
		for(int i = 0;i < this.populationSize;i++) {
			double cost = getCostOfRoute(this.population.get(i));
			if(cost < lowestCost) {
				this.lowestCost = cost;
				this.shortestRoute.clear();
				this.shortestRoute.addAll(this.population.get(i));
			}
			this.fitness[i] = cost;
		}
	}
	
	/**
	 * Tournament Selection of parents to crossover and mutate to produce offspring
	 * @param population
	 * @return parentOne Parent selected 
	 * OR
	 * @return parentTwo
	 */
	public List<Integer> parentSelection(List<List<Integer>> population) {
		int indexOne = random.nextInt(population.size());
		int indexTwo = random.nextInt(population.size());
		List<Integer> parentOne = population.get(indexOne);
		List<Integer> parentTwo = population.get(indexTwo);
		
		while(indexOne == indexTwo) {
			indexTwo = random.nextInt(populationSize);
			parentTwo = population.get(indexTwo);
		}
		
		if(getCostOfRoute(parentOne) > getCostOfRoute(parentTwo)) {
			return parentOne;
		}
		else {
			return parentTwo;
		}
	}
	
	/**
	 * Selects most fit individuals between both parents and children
	 * @param population
	 * @return result Most fit individuals out of all the parents and offspring
	 */
	public List<List<Integer>> survivorSelection(List<List<Integer>> population) {
		int count = 0;
		//create a HasMap containing all routes in population and cost of each route
		HashMap<List<Integer>, Double> routeAndFitness = new HashMap<>();
		List<List<Integer>> result = new ArrayList<>();
		for(List<Integer> p: population) {
			routeAndFitness.put(p, getCostOfRoute(p));
		}
		//Sort routes in ascending order
		HashMap<List<Integer>, Double> sortedMap = sortByValue(routeAndFitness);
		//Select most fit individuals by selecting top [populationSize] indiviudals/routes
		for(Map.Entry<List<Integer>, Double> entry : sortedMap.entrySet()) {
			if (count < this.populationSize) {
				result.add(entry.getKey());
			}
			count++;
		}
		return result;
	}
	
	/**
	 * Produces next generation i.e. offspring of parents
	 * @param population
	 */
	public void doGeneration(List<List<Integer>> population) {
		//Ensure new population in each generation remains the same size forever
		while(offspringPopulation.size() < this.populationSize) {
			//Generate offspring
			List<Integer> parentOne = parentSelection(population);
			List<Integer> parentTwo = parentSelection(population);
			List<Integer> childOne = crossover(parentOne, parentTwo, this.recombinationProbability);
			List<Integer> childTwo = crossover(parentTwo, parentOne, this.recombinationProbability);
			mutate(childOne, this.mutationProbability);
			mutate(childTwo, this.mutationProbability);
			offspringPopulation.add(childOne);
			offspringPopulation.add(childTwo);
		}
		//add offspring to population
		this.population.addAll(offspringPopulation);
		//set population to contain only most fit individuals between parents and offspring
		this.population = survivorSelection(this.population);
		//empty offspring list, read for next generation
		this.offspringPopulation.clear();
	} 
	
	/**
	 * Produces one child by applying crossover between two parents
	 * @param parentOne Randomly selected parent
	 * @param parentTwo Randomly selected parent
	 * @param recombinationProbability Probability of applying crossover to produce variation
	 * @return child One child produced by crossover of parents
	 */
	public List<Integer> crossover(List<Integer> parentOne, List<Integer> parentTwo, double recombinationProbability) {
		double probability = random.nextDouble();
		List<Integer> child  = null;
		if(probability <= recombinationProbability) {
			int crossoverPosition = random.nextInt(parentOne.size() - 1) + 1;
			List<Integer> parentCopy = new ArrayList<>(parentOne);
			child = parentCopy.subList(0, crossoverPosition);
			for(int i = 0;i < parentTwo.size();i++) {
				int city = parentTwo.get(i);
				if(!child.contains(city)) {
					child.add(city);
				}
			}
		}
		else {
			child = parentOne;
		}
		return child;
	}
	
	//Swap two elements at specific indexes in a given list
	public void mutate(List<Integer> list, double mutationProbability) {
		double probability = random.nextDouble();
		if(probability < mutationProbability) {
				int i = random.nextInt(list.size());
				int j = random.nextInt(list.size());
				int[] arr = new int[list.size()];
				for(int index = 0;index < list.size();index++) {
					arr[index] = list.get(index);
				}
				int temp = arr[i];
				arr[i] = arr[j];
				arr[j] = temp;
				list.clear();
				for(int index = 0;index < arr.length;index++) {
					list.add(index, arr[index]);
			}
		}
	}
	
	//Generate a random legal route through the graph/matrix
	public List<Integer> generateRandomRoute() {
		List<Integer> route = new ArrayList<>();
		int routeSize = route.size();
		while(routeSize < this.graph.length) {
			int city = this.random.nextInt(this.graph.length) + 1;
			if(!route.contains(city)) {
				route.add(city);
				routeSize++;
			}
		}
		return route;
	}
	
	/**
	 * Get cost of a single route
	 * @param routeList List containing a single route
	 * @return cost Cost of single route given as a parameter
	 */
	public double getCostOfRoute(List<Integer> routeList) {
		//starting point/city
		int start = routeList.get(0)-1;
		//last city before you have to return to the starting point
		int finish = routeList.get(routeList.size()-1)-1;
		//cost of route
		double cost = 0;
		//start FROM one city TO another
		//because each city ranges from index 0..n it is easy to just iterate through 2D array of all the cities
		for(int i = 0;i < routeList.size()-1;i++) {
			int from = routeList.get(i)-1;
			int to = routeList.get(i+1)-1;
			//keep track of running cost of the route as you go along the route
			cost += this.graph[from][to];
		}
		//Finally, add the cost of returning back to the starting point/city
		cost += this.graph[finish][start];
		//return total cost of route
		return cost;
	}
	
	/**
	 * Run evolutionary algorithm
	 * @param tsp TSP instance
	 * @param numOfIterations Number of generations produced
	 */
	public void run(Tsp tsp, int numOfIterations) {
		//Fill population with given number of particles
		for(int i = 0; i < tsp.getPopulationSize();i++) {
			tsp.getPopulation().add(tsp.generateRandomRoute());
		}	
		//Calculate fitness of population, then generate new population
		for(int i = 0;i < numOfIterations;i++) {
			System.out.println("ITERATION: " + i);
			for(int j = 0;j < tsp.getPopulationSize();j++) {
				System.out.println(tsp.getPopulation().get(j) + " = " + tsp.getCostOfRoute(tsp.getPopulation().get(j)));
			}
			tsp.calculateFitness();
			tsp.doGeneration(tsp.getPopulation());
			//Print current shortest route
			System.out.println("Record route: " + tsp.getShortestRoute() + " = " + tsp.getLowestCost());
		}
		//Print final answer
		System.out.println("BEST Record route: " + tsp.getShortestRoute() + " = " + tsp.getLowestCost());
	}
	
	public double[][] createGraph() throws FileNotFoundException {
		//Read in the file
		String filename = getFileName();
		File file = new File(filename);
		Scanner inputStream = new Scanner(file);
		List<String> data = new ArrayList<>();
		//keep track of which line is being read in
		int lines = 0;
		//add each line(city and its x and y coordinates) to an ArrayList
		while(inputStream.hasNext()) {
			String value = inputStream.next();
			data.add(value);
			lines++;
		}
		inputStream.close();
		//3D Array that contains each city and its coordinates
		String[][][] values = new String[lines][3][1];
		//this fills the 3D array such that the 3D array contains
		//a bunch of arrays, each of which contains 2 arrays of length 1
		//these 2 arrays contain the x coordinate and the y coordinate of a city
		for(int i = 0;i < data.size();i++) {
			String[] splitData= data.get(i).split(",");
			for(int j = 1;j<3;j++) {
				values[i][0][0] = splitData[1];
				values[i][1][0] = splitData[2];
			}
		}
		//Fill the 2D routes matrix with the distances from one city and another
		//E.g. routes[0][1] contains the distance from city 1 to city 2
		double[][] routes = new double[lines][lines];
		for(int i = 0;i < routes.length;i++) {
			for(int j = 0;j < routes.length;j++) {
				if(i == j) {
					
				} else {
					double city1x = Double.parseDouble(values[i][0][0]);
					double city1y = Double.parseDouble(values[i][1][0]);
					double city2x = Double.parseDouble(values[j][0][0]);
					double city2y = Double.parseDouble(values[j][1][0]);
					double distanceBtwnCities = Math.sqrt(Math.pow((city2y-city1y), 2) + Math.pow((city2x-city1x), 2));
					routes[i][j] = distanceBtwnCities;
				}
			}
		}
		//Return full, complete matrix read for use
		return routes;
	}
	
    // method to sort hashmap by values 
    public static HashMap<List<Integer>, Double> sortByValue(HashMap<List<Integer>, Double> hm) 
    { 
        // Create a list from elements of HashMap 
        List<Map.Entry<List<Integer>, Double> > list = 
               new LinkedList<Map.Entry<List<Integer>, Double> >(hm.entrySet()); 
  
        // Sort the list 
        Collections.sort(list, new Comparator<Map.Entry<List<Integer>, Double> >() { 
            public int compare(Map.Entry<List<Integer>, Double> o1,  
                               Map.Entry<List<Integer>, Double> o2) 
            { 
                return (o1.getValue()).compareTo(o2.getValue()); 
            } 
        }); 
          
        // put data from sorted list to HashMap  
        HashMap<List<Integer>, Double> temp = new LinkedHashMap<List<Integer>, Double>(); 
        for (Map.Entry<List<Integer>, Double> aa : list) { 
            temp.put(aa.getKey(), aa.getValue()); 
        } 
        return temp; 
    } 
}
