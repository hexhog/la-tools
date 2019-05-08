package main

import (
	"bufio"
	"container/list"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"strconv"
	"strings"
)

var workspace *WorkSpace // // static member of the class

// loader section
func loadResponseVector(response *VectorXf, directory string, column string, performLog bool, noise *Noise) {
	// struct dirent *dp;

	files, err := ioutil.ReadDir("directory")
	if err != nil {
		log.Fatal(err)
	}

	var cols, col_i int
	var header string

	rows := response.getLength()

	var tempString string
	var tempFloat float32
	var tempInt int

	// the response count vector
	responseCount := NewVectorXf(rows)

	// grab data vectors
	data = response.getData()
	responseData = responseCount.getData()

	// initialize vectors
	for row_i := 0; row_i < rows; row_i++ {
		data[row_i] = 0
		responseData[row_i] = 0
	}

	for _, file := range files {
		// go through not hidden files
		// if (dp->d_name[0] != '.') {

		file, err := os.Open(directory + "/" + file.Name())
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()

		scanner := bufio.NewScanner(file)
		i := 0
		col_i = -1 // this is the column index we are lookin for
		for scanner.Scan() {
			if i == 0 {
				// read the headers
				tempInt, err = strconv.Atoi(scanner.Text())
				if err != nil {
					log.Fatal(err)
				}
				// make sure the rows match
				if tempInt != rows {
					log.Fatal("Row mismatch (LA vs RE): ", file, "Expected ", rows, " but received ", tempInt)
				}
			} else if i == 1 {
				// find the relevant column
				line := ""
				value := ""
				words := strings.Fields(scanner.Text()) // terminating new line of rows // headers

				for tcol_i, value := range words {
					if value == column {
						col_i = tcol_i
						break
					}
				}
			} else if 2 <= i < 2+rows {
				// ensure we found the relevant column
				if col_i == -1 {
					fmt.Println("Could not find column \"", column, "\" in \"", file, "\"")
				} else {
					// read into the response vector struct
					data := response.getData()
					words := strings.Fields(scanner.Text())

					for tcol_i, value := range words {
						if tcol_i == col_i {
							if value == "" {
								fmt.Println("No value found")
							} else {
								tempFloat, err := strconv.ParseFloat(word, 32)
								if err != nil {
									log.Fatal(err)
								}

								if tempFloat == 0 {
									fmt.Println("Warning: Float value response was 0 in ", file, " in row ", row_i)
								}

								// process the data
								responseData[row_i] += 1
								if performLog {
									data[row_i] += log(tempFloat)
								} else {
									data[row_i] += tempFloat
								}
							}
						}
					}
				}
			}
		}
	}

	// average the responses
	data = response.getData()
	var minResponse float32 = data[0] / responseData[0]
	var maxResponse float32 = data[0] / responseData[0]
	for row_i := 0; row_i < rows; row_i++ {
		data[row_i] /= responseData[row_i]
		//cout << data[row_i] << " from " << responseData[row_i] << " responses " << endl;

		if data[row_i] < minResponse {
			minResponse = data[row_i]
		}
		if data[row_i] > maxResponse {
			maxResponse = data[row_i]
		}
	}

	// add noise as necessary
	var _range float32 = maxResponse - minResponse
	if noise != NULL {
		for row_i := 0; row_i < rows; row_i++ {
			data[row_i] = noise.addNoise(data[row_i], _range)
		}
	}

	// calculate SStot for r-squared calculations
	response.calculateSStot()

	fmt.Println("Loaded responses")
}

// comparison, not case sensitive.
func compareOccurrence(first, second *Occurrence) bool {
	return (first.count < second.count)
}

func allocateOccurrences(occurrence *Occurrence, t int, factors int, occurrenceLists *list) {

	// no more factors or interactions to allocate for
	if factors == 0 || t == 0 {
		// occurrence.list = nil;
		return
	}

	// allocate memory for the next dimension of occurrences
	occurrence.list = make([]Occurrence, factors)

	for factor_i := 0; factor_i < factors; factor_i++ {
		// initialize occurrence
		occurrence.list[factor_i].factorList = make([]int, occurrence.factorList_n+1)
		occurrence.list[factor_i].factorList = make([]int, occurrence.factorList_n+1)
		occurrence.list[factor_i].factorList_n = occurrence.factorList_n + 1
		occurrence.list[factor_i].count = 0
		occurrence.list[factor_i].magnitude = 0

		// populate the factor list
		for factorList_i := 0; factorList_i < occurrence.factorList_n; factorList_i++ {
			occurrence.list[factor_i].factorList[factorList_i] = occurrence.factorList[factorList_i]
		}
		occurrence.list[factor_i].factorList[occurrence.factorList_n] = factor_i

		// push new occurrence into list
		occurrenceLists[occurrence.factorList_n].PushBack(&occurrence.list[factor_i])

		// allocate the next dimension
		allocateOccurrences(&occurrence.list[factor_i], t-1, factor_i, occurrenceLists)
	}
}

func deallocateOccurrences(occurrence *Occurrence, int factors) {

	if occurrence != nil && occurrence.list != nil {

		for factor_i := 0; factor_i < factors; factor_i++ {
			// deallocate list of factors for occurrences
			occurrence.list[factor_i].factorList = nil

			// deallocate the next dimension
			deallocateOccurrences(&occurrence.list[factor_i], factor_i)
		}
		occurrence.list = nil
	}

}

/*
	ANALYSIS Procedure:

	Priority queue begins with 1 model: the model with no terms but an intercept.
	At every iteration, the models are pulled out of the queue, one by one, and for each,
	the top (n) terms are taken, one at a time, based on the distance to residuals.
	For every term, it is added to the current model and r-squared is calculated.
	The new model is then placed in the next priority queue. This process stops
	when all models have as many terms as maxTerms. Each model taken out of the queue
	produces newModels_n new models that are added to the next queue. The queues are
	then priority queues (by model R^2) with a maximum number of models. The same model
	(with the same terms) can be generated in multiple ways, and in this case, duplicates
	will not be added and the function below prints "Duplicate Model!!!".

*/
// struct for linked list of top columns of CS matrix
type ColDetails struct {
	termIndex  int
	dotProduct float32
	used       bool
} // *colDetails;

func createModels(locatingArray *LocatingArray, response *VectorXf, csMatrix *CSMatrix,
	maxTerms int, models_n int, newModels_n int) {
	fmt.Println("Creating Models...")
	setupWorkSpace(response.getLength(), maxTerms)

	// work variables
	var model *Model // current model we are working on

	topModels := make([]*Model, models_n)
	nextTopModels := make([]*Model, models_n)

	// populate initial top models
	topModels[0] = newModel(response, maxTerms, csMatrix)
	for model_i = 1; model_i < models_n; model_i++ {
		topModels[model_i] = nil
	}

	// allocate memory for top columns
	colDetails := make([]ColDetails, csMatrix.getCols())

	for topModels[0] != nil && topModels[0].getTerms() < maxTerms {
		// LOOP HERE

		// make sure all next top models are NULL
		for model_i := 0; model_i < models_n; model_i++ {
			nextTopModels[model_i] = nil
		}

		// grab the models from the topModels priority queue
		for model_i := 0; model_i < models_n; model_i++ {

			// we are done finding the next top models if we hit a NULL model
			if topModels[model_i] == nil {
				break
			}

			// grab the model from top models queue
			model = topModels[model_i]

			//			cout << "Grabbing model" << endl;
			//			model->printModelFactors();
			//			cout << endl;

			// grab the distances to columns in cs matrix
			for col_i := 0; col_i < csMatrix.getCols(); col_i++ {
				colDetails[col_i].dotProduct =
					csMatrix.getProductWithCol(col_i, model.getResiVec())
				colDetails[col_i].termIndex = col_i
				colDetails[col_i].used = model.termExists(col_i)
			}

			// find the columns with the largest dot products (at most as many as we have models)
			for colsUsed := 0; colsUsed < newModels_n; colsUsed++ {

				var largestDotProduct float32
				bestCol_i := -1

				// find the unused column with smallest distance
				for col_i := 0; col_i < csMatrix.getCols(); col_i++ {

					// check if this column has larger dot product and should be marked as best
					if !colDetails[col_i].used &&
						(bestCol_i == -1 || colDetails[col_i].dotProduct > largestDotProduct) {
						bestCol_i = col_i
						largestDotProduct = colDetails[col_i].dotProduct
					}

				}

				// check if all columns have been used already
				if bestCol_i == -1 {
					break
				}

				// mark the term as used
				colDetails[bestCol_i].used = true

				// add the term to the model temporarily
				model.addTerm(bestCol_i)
				model.leastSquares()

				//				cout << "Ran least squares on:" << endl;
				//				model->printModelFactors();
				//				cout << endl;

				// check if the model is a duplicate
				isDuplicate := false
				for model_i := 0; model_i < models_n; model_i++ {
					if nextTopModels[model_i] == nil {
						break
					} else if nextTopModels[model_i].isDuplicate(model) {
						fmt.Println("Duplicate Model!!!")
						isDuplicate = true
						break
					}
				}

				if !isDuplicate {
					// find a possible next top model to replace
					if nextTopModels[models_n-1] == nil ||
						nextTopModels[models_n-1].getRSquared() < model.getRSquared() {

						// make sure we deallocate the older next top model
						if nextTopModels[models_n-1] != nil {
							nextTopModels[models_n-1] = nil
							nextTopModels[models_n-1] = nil
						}

						// insert the new next top model
						nextTopModels[models_n-1] = newModel(model)
					}

					// perform swapping to maintain sorted list
					for model_i := models_n - 2; model_i >= 0; model_i-- {
						if nextTopModels[model_i] == nil ||
							nextTopModels[model_i].getRSquared() < nextTopModels[model_i+1].getRSquared() {

							temp := nextTopModels[model_i]
							nextTopModels[model_i] = nextTopModels[model_i+1]
							nextTopModels[model_i+1] = temp
						}
					}

				}

				// remove the temporarily added term
				model.removeTerm(bestCol_i)

			}

			// delete the model from the top models queue since it has been processed
			// delete model;

		}

		// copy next top models to top models
		for model_i := 0; model_i < models_n; model_i++ {
			topModels[model_i] = nextTopModels[model_i]
		}

		// find the top model
		if topModels[0] != NULL {
			fmt.Println("Top Model (", topModels[0].getRSquared(), "):")
			topModels[0].printModelFactors()
		} else {
			fmt.Println("No Model")
		}

	}

	// delete[] colDetails;

	// count occurrences
	occurrence = *Occurrence{}
	occurrence.factorList = make([]int)
	occurrence.factorList_n = 0
	occurrence.count = 0
	occurrence.magnitude = 0
	// list<Occurrence*> *occurrenceLists = new list<Occurrence*>[locatingArray->getT()];
	occurrenceLists = list.New()
	allocateOccurrences(occurrence, locatingArray.getT(), locatingArray.getFactors(), occurrenceLists)

	fmt.Println("")
	fmt.Println("Final Models Ranking: ")
	for model_i := 0; model_i < models_n; model_i++ {
		if topModels[model_i] != nil {
			fmt.Println("Model ", (model_i + 1), " (", topModels[model_i].getRSquared(), "):")
			topModels[model_i].printModelFactors()
			fmt.Println("")

			topModels[model_i].countOccurrences(occurrence)
			// delete topModels[model_i];
			topModels[model_i] = nil
		} else {
			break
		}
	}

	// sort number of occurrences and print
	fmt.Println("Occurrence Counts")
	for t := 0; t < locatingArray.getT(); t++ {

		occurrenceLists[t].sort(compareOccurrence)
		occurrenceLists[t].reverse()

		// print out the headers
		// cout << setw(10) << right << "Count" << " | " << setw(15) << "Magnitude" << " | " << "Factor Combination" << endl;
		// cout << setw(10) << setfill('-') << "" << " | " <<
		// 		setw(15) << setfill('-') << "" << " | " <<
		// 		setw(20) << setfill('-') << "" << setfill(' ') << endl;
		fmt.Println(right, "Count", " | ", "Magnitude", " | ", "Factor Combination")
		fmt.Println("", " | ", "", " | ", "")

		// print out each count
		for e := occurrenceLists[t].Front(); e != nil; e = e.Next() {
			// only print the count if it has a count
			if e.count > 0 {
				// cout << setw(10) << right << (*it)->count << " | ";
				// cout << setw(15) << right << (*it)->magnitude << " | ";
				fmt.Printf(e.count, " | ")
				fmt.Printf(e.magnitude, " | ")

				// print out the factor combination names
				for factorList_i := 0; factorList_i < e.factorList_n; factorList_i++ {
					if factorList_i != 0 {
						fmt.Printf(" & ")
					}
					fmt.Printf(locatingArray.getFactorData().getFactorName(e.factorList[factorList_i]))
				}
				fmt.Println("")
			}

		}

		fmt.Println("")

	}

	deallocateOccurrences(occurrence, locatingArray.getFactors())
	occurrence.factorList = nil
	occurrence = nil

	occurrenceLists = nil
	topModels = nil
	nextTopModels = nil

}

func main() {

	//LocatingArray *array = new LocatingArray("LA_SMALL.tsv");
	//FactorData *factorData = new FactorData("Factors_SMALL.tsv");
	//VectorXf *response = loadResponseVector("responses_SMALL", "Exposure", false);
	//LocatingArray *array = new LocatingArray("LA_LARGE.tsv");
	//FactorData *factorData = new FactorData("Factors_LARGE.tsv");
	//VectorXf *response = loadResponseVector("responses_LARGE", "Throughput", true);

	// ./Search LA_LARGE.tsv Factors_LARGE.tsv analysis responses_LARGE Throughput 1 13 50 50

	// long long int seed = time(NULL);
	// cout << "Seed:\t" << seed << endl;
	// srand(seed);

	args := os.Args[1:]

	// 	if (argc < 3) {
	// 		cout << "Usage: " << argv[0] << " [LocatingArray.tsv] ([FactorData.tsv]) ..." << endl;
	// 		return 0;
	// 	}

	var noise *Noise
	array := newLocatingArray(args[0], args[1])

	matrix = newCSMatrix(array)

	if args[2] == "memchk" {
		fmt.Println("Check memory and press ENTER")
	} else if args[2] == "analysis" {
		performLog, err := strconv.ParseBool(args[5])
		if err != nil {
			log.Fatal(err)
		}

		terms_n, err := strconv.Atoi(args[6])
		if err != nil {
			log.Fatal(err)
		}

		models_n, err := strconv.Atoi(args[7])
		if err != nil {
			log.Fatal(err)
		}

		newModels_n, err := strconv.Atoi(args[8])
		if err != nil {
			log.Fatal(err)
		}

		response = newVectorXf(array.getTests())
		loadResponseVector(response, args[3], args[4], performLog, noise)
		fmt.Println("Response range: ", response.getData()[0], " to ", response.getData()[response.getLength()-1])

		createModels(array, response, matrix, terms_n, models_n, newModels_n)
		response = nil

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [ResponsesDirectory] [response_column] [1/0 - perform log on responses] [nTerms] [nModels] [nNewModels]" << endl;
	} else if args[2] == "autofind" {
		k, err := strconv.Atoi(args[3])
		if err != nil {
			log.Fatal(err)
		}

		c, err := strconv.Atoi(args[4])
		if err != nil {
			log.Fatal(err)
		}

		startRows, err := strconv.Atoi(args[5])
		if err != nil {
			log.Fatal(err)
		}

		matrix.autoFindRows(k, c, startRows)

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [k Separation] [c Minimum Count] [Start Rows]" << endl;
	} else if args[2] == "fixla" {
		matrix.exactFix()
		array.writeToFile(args[3])

		// cout << " [FixedOutputLA.tsv]" << endl;
	} else if args[2] == "model" {
		responseDir := args[3]
		responseCol := args[4]
		terms, err := strconv.Atoi(args[5])
		if err != nil {
			log.Fatal(err)
		}

		coefficients := make([]float32, terms)
		columns := make([]int, terms)

		i := 5
		for term_i := 0; term_i < terms; term_i++ {
			if i+2 < len(args) {
				coefficients[term_i], err = strconv.ParseFloat(args[i+1], 32)
				if err != nil {
					log.Fatal(err)
				}

				columns[term_i], err = strconv.ParseFloat(args[i+2], 32)
				if err != nil {
					log.Fatal(err)
				}
			} else {
				fmt.Println("Terms are missing")
				coefficients[term_i] = 0
				columns[term_i] = 0
			}
			i += 2
		}

		matrix.writeResponse(responseDir, responseCol, terms, coefficients, columns)

		coefficients = nil
		columns = nil

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [ResponsesDirectory] [response_column] [Terms] [Term 0 coefficient] [Term 0 column] ..." << endl;
	} else if args[2] == "mtfixla" {
		k, err := strconv.Atoi(args[3])
		if err != nil {
			log.Fatal(err)
		}

		c, err := strconv.Atoi(args[4])
		if err != nil {
			log.Fatal(err)
		}

		totalRows, err := strconv.Atoi(args[5])
		if err != nil {
			log.Fatal(err)
		}

		matrix.randomFix(k, c, totalRows)
		array.writeToFile(args[6])

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [k Separation] [c Minimum Count] [Total Rows] [FixedOutputLA.tsv]" << endl;
	} else if args[2] == "sysfixla" {
		k, err := strconv.Atoi(args[3])
		if err != nil {
			log.Fatal(err)
		}

		c, err := strconv.Atoi(args[4])
		if err != nil {
			log.Fatal(err)
		}

		initialRows, err := strconv.Atoi(args[5])
		if err != nil {
			log.Fatal(err)
		}

		minChunk, err := strconv.Atoi(args[6])
		if err != nil {
			log.Fatal(err)
		}

		matrix.systematicRandomFix(k, c, initialRows, minChunk)
		array.writeToFile(args[7])

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [k Separation] [c Minimum Count] [Initial Rows] [Minimum Chunk] [FixedOutputLA.tsv]" << endl;
	} else if args[2] == "checkla" {
		k, err := strconv.Atoi(args[3])
		if err != nil {
			log.Fatal(err)
		}

		c, err := strconv.Atoi(args[4])
		if err != nil {
			log.Fatal(err)
		}

		matrix.performCheck(k, c)

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [k Separation] [c Minimum Count]" << endl;
	} else if args[2] == "noise" {
		ratio, err := strconv.ParseFloat(args[3], 32)
		if err != nil {
			log.Fatal(err)
		}

		noise = newNoise(ratio)

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [ratio]" << endl;
	} else if args[2] == "printcs" {
		fmt.Println("CS Matrix:")
		matrix.print()
	} else if args[2] == "reorderrowsla" {
		k, err := strconv.Atoi(args[3])
		if err != nil {
			log.Fatal(err)
		}

		c, err := strconv.Atoi(args[4])
		if err != nil {
			log.Fatal(err)
		}

		matrix.reorderRows(k, c)
		array.writeToFile(argv[5])

		// cout << "Usage: ... " << argv[arg_i];
		// cout << " [k Separation] [c Minimum Count] [ReorderedOutputLA.tsv]" << endl;
	}

	fmt.Println("")
	fmt.Println("Other Stuff:")
	fmt.Println("Columns in CSMatrix: ", matrix.getCols())

	matrix = nil
	array = nil

	// 	return 0;

}
