package main

import (
	"bufio"
	"container/list"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"time"
	"unsafe"
)

const ENTRY_A = 1
const ENTRY_B = -1

type FactorSetting struct {
	grouped       bool
	factor_i      int8
	index         int8 // level index (1st level if part of group)
	levelsInGroup int8 // levels in group (1 if not a group)
}

type Mapping struct {
	mappedTo int
	mapping  []*Mapping
}

type Path struct {
	min    int
	max    int
	entryA *Path
	entryB *Path
}

type CSCol struct {
	dataVector []float64       // column data
	dataP      []float64       // pointer to the 1st vector element
	factors    int             // number of contributing factors
	setting    []FactorSetting // column headers
	coverable  bool
}

type CSMatrix struct {
	rows          int
	factorData    *FactorData
	locatingArray *LocatingArray
	groupingInfo  []*GroupingInfo
	data          []*CSCol // m by n // vector <CSCol*>*data;

	/* The Mapping structure helps map (factor + level) interactions to their
	appropriate column in the CS Matrix. First, a (factor + level) combination
	is mapped to an index using 'factorLevelMap'. The index for factor = 1,
	level = 2 can be found by: factorLevelMap[1][2]; This index represents a
	(factor + level) combination. A t-way interaction of these indeces is
	mapped to a column of the CS Matrix using 'mapping'. Indeces are placed in
	decreasing order. For example 10, 7, 2 represents a 3-way interaction and
	its column index would be found in the following way:
	mapping->mapping[10]->mapping[7]->mapping[2]->mappedTo;
	*/
	factorLevelMap [][]int
	mapping        *Mapping
}

func newCSMatrix(locatingArray *LocatingArray) *CSMatrix {
	csm := &CSMatrix{}

	csm.locatingArray = locatingArray

	// assign factor data variable for use with grabbing column names
	csm.factorData = locatingArray.getFactorData()

	// to be used when populating CSMatrix
	var csCol *CSCol
	var sum float64

	// get number of factors in locating array
	factors := locatingArray.getFactors()

	// get level counts (how many levels for each factor)
	csm.groupingInfo = locatingArray.getGroupingInfo()

	// get the level matrix from locating array (this is the main data)
	levelMatrix := locatingArray.getLevelMatrix()

	// a row for every test
	csm.rows = locatingArray.getTests()

	// initialize the data column vector
	csm.data = make([]*CSCol, 0)

	// initialize (factor + level)-to-index map
	csm.factorLevelMap = make([][]int, factors)
	for factor_i := 0; factor_i < factors; factor_i++ {
		csm.factorLevelMap[factor_i] = make([]int, csm.groupingInfo[factor_i].levels)
	}

	// initialize index-to-column map (to be populated later)
	csm.mapping = &Mapping{}
	csm.mapping.mappedTo = 0

	// sum of squares of each column for normalization
	sumOfSquares := make([]float64, 0)

	// populate compressive sensing matrix column by column
	col_i := 0

	// Add the INTERCEPT
	csCol = &CSCol{}
	sum = 0

	csCol.factors = 0
	csCol.setting = make([]FactorSetting, 0)

	// populate the intercept
	for row_i := 0; row_i < csm.rows; row_i++ {
		csm.addArray(csCol)
		csCol.dataP[row_i] = ENTRY_A
		sum += 1
	}

	// push into vector
	csm.data = append(csm.data, csCol)
	sumOfSquares = append(sumOfSquares, sum)
	col_i++

	// go over all 1-way interactions
	for factor_i := 0; factor_i < factors; factor_i++ {

		for level_i := 0; level_i < csm.groupingInfo[factor_i].levels; level_i++ {

			csm.addOneWayInteraction(factor_i, byte(level_i), levelMatrix, &sumOfSquares)

			// set factor level map
			csm.factorLevelMap[factor_i][level_i] = col_i - 1 // subtract 1 for INTERCEPT

			// increment column index
			col_i++
		}

	}

	lastOneWay_i := col_i

	// initialize 1st level mapping
	csm.mapping.mapping = make([]*Mapping, col_i-1)

	fmt.Println("Adding t-way interactions")
	csm.addTWayInteractions(csCol, col_i-1, &col_i, csm.locatingArray.getT(),
		csm.mapping.mapping, &sumOfSquares, csm.groupingInfo, levelMatrix)

	fmt.Println("Went over ", col_i, " columns")

	// check coverability
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		csCol = csm.data[col_i]
		csCol.coverable = csm.checkColumnCoverability(csCol)

		if !csCol.coverable {
			fmt.Println("Not coverable: ", csm.getColName(csCol))
		}
	}

	// perform sqaure roots on sum of squares
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		sumOfSquares[col_i] = math.Sqrt(sumOfSquares[col_i])
	}

	// normalize all columns
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		csCol = csm.data[col_i]
		for row_i := 0; row_i < csm.getRows(); row_i++ {
			//	csCol->data[row_i] = csCol->data[row_i] / sumOfSquares[col_i];
		}
		csm.data[col_i] = csCol
	}

	fmt.Println("Finished constructing CS Matrix")

	return csm
}

func (csm *CSMatrix) addArray(csCol *CSCol) {
	csCol.dataVector = append(csCol.dataVector, 0)
	csCol.dataP = csCol.dataVector
}

func (csm *CSMatrix) remArray(csCol *CSCol) {
	csCol.dataVector = csCol.dataVector[:len(csCol.dataVector)-1]
	csCol.dataP = csCol.dataVector
}

func (csm *CSMatrix) addTWayInteractions(csColA *CSCol, colBMax_i int, col_i *int, t int,
	mapping []*Mapping, sumOfSquares *[]float64, groupingInfo []*GroupingInfo, levelMatrix [][]int8) {
	var csCol, csColB, csColC *CSCol
	var sum float64

	// the offset is 1 for the INTERCEPT
	colBOffset := 1

	colCMax_i := 0

	// start after INTERCEPT, iterate over all main effects until colBMax_i
	for colB_i := 0; colB_i < colBMax_i; colB_i++ {

		// grab pointer to column
		csColB = csm.data[colB_i+colBOffset]

		// the next 2 lines move colCMax_i to the 1st column with the same factor as csColB
		csColC = csm.data[colCMax_i+colBOffset]
		if csColB.setting[0].factor_i > csColC.setting[0].factor_i {
			colCMax_i = colB_i
		}

		// set mapping for this column
		mapping[colB_i] = &Mapping{}
		if t > 1 {
			mapping[colB_i].mapping = make([]*Mapping, colCMax_i)
		}
		if csColA.factors > 0 {
			mapping[colB_i].mappedTo = *col_i
		} else {
			mapping[colB_i].mappedTo = colB_i + colBOffset
		}
		groupMapping := mapping[colB_i]

		var colsInGroupB int8 = 1
		var groupIndexB int8 = -1
		levelIndexB := csColB.setting[0].index

		// check if this column is grouped
		if groupingInfo[csColB.setting[0].factor_i].grouped {
			groupIndexB = groupingInfo[csColB.setting[0].factor_i].levelGroups[levelIndexB]

			for level_i := levelIndexB + 1; int(level_i) < groupingInfo[csColB.setting[0].factor_i].levels &&
				groupingInfo[csColB.setting[0].factor_i].levelGroups[level_i] == groupIndexB; level_i++ {

				colsInGroupB++
				colB_i++

				// set to group mapping unless main effect (no grouping on main effects)
				if csColA.factors > 0 {
					mapping[colB_i] = groupMapping
				} else {
					mapping[colB_i] = &Mapping{}
					mapping[colB_i].mapping = groupMapping.mapping
					mapping[colB_i].mappedTo = colB_i + colBOffset
				}
			}

			csColB = csm.data[colB_i+colBOffset]
		}

		// create new column for CS Matrix
		csCol = &CSCol{}
		for row_i := 0; row_i < csm.rows; row_i++ {
			csm.addArray(csCol)
		}

		// set the headers from the combining columns
		csCol.factors = csColA.factors + 1
		csCol.setting = make([]FactorSetting, csCol.factors)

		// copy previous factors
		for setting_i := 0; setting_i < csColA.factors; setting_i++ {
			csCol.setting[setting_i] = csColA.setting[setting_i]
		}

		// assign new setting
		csCol.setting[csColA.factors].grouped = (groupIndexB != -1)         // is 1st factor grouped
		csCol.setting[csColA.factors].factor_i = csColB.setting[0].factor_i // 1st factor
		csCol.setting[csColA.factors].index = levelIndexB                   // set level of 1st factor
		csCol.setting[csColA.factors].levelsInGroup = colsInGroupB          // set last level if group

		// populate the actual column data of CS matrix
		sum = csm.populateColumnData(csCol, levelMatrix, 0, csm.rows)

		// add new CS column to matrix if needed
		colAddedToMatrix := false
		if csCol.factors > 1 {
			// push into vector
			csm.data = append(csm.data, csCol)
			*sumOfSquares = append(*sumOfSquares, sum)

			*col_i++

			colAddedToMatrix = true
		}

		if t > 1 {
			// recursive call
			csm.addTWayInteractions(csCol, colCMax_i, col_i, t-1,
				mapping[colB_i].mapping, sumOfSquares, groupingInfo, levelMatrix)
		}

	}
}

func (csm *CSMatrix) populateColumnData(csCol *CSCol, levelMatrix [][]int8, row_top int, row_len int) float64 {
	var sum float64

	// populate every row
	var rowData bool // start with true, perform AND operation
	for row_i := row_top; row_i < row_top+row_len; row_i++ {

		rowData = true // start with true, perform AND operation

		for setting_i := 0; setting_i < csCol.factors; setting_i++ {

			rowData = rowData && (levelMatrix[row_i][csCol.setting[setting_i].factor_i] >=
				csCol.setting[setting_i].index)

			rowData = rowData && levelMatrix[row_i][csCol.setting[setting_i].factor_i] <
				csCol.setting[setting_i].index+csCol.setting[setting_i].levelsInGroup

		}

		// AND operation
		if rowData {
			csCol.dataP[row_i] = ENTRY_A
		} else {
			csCol.dataP[row_i] = ENTRY_B
		}

		// add to sum of squares
		sum += csCol.dataP[row_i] * csCol.dataP[row_i]
	}

	return sum
}

func (csm *CSMatrix) reorderRows(k, c int) {
	var settingToResample []*FactorSetting
	var nPaths int
	var score int64
	cols := csm.getCols()

	// check advanced
	array := make([]*CSCol, cols)
	for col_i := 0; col_i < cols; col_i++ {
		array[col_i] = csm.data[col_i]
	}

	coverableMin := csm.sortByCoverable(array, 0, csm.getCols()-1)
	fmt.Println("Coverable columns begin at: ", coverableMin)

	tWayMin := csm.sortByTWayInteraction(array, coverableMin, csm.getCols()-1)
	fmt.Println("t-way interactions begin at: ", tWayMin)

	rowContributions := make([]int64, csm.rows)

	path := &Path{}
	path.min = coverableMin
	path.max = csm.getCols() - 1

	for {
		nPaths = 0
		for row_i := 0; row_i < csm.rows; row_i++ {
			rowContributions[row_i] = 0
		}
		csm.pathSort(array, path, 0, nPaths, nil)

		score = 0
		var settingToResample *FactorSetting
		csm.pathLAChecker(array, path, path, 0, k, &score, &settingToResample, rowContributions)
		csm.minCountCheck(array, c, &score, &settingToResample, rowContributions)

		fmt.Println("Score: ", score)

		swaps := 0
		for {

			// find the max contribution out of place
			outOrderRow_i := -1
			for row_i := 1; row_i < csm.rows; row_i++ {
				if rowContributions[row_i-1] < rowContributions[row_i] {
					outOrderRow_i = row_i
					break
				}
			}

			if outOrderRow_i != -1 {
				row_i2 := outOrderRow_i
				for row_i := row_i2 + 1; row_i < csm.rows; row_i++ {
					if rowContributions[row_i] >= rowContributions[row_i2] {
						row_i2 = row_i
					}
				}

				row_i1 := -1
				for row_i := 0; row_i < csm.rows; row_i++ {
					if rowContributions[row_i2] > rowContributions[row_i] {
						row_i1 = row_i
						break
					}
				}

				// perform swap if necessary
				csm.swapRows(array, row_i1, row_i2)
				tempContribution := rowContributions[row_i1]
				rowContributions[row_i1] = rowContributions[row_i2]
				rowContributions[row_i2] = tempContribution
				swaps++
			} else {
				break
			}

		}

		fmt.Println("\tSwaps: ", swaps)
		if swaps == 0 {
			break
		}
	}

	for row_i := 0; row_i < csm.rows; row_i++ {
		fmt.Println(row_i, "\t", rowContributions[row_i])
	}

}

func (csm *CSMatrix) minCountCheck(array []*CSCol, c int, score *int64, settingToResample **FactorSetting,
	rowContributions []int64) {
	cols := csm.getCols()

	count := make([]int, cols)

	for col_i := 0; col_i < cols; col_i++ {
		count[col_i] = 0

		for row_i := 0; row_i < csm.rows && count[col_i] < c; row_i++ {
			if array[col_i].coverable && array[col_i].dataP[row_i] == ENTRY_A {
				count[col_i]++

				// add row contributions
				if rowContributions != nil {
					rowContributions[row_i]++
				}
			}
		}
	}

	for col_i := 0; col_i < cols; col_i++ {
		if array[col_i].coverable && count[col_i] < c {
			//			cout << "Below c requirement (" << count[col_i] << "): " << getColName(array[col_i]) << endl;
			*score += int64(c - count[col_i])
			if *settingToResample == nil {
				// randomly choose a setting in the column to resample
				*settingToResample = &array[col_i].setting[rand.Intn(array[col_i].factors)]
			}
		}
	}

}

func (csm *CSMatrix) exactFix() {
	// create a work array
	array := make([]*CSCol, csm.getCols())
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		array[col_i] = csm.data[col_i]
	}

	var score int64

	csm.smartSort(array, 0)
	score = csm.getArrayScore(array)

	fmt.Println("Original linear LA Score: ", score)

	if csm.locatingArray.getNConGroups() == 0 {
		for score > 0 {
			csm.addRowFix(array, &score)
		}

		fmt.Println("Complete LA created with score: ", score)
		fmt.Println("Rows: ", csm.getRows())
	} else {
		fmt.Println("Unable to perform fixla because constraints were found. Remove the constraints to perform this operation!")
	}

}

func (csm *CSMatrix) performCheck(k, c int) {
	nConGroups := csm.locatingArray.getNConGroups()
	conGroups := csm.locatingArray.getConGroups()
	var score int64

	// create a work array
	array := make([]*CSCol, csm.getCols())
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		array[col_i] = csm.data[col_i]
	}

	coverableMin := csm.sortByCoverable(array, 0, csm.getCols()-1)
	fmt.Println("Coverable columns begin at: ", coverableMin)

	tWayMin := csm.sortByTWayInteraction(array, coverableMin, csm.getCols()-1)
	fmt.Println("t-way interactions begin at: ", tWayMin)

	path := &Path{}
	path.min = coverableMin
	path.max = csm.getCols() - 1

	pathList := list.New()

	nPaths := 0

	// grab initial time
	start := time.Now()
	csm.pathSort(array, path, 0, nPaths, &pathList)
	// check current time
	elapsedTime := time.Since(start).Seconds()
	fmt.Println("Elapsed for path sort: ", elapsedTime)

	fmt.Println("Memory check: nPaths = ", nPaths, " of size ", unsafe.Sizeof(path))
	fmt.Println("Unfinished paths: ", pathList.Len())

	// grab initial time
	start = time.Now()
	var settingToResample *FactorSetting
	csm.pathLAChecker(array, path, path, 0, k, &score, &settingToResample, nil)
	csm.minCountCheck(array, c, &score, &settingToResample, nil)
	// check current time
	elapsedTime = time.Since(start).Seconds()

	fmt.Println("Path and min count LA Score: ", score)
	fmt.Println("Elapsed for path and min count check: ", elapsedTime)

	csm.deletePath(path)

	// grab initial time
	start = time.Now()
	csm.smartSort(array, 0)
	// check current time
	elapsedTime = time.Since(start).Seconds()

	fmt.Println("Elapsed for smart sort: ", elapsedTime)

	// grab initial time
	start = time.Now()
	score = csm.getArrayScore(array)
	csm.minCountCheck(array, c, &score, &settingToResample, nil)
	// check current time
	elapsedTime = time.Since(start).Seconds()

	fmt.Println("Weird linear check LA Score (should not match other scores): ", score)
	fmt.Println("Elapsed for linear check: ", elapsedTime)

	/* BRUTE FORCE */
	fmt.Println("Performing brute force check... this could take awhile")
	// grab initial time
	start = time.Now()
	score = csm.getBruteForceArrayScore(array, k)
	csm.minCountCheck(array, c, &score, &settingToResample, nil)
	// check current time
	elapsedTime = time.Since(start).Seconds()

	fmt.Println("Brute force and min count LA Score (should match path score): ", score)
	fmt.Println("Elapsed for brute force check: ", elapsedTime)

	fmt.Println("Checking for constraint violations...")
	for row_i := 0; row_i < csm.rows; row_i++ {
		for conGroup_i := 0; conGroup_i < nConGroups; conGroup_i++ {
			if !conGroups[conGroup_i].getResult(row_i) {
				fmt.Println("Constraint group ", conGroup_i, " violated in row ", row_i)
			}
		}

	}
	fmt.Println("Done!")

}

func (csm *CSMatrix) autoFindRows(k, c, startRows int) {
	iters := 1000
	cols := csm.getCols()

	upperBound := startRows
	lowerBound := 1

	// check advanced
	array := make([]*CSCol, csm.getCols())
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		array[col_i] = csm.data[col_i]
	}

	var twoWayMin int
	for twoWayMin := 0; twoWayMin < cols; twoWayMin++ {
		if array[twoWayMin].factors > 1 {
			break
		}
	}
	fmt.Println("Two-Way Min: ", twoWayMin)

	factors := csm.locatingArray.getFactors()
	nPaths := 0
	var score int64
	var settingToResample *FactorSetting

	path := &Path{}
	//	path.min = twoWayMin;
	path.min = 0
	path.max = csm.getCols() - 1

	// add more rows to reach total count
	csm.resizeArray(array, startRows)

	// use a binary search to find the correct value
	for {

		// check if it finds a proper array once in 5 times
		testPassed := false
		for i := 0; i < 5; i++ {

			csm.randomizeArray(array)

			pathList := list.New()
			pathList.PushFront(path)

			score = 0
			csm.randomizePaths(array, &settingToResample, path, 0, k, c, &score, &pathList, iters)

			fmt.Println("Score: ", score)

			if settingToResample == nil {
				testPassed = true
				break
			} else if score > 100 {
				// do not try 5 times because the score is greater than 100
				break
			}
		}

		// reset upper / lower bounds
		if testPassed {
			upperBound = rows
		} else {
			lowerBound = rows + 1
		}

		// check if the upper and lower bounds match
		if upperBound < lowerBound {
			fmt.Println("Our bounds messed up :(")
		}
		if upperBound <= lowerBound {
			break
		}

		// calculate a median row count
		medRowCount := 1 * int((lowerBound+upperBound)/2.0)

		// resize the array
		csm.resizeArray(array, medRowCount)

	}

	fmt.Println("Finished with array containing (", lowerBound, ":", upperBound, ") rows!")

	csm.deletePath(path)
}

func (csm *CSMatrix) randomFix(k, c, totalRows int) {
	iters := 1000
	cols := csm.getCols()

	// check advanced
	array = make([]*CSCol, cols)
	for col_i := 0; col_i < cols; col_i++ {
		array[col_i] = csm.data[col_i]
	}

	coverableMin := csm.sortByCoverable(array, 0, cols-1)
	fmt.Println("Coverable columns begin at: ", coverableMin)

	tWayMin := csm.sortByTWayInteraction(array, coverableMin, cols-1)
	fmt.Println("t-way interactions begin at: ", tWayMin)

	factors := csm.locatingArray.getFactors()
	nPaths := 0
	var score int64
	settingToResample = *FactorSetting

	path := &Path{}
	path.min = coverableMin
	path.max = csm.getCols() - 1

	// add more rows to reach total count
	csm.resizeArray(array, totalRows)

	pathList := list.New()
	pathList.PushFront(path)

	score = 0
	csm.randomizePaths(array, settingToResample, path, 0, k, c, score, &pathList, iters)

	csm.minCountCheck(array, c, score, settingToResample, NULL)

	fmt.Println("Score: ", score)

	csm.deletePath(path)
}

func (csm *CSMatrix) systematicRandomFix(k, c, initialRows, minChunk int) {
	chunk := initialRows
	finalizedRows := rows
	totalRows := (finalizedRows + chunk)
	cols := csm.getCols()

	// check advanced
	array = make([]*CSCol, cols)
	for col_i := 0; col_i < cols; col_i++ {
		array[col_i] = csm.data[col_i]
	}

	coverableMin := csm.sortByCoverable(array, 0, cols-1)
	fmt.Println("Coverable columns begin at: ", coverableMin)

	tWayMin := csm.sortByTWayInteraction(array, coverableMin, cols-1)
	fmt.Println("t-way interactions begin at: ", tWayMin)

	path = &Path{}
	path.min = coverableMin
	path.max = csm.getCols() - 1

	pathList := list.New()

	nPaths := 0
	csm.pathSort(array, path, 0, nPaths, &pathList)
	fmt.Println("nPaths: ", nPaths, " of size ", sizeof(Path))
	fmt.Println("Unfinished paths: ", pathList.size())

	var score int64

	// grab initial time
	start := time.Now()

	settingToResample * FactorSetting
	csm.pathLAChecker(array, path, path, 0, k, score, settingToResample, nil)
	csm.minCountCheck(array, c, score, settingToResample, nil)

	// check current time
	elapsedTime := time.Since(start).Seconds()
	fmt.Println("Elapsed: ", elapsedTime)

	fmt.Println("Score: ", score)
	factors := csm.locatingArray.getFactors()

	for settingToResample != nil {
		csm.resizeArray(array, totalRows)

		csm.randomizePaths(array, settingToResample, path, finalizedRows, k, c, score, &pathList, 1000)
		finalizedRows = rows

		chunk -= chunk / 2
		if chunk < minChunk {
			chunk = minChunk
		}
		totalRows += chunk
	}

}

func (csm *CSMatrix) randomizePaths(array []*CSCol, settingToResample **FactorSetting, path *Path, row_top int,
	k int, c int, score *int64, pathList *List, iters int) {
	cols := getCols()
	var newScore int64
	var factor_i, nPaths int
	var conGroup *ConstraintGroup
	levelMatrix := csm.locatingArray.getLevelMatrix()
	factors := csm.locatingArray.getFactors()

	// custom factors to resample
	nCustomFactors := 0
	customFactorIndeces = make([]int, nCustomFactors)
	//	customFactorIndeces[0] = 16;
	//	customFactorIndeces[1] = 17;

	// allocate memory for saving old levels of locating array
	oldLevels = make([][]byte, rows)
	for row_i := 0; row_i < rows; row_i++ {
		oldLevels[row_i] = make([]byte, factors)
	}

	// sort paths
	for e := pathList.Front(); e != nil; e = e.Next() {
		nPaths := 0
		csm.pathSort(array, *e, row_top, nPaths, nil)
	}

	// run initial checker
	*score = 0
	settingToResample = nil
	csm.pathLAChecker(array, path, path, 0, k, score, settingToResample, nil)
	csm.minCountCheck(array, c, score, settingToResample, NULL)
	fmt.Println("Score: ", score, endl)

	for iter := 0; iter < iters && score > 0; iter++ {

		// ensure we recieved an actual setting
		if settingToResample == nil {
			fmt.Println("No resampleable setting was found", endl)
			break
		}

		// get factors to resample (all factors in constraint group if one exists)
		conGroup = csm.groupingInfo[settingToResample.factor_i].conGroup
		if conGroup == nil {
			// get factor to resample
			factor_i = settingToResample.factor_i

			// resample locating array
			for row_i := row_top; row_i < csm.rows; row_i++ {
				oldLevels[row_i][factor_i] = levelMatrix[row_i][factor_i]
				if rand.Intn(100) < 100 {
					levelMatrix[row_i][factor_i] = rand.Intn(csm.groupingInfo[factor_i].levels)
				}
			}

			// repopulate columns of CS matrix
			for level_i := 0; level_i < csm.groupingInfo[factor_i].levels; level_i++ {
				csm.repopulateColumns2(factor_i, level_i, row_top, rows-row_top)
			}
		} else if nCustomFactors > 0 {
			for row_i := row_top; row_i < csm.rows; row_i++ {
				for factor_i := 0; factor_i < nCustomFactors; factor_i++ {
					oldLevels[row_i][customFactorIndeces[factor_i]] = levelMatrix[row_i][customFactorIndeces[factor_i]]
				}
				if rand.Intn(100) < 100 {
					for factor_i := 0; factor_i < nCustomFactors; factor_i++ {
						levelMatrix[row_i][customFactorIndeces[factor_i]] = rand.Intn(groupingInfo[customFactorIndeces[factor_i]].levels)
					}
				}
			}

			// repopulate columns of CS matrix for every factor in constraint group and for each level
			for factor_i := 0; factor_i < nCustomFactors; factor_i++ {
				for level_i := 0; level_i < csm.groupingInfo[customFactorIndeces[factor_i]].levels; level_i++ {
					csm.repopulateColumns2(customFactorIndeces[factor_i], level_i, row_top, rows-row_top)
				}
			}
		} else {
			for row_i := row_top; row_i < csm.rows; row_i++ {
				for factor_i := 0; factor_i < conGroup.factors; factor_i++ {
					oldLevels[row_i][conGroup.factorIndeces[factor_i]] = levelMatrix[row_i][conGroup.factorIndeces[factor_i]]
				}
				if rand.Intn(100) < 100 {
					conGroup.randPopulateLevelRow(levelMatrix[row_i])
				}
			}

			// repopulate columns of CS matrix for every factor in constraint group and for each level
			for factor_i := 0; factor_i < conGroup.factors; factor_i++ {
				for level_i := 0; level_i < groupingInfo[conGroup.factorIndeces[factor_i]].levels; level_i++ {
					csm.repopulateColumns2(conGroup.factorIndeces[factor_i], level_i, row_top, rows-row_top)
				}
			}
		}

		// grab initial time
		start := time.Now()
		// sort paths and recheck score
		for e := pathList.Front(); e != nil; e = e.Next() {
			nPaths := 0
			csm.pathSort(array, *e, row_top, nPaths, nil)
		}

		elapsedTime := time.Since(start).Seconds()
		//		cout << "Elapsed After Sort: " << elapsedTime << endl;

		newScore = 0
		newSettingToResample = *FactorSetting

		// grab initial time
		start = time.Now()
		csm.pathLAChecker(array, path, path, 0, k, newScore, newSettingToResample, nil)
		csm.minCountCheck(array, c, newScore, newSettingToResample, nil)
		elapsedTime = time.Since(start).Seconds()
		//		cout << "Elapsed After Checker: " << elapsedTime << endl;

		if newScore <= score { // add "|| true" to cause every change to be implemented, not just improving changes
			fmt.Println("Rows: ", rows, " Iter: ", iter, ": ")
			fmt.Println(newScore, ": \t", score, " \tAccepted ")
			settingToResample = newSettingToResample
			score = newScore
		} else {
			fmt.Println("Rows: ", rows, " Iter: ", iter, ": ")
			fmt.Println(score, " \tMaintained ")

			if conGroup == nil {
				// get factor to resample
				factor_i = settingToResample.factor_i

				// rollback the change
				for row_i := row_top; row_i < csm.rows; row_i++ {
					levelMatrix[row_i][factor_i] = oldLevels[row_i][factor_i]
				}

				// repopulate columns of CS matrix
				for level_i = 0; level_i < csm.groupingInfo[factor_i].levels; level_i++ {
					csm.repopulateColumns2(factor_i, level_i, row_top, rows-row_top)
				}
			} else if nCustomFactors > 0 {
				for row_i := row_top; row_i < csm.rows; row_i++ {
					for factor_i := 0; factor_i < nCustomFactors; factor_i++ {
						levelMatrix[row_i][customFactorIndeces[factor_i]] = oldLevels[row_i][customFactorIndeces[factor_i]]
					}
				}

				// repopulate columns of CS matrix for every factor in constraint group and for each level
				for factor_i := 0; factor_i < nCustomFactors; factor_i++ {
					for level_i := 0; level_i < csm.groupingInfo[customFactorIndeces[factor_i]].levels; level_i++ {
						csm.repopulateColumns2(customFactorIndeces[factor_i], level_i, row_top, rows-row_top)
					}
				}
			} else {
				// rollback the change
				for row_i := row_top; row_i < rows; row_i++ {
					for factor_i := 0; factor_i < conGroup.factors; factor_i++ {
						levelMatrix[row_i][conGroup.factorIndeces[factor_i]] = oldLevels[row_i][conGroup.factorIndeces[factor_i]]
					}
				}

				// repopulate columns of CS matrix for every factor in constraint group and for each level
				for factor_i := 0; factor_i < conGroup.factors; factor_i++ {
					for level_i := 0; level_i < groupingInfo[conGroup.factorIndeces[factor_i]].levels; level_i++ {
						csm.repopulateColumns2(conGroup.factorIndeces[factor_i], level_i, row_top, rows-row_top)
					}
				}
			}
			//			cout << score << ": \t" << newScore << " \tRejected" << endl;
		}
	}

	// perform one last sort and save final unfinished paths to list
	newPathList := List.New() // list <Path*>newPathList(*pathList);
	pathList.init()           // newPathList

	// sort paths and recheck score
	for e := newPathList.Front(); e != nil; e = e.Next() {
		nPaths := 0
		csm.pathSort(array, *e, row_top, nPaths, pathList)
	}

	*score = 0
	settingToResample = nil
	csm.pathLAChecker(array, path, path, 0, k, score, settingToResample, nil)
	csm.minCountCheck(array, c, score, settingToResample, nil)
}

// LEGACY
func (csm *CSMatrix) randomizeRows(backupArray []*CSCol, array []*CSCol, csScore *int64, row_top int, row_len int) {
	cols := csm.getCols()
	var newCsScore int64
	var factor_i, resampleFactor int

	var settingToResample *FactorSetting

	levelMatrix := csm.locatingArray.getLevelMatrix()

	oldLevels := make([]byte, rows)

	csm.smartSort(array, row_top)
	csScore = csm.getArrayScore(array)

	fmt.Println("Score: ", csScore)

	for iter := 0; iter < 1000; {

		if csScore <= 0 {
			break
		}

		for col_i := 0; col_i < cols-1; col_i++ {

			// check if the streak ended
			if csm.compare(array[col_i], array[col_i+1], 0, rows) == 0 {

				// copy to backup array
				copy(backupArray, array)

				for iter < 1000 {
					iter++

					resampleFactor = rand.Intn(array[col_i].factors + array[col_i+1].factors)

					if resampleFactor < array[col_i].factors {
						factor_i = array[col_i].setting[resampleFactor].factor_i
					} else {
						factor_i = array[col_i+1].setting[resampleFactor-array[col_i].factors].factor_i
					}

					for row_i := row_top; row_i < row_top+row_len; row_i++ {
						oldLevels[row_i] = levelMatrix[row_i][factor_i]
						if rand.Intn(100) < 100 {
							levelMatrix[row_i][factor_i] = rand.Intn(csm.groupingInfo[factor_i].levels)
						}
					}

					for level_i := 0; level_i < csm.groupingInfo[factor_i].levels; level_i++ {
						csm.repopulateColumns2(factor_i, level_i, row_top, row_len)
					}

					csm.smartSort(array, row_top)
					newCsScore = csm.getArrayScore(array)

					var likelihood float64 = 10 / Pow(newCsScore/csScore, 10)
					fmt.Println("Rows: ", rows, " Iter: ", iter, ": ")
					if newCsScore <= csScore { // || rand() % 100 < likelihood) {
						fmt.Println(newCsScore, ": \t", csScore, " \tAccepted ", likelihood, "%")
						csScore = newCsScore
						break
					} else {
						// rollback the change
						for row_i := row_top; row_i < row_top+row_len; row_i++ {
							levelMatrix[row_i][factor_i] = oldLevels[row_i]
						}

						for level_i := 0; level_i < groupingInfo[factor_i].levels; level_i++ {
							csm.repopulateColumns2(factor_i, level_i, row_top, row_len)
						}

						// restore order from backup array
						copy(array, backupArray)

						fmt.Println(csScore, ": \t", newCsScore, " \tRejected")
					}
				}

				break

			} else if csm.compare(array[col_i], array[col_i+1], 0, rows) > 0 {
				fmt.Println("Mistake in array")
			}
		}

	}

}

func (csm *CSMatrix) repopulateColumns2(setFactor_i int, setLevel_i int, row_top int, row_len int) {
	lastCol_i := -1
	csm.repopulateColumns(setFactor_i, setLevel_i, csm.locatingArray.getFactors()-1, csm.locatingArray.getT(),
		csm.mapping, csm.locatingArray.getLevelMatrix(), &lastCol_i, row_top, row_len)
}

func (csm *CSMatrix) repopulateColumns(setFactor_i int, setLevel_i int, maxFactor_i int, t int,
	mapping *Mapping, levelMatrix [][]byte, lastCol_i *int, row_top int, row_len int) {

	if setFactor_i > maxFactor_i && mapping.mappedTo != *lastCol_i {
		populateColumnData(data[mapping.mappedTo], levelMatrix, row_top, row_len)

		lastCol_i = mapping.mappedTo
	}

	if t == 0 {
		return
	}

	var minFactor_i int
	var minLevel_i, maxLevel_i int

	if setFactor_i > maxFactor_i {
		minFactor_i = t - 1
	} else {
		minFactor_i = setFactor_i

		// if this is the last interaction, make sure setFactor shows up
		if t == 1 {
			maxFactor_i = minFactor_i
		}
	}

	for factor_i := minFactor_i; factor_i <= maxFactor_i; factor_i++ {

		if factor_i == setFactor_i {
			csm.repopulateColumns(setFactor_i, setLevel_i, factor_i-1, t-1,
				mapping.mapping[factorLevelMap[factor_i][setLevel_i]], levelMatrix, lastCol_i, row_top, row_len)
		} else {
			for level_i := 0; level_i < groupingInfo[factor_i].levels; level_i++ {
				csm.repopulateColumns(setFactor_i, setLevel_i, factor_i-1, t-1,
					mapping.mapping[factorLevelMap[factor_i][level_i]], levelMatrix, lastCol_i, row_top, row_len)
			}
		}

	}
}
func (csm *CSMatrix) getRows() int {
	return csm.rows
}

func (csm *CSMatrix) getCols() int {
	return csm.data.size()
}

func (csm *CSMatrix) getCol(col_i int) *CSCol {
	return csm.data[col_i]
}

func (csm *CSMatrix) getDistanceToCol(col_i int, residuals *float64) float64 {
	var distanceSum, subtractionResult float64
	csCol := data[col_i]

	for row_i := 0; row_i < csm.rows; row_i++ {
		subtractionResult = csCol.dataP[row_i] - residuals[row_i]
		distanceSum += subtractionResult * subtractionResult
	}

	return distanceSum
}

func (csm *CSMatrix) getProductWithCol(col_i int, residuals *float64) float64 {
	var dotSum float64
	csCol := data[col_i]

	for row_i := 0; row_i < csm.rows; row_i++ {
		dotSum += csCol.dataP[row_i] * residuals[row_i]
	}

	return Abs(dotSum)
}

func (csm *CSMatrix) getColName(csCol *CSCol) string {
	colName = ""

	if csCol.factors == 0 {
		colName += "INTERCEPT"
	}

	// add all factor strings
	for setting_i := 0; setting_i < csCol.factors; setting_i++ {
		if setting_i != 0 {
			colName += " & "
		}
		colName += getFactorString(csCol.setting[setting_i])
	}

	return colName
}

func (csm *CSMatrix) getFactorString(setting FactorSetting) string {
	factorString = ""

	factorString += factorData.getFactorName(setting.factor_i)
	factorString += "="

	if setting.grouped {
		factorString += "GROUP("
	}

	for index := setting.index; index < setting.index+setting.levelsInGroup; index++ {
		if index != setting.index {
			factorString += "|"
		}
		factorString += getFactorLevelName(setting.factor_i, index)
	}

	if setting.grouped {
		factorString += ")"
	}

	return factorString
}

func (csm *CSMatrix) getFactorLevelName(factor_i, level_i int) string {
	return csm.factorData.getFactorLevelName(factor_i, level_i)
}

func (csm *CSMatrix) addOneWayInteraction(factor_i int, level_i byte, levelMatrix [][]int8, sumOfSquares *[]float64) {
	// create new column for CS Matrix
	csCol := &CSCol{}
	sum := 0

	// assign the headers
	csCol.factors = 1 // 1 contributing factor (1-way interaction)
	csCol.setting = make([]FactorSetting, 1)
	csCol.setting[0].grouped = false
	csCol.setting[0].factor_i = factor_i // single factor
	csCol.setting[0].index = level_i     // level of single factor
	csCol.setting[0].levelsInGroup = 1

	// populate every row
	for row_i := 0; row_i < csm.rows; row_i++ {
		csm.addRow(csCol)

		// 1 or -1 depending on if the factor levels matched
		if level_i == levelMatrix[row_i][factor_i] {
			csCol.dataP[row_i] = ENTRY_A
		} else {
			csCol.dataP[row_i] = ENTRY_B
		}

		// add to sum of squares
		sum += csCol.dataP[row_i] * csCol.dataP[row_i]
	}

	// push into vector
	csm.data.push_back(csCol)
	*sumOfSquares = append(*sumOfSquares, sum)
}

func (csm *CSMatrix) print() {
	fmt.Println("CSMatrix:")
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		fmt.Printf("%s\t", csm.getColName(data.at(col_i)))
	}
	fmt.Println("")

	for row_i := 0; row_i < csm.rows; row_i++ {
		for col_i := 0; col_i < csm.getCols(); col_i++ {
			fmt.Printf("%v\t", csm.data[col_i].dataP[row_i])
		}
		fmt.Println("")
	}

	// verify that the column mappings are correct
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		if col_i != csm.getColIndex(csm.data[col_i]) {
			fmt.Println("Invalid Mapping!!")
		}
	}
}

func (csm *CSMatrix) countOccurrences(csCol *CSCol, occurrence *Occurrence, minSetting_i int, magnitude float64) {
	if occurrence.list == nil {
		return
	}

	for setting_i := minSetting_i; setting_i < csCol.factors; setting_i++ {
		occurrence.list[csCol.setting[setting_i].factor_i].count++
		occurrence.list[csCol.setting[setting_i].factor_i].magnitude += Abs(magnitude)

		csnm.countOccurrences(csCol, &occurrence.list[csCol.setting[setting_i].factor_i], setting_i+1, magnitude)
	}
}

func (csm *CSMatrix) getColIndex(csCol *CSCol) int {
	m := csm.mapping
	for setting_i := 0; setting_i < csCol.factors; setting_i++ {
		if m == nil {
			return -1
		}
		m = m.mapping[factorLevelMap[csCol.setting[setting_i].factor_i][csCol.setting[setting_i].index]]
	}
	return m.mappedTo
}

func (csm *CSMatrix) swapColumns(array []*CSCol, col_i1 int, col_i2 int) {
	tempCol := array[col_i1]
	array[col_i1] = array[col_i2]
	array[col_i2] = tempCol
}

func (csm *CSMatrix) swapRows(array []*CSCol, row_i1 int, row_i2 int) {
	var tempData float64
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		tempData = array[col_i].dataP[row_i1]
		array[col_i].dataP[row_i1] = array[col_i].dataP[row_i2]
		array[col_i].dataP[row_i2] = tempData
	}

	levelMatrix := csm.locatingArray.getLevelMatrix()
	tempLevel := levelMatrix[row_i1]
	levelMatrix[row_i1] = levelMatrix[row_i2]
	levelMatrix[row_i2] = tempLevel
}

// sort the array, given that some rows are already sorted
func (csm *CSMatrix) smartSort(array []*CSCol, sortedRows int) {
	var streakMin, streakMax int
	streakMin = 0

	for col_i := 1; col_i < csm.getCols(); col_i++ {
		// check if the streak ended
		if csm.compare(array[col_i-1], array[col_i], 0, sortedRows) < 0 {
			streakMax = col_i - 1
			if streakMin < streakMax {
				csm.rowSort(array, streakMin, streakMax, sortedRows, rows-sortedRows)
			}
			streakMin = col_i
		}
	}

	// add the final streak
	streakMax = csm.getCols() - 1
	if streakMin < streakMax {
		csm.rowSort(array, streakMin, streakMax, sortedRows, rows-sortedRows)
	}
}

func (csm *CSMatrix) rowSort(array []*CSCol, min int, max int, row_i int, row_len int) {
	if min >= max || row_len <= 0 {
		return
	}

	tempMin := min - 1
	tempMax := max + 1

	for {
		for tempMin < max && array[tempMin+1].dataP[row_i] == ENTRY_A {
			tempMin++
		}
		for tempMax > min && array[tempMax-1].dataP[row_i] == ENTRY_B {
			tempMax--
		}

		if tempMax-1 > tempMin+1 {
			csm.swapColumns(array, tempMin+1, tempMax-1)
		} else {
			break
		}
	}

	csm.rowSort(array, min, tempMin, row_i+1, row_len-1)
	csm.rowSort(array, tempMax, max, row_i+1, row_len-1)
}

func (csm *CSMatrix) sortByCoverable(array []*CSCol, min int, max int) int {
	if min >= max {
		return min
	}

	tempMin := min - 1
	tempMax := max + 1

	for {
		for tempMin < max && !array[tempMin+1].coverable {
			tempMin++
		}
		for tempMax > min && array[tempMax-1].coverable {
			tempMax--
		}

		if tempMax-1 > tempMin+1 {
			csm.swapColumns(array, tempMin+1, tempMax-1)
		} else {
			break
		}
	}

	// verification
	for col_i := min; col_i <= tempMin; col_i++ {
		if array[col_i].coverable {
			fmt.Println("mistake")
		}
	}

	return tempMin + 1
}

func (csm *CSMatrix) sortByTWayInteraction(array []*CSCol, min int, max int) int {
	if min >= max {
		return min
	}

	t := csm.locatingArray.getT()
	tempMin := min - 1
	tempMax := max + 1

	for {
		for tempMin < max && array[tempMin+1].factors < t {
			tempMin++
		}
		for tempMax > min && array[tempMax-1].factors == t {
			tempMax--
		}

		if tempMax-1 > tempMin+1 {
			csm.swapColumns(array, tempMin+1, tempMax-1)
		} else {
			break
		}
	}

	// verification
	for col_i := min; col_i <= tempMin; col_i++ {
		if array[col_i].factors == t {
			fmt.Println("mistake")
		}
	}

	return tempMin + 1
}

func (csm *CSMatrix) pathSort(array []*CSCol, path *Path, row_i int, nPaths *int, pathList *List) {
	if path.min == path.max {
		deletePath(path.entryA)
		deletePath(path.entryB)
		path.entryA = nil
		path.entryB = nil
		return
	} else if row_i >= rows {
		// add to list
		if pathList != nil {
			pathList.PushFront(path)
		}

		return
	}

	tempMin := path.min - 1
	tempMax := path.max + 1

	for {
		for tempMin < path.max && array[tempMin+1].dataP[row_i] == ENTRY_A {
			tempMin++
		}
		for tempMax > path.min && array[tempMax-1].dataP[row_i] == ENTRY_B {
			tempMax--
		}

		if tempMax-1 > tempMin+1 {
			csm.swapColumns(array, tempMin+1, tempMax-1)
		} else {
			break
		}
	}

	// verification
	for col_i := path.min; col_i <= tempMin; col_i++ {
		if array[col_i].dataP[row_i] != ENTRY_A {
			fmt.Println("mistake")
		}
	}
	for col_i = tempMax; col_i <= path.max; col_i++ {
		if array[col_i].dataP[row_i] != ENTRY_B {
			fmt.Println("mistake")
		}
	}
	if tempMin != tempMax-1 {
		fmt.Println("mistake")
	}

	if path.min <= tempMin {
		nPaths++
		// allocate memory if none exists
		if path.entryA == nill {
			path.entryA = &Path{}
			path.entryA.entryA = nill
			path.entryA.entryB = nill
		}

		// populate path for entryA
		path.entryA.min = path.min
		path.entryA.max = tempMin

		// sort path for entryA
		csm.pathSort(array, path.entryA, row_i+1, nPaths, pathList)
	} else {
		// delete unnecessary path for entryA
		csm.deletePath(path.entryA)
		path.entryA = nil
	}

	if tempMax <= path.max {
		nPaths++
		// allocate memory if none exists
		if path.entryB == nil {
			path.entryB = &Path{}
			path.entryB.entryA = nil
			path.entryB.entryB = nil
		}

		// populate path for entryB
		path.entryB.min = tempMax
		path.entryB.max = path.max

		// sort path for entryB
		csm.pathSort(array, path.entryB, row_i+1, nPaths, pathList)
	} else {
		csm.deletePath(path.entryB)
		path.entryB = nil
	}
}

func (csm *CSMatrix) deletePath(path *Path) {
	if path != nil {
		csm.deletePath(path.entryA)
		csm.deletePath(path.entryB)
		path = nil
	}
}

// Locating Array Checker
func (csm *CSMatrix) pathLAChecker(array []*CSCol, pathA *Path, pathB *Path, row_i int, k int,
	score *int64, settingToResample **FactorSetting, rowContributions []int64) {
	if k == 0 || pathA == nil || pathB == nil || pathA.min == pathB.max {
		return
	} else if row_i == rows {
		//		cout << "Issue " << getColName(array[pathA->min]) << " vs " << getColName(array[pathB->max]) << endl;
		if pathA == pathB {
			score += int64(k) * (int64(pathA.max-pathA.min+1) * int64(pathA.max-pathA.min)) / 2
		} else {
			score += int64(k) * int64(pathA.max-pathA.min+1) * int64(pathB.max-pathB.min+1)
		}

		// set a setting to resample
		if settingToResample == NULL {
			columnToResample := -1

			/*
				int offset;
				do {
					offset = rand() % (pathA->max - pathA->min + 1 + pathB->max - pathB->min + 1);
					if (offset <= pathA->max - pathA->min) {
						columnToResample = pathA->min + offset;
					} else {
						offset -= (pathA->max - pathA->min + 1);
						columnToResample = pathB->min + offset;
					}
				} while (array[columnToResample]->factors <= 0);
			*/

			// find a pair of columns that is distinguishable
			for i_a := pathA.min; i_a <= pathA.max && columnToResample == -1; i_a++ {
				var int i_b
				if pathA == pathB {
					i_b = i_a + 1
				}
				for i_b; i_b <= pathB.max && columnToResample == -1; i_b++ {
					if csm.checkDistinguishable(array[i_a], array[i_b]) {
						if rand.Intn(2) && array[i_a].factors > 0 {
							columnToResample = i_a
						} else if array[i_b].factors > 0 {
							columnToResample = i_b
						}
					} else {
						//						score--;
						//						cout << "Not distinguishable: " << getColName(array[i_a]) << " vs " << getColName(array[i_b]) << endl;
					}
				}
			}

			// randomly choose a setting in the column to resample
			if columnToResample != -1 {
				settingToResample = &array[columnToResample].setting[rand.Intn(array[columnToResample].factors)]
			}
		}

		return
	}

	var pathAentryA, pathAentryB, pathBentryA, pathBentryB *Path

	if pathA.min == pathA.max {
		if array[pathA.min].dataP[row_i] == ENTRY_A {
			pathAentryA = pathA
		}
		if array[pathA.min].dataP[row_i] == ENTRY_B {
			pathAentryB = pathA
		}
	} else {
		pathAentryA = pathA.entryA
		pathAentryB = pathA.entryB
	}

	if pathB.min == pathB.max {
		if array[pathB.min].dataP[row_i] == ENTRY_A {
			pathBentryA = pathB
		}
		if array[pathB.min].dataP[row_i] == ENTRY_B {
			pathBentryB = pathB
		}
	} else {
		pathBentryA = pathB.entryA
		pathBentryB = pathB.entryB
	}

	csm.pathLAChecker(array, pathAentryA, pathBentryA, row_i+1, k, score, settingToResample, rowContributions)
	csm.pathLAChecker(array, pathAentryB, pathBentryB, row_i+1, k, score, settingToResample, rowContributions)
	csm.pathLAChecker(array, pathAentryA, pathBentryB, row_i+1, k-1, score, settingToResample, rowContributions)

	// add row contributions
	if rowContributions != nil && pathAentryA != nil && pathBentryB != nil {
		rowContributions[row_i] += (pathAentryA.max - pathAentryA.min + 1) * (pathBentryB.max - pathBentryB.min + 1)
	}

	if pathA != pathB {
		csm.pathLAChecker(array, pathAentryB, pathBentryA, row_i+1, k-1, score, settingToResample, rowContributions)

		// add row contributions
		if rowContributions != nil && pathAentryB != nil && pathBentryA != nil {
			rowContributions[row_i] += (pathAentryB.max - pathAentryB.min + 1) * (pathBentryA.max - pathBentryA.min + 1)
		}

	}
}

// Detecting Array Checker
func (csm *CSMatrix) pathDAChecker(array []*CSCol, pathA *Path, pathB *Path, row_i int, k int,
	score *int64, settingToResample **FactorSetting, rowContributions []int64) {
	if k == 0 || pathA == nil || pathB == nil || pathA.min == pathB.max {
		return
	} else if row_i == rows {
		//		cout << getColName(array[pathA->min]) << " vs " << getColName(array[pathB->max]) << endl;
		if pathA == pathB {
			score += int64(k) * (int64(pathA.max-pathA.min+1) * int64(pathA.max-pathA.min)) / 2
		} else {
			score += int64(k) * int64(pathA.max-pathA.min+1) * int64(pathB.max-pathB.min+1)
		}

		// set a setting to resample
		if settingToResample == NULL {
			var columnToResample int
			var offset int
			for ok := true; ok; ok = array[columnToResample].factors <= 0 {
				offset = rand.Intn(pathA.max - pathA.min + 1 + pathB.max - pathB.min + 1)
				if offset <= pathA.max-pathA.min {
					columnToResample = pathA.min + offset
				} else {
					offset -= (pathA.max - pathA.min + 1)
					columnToResample = pathB.min + offset
				}
			}

			settingToResample = &array[columnToResample].setting[rand.Intn(array[columnToResample].factors)]
		}

		return
	}

	var pathAentryA, pathAentryB, pathBentryA, pathBentryB *Path

	if pathA.min == pathA.max {
		if array[pathA.min].dataP[row_i] == ENTRY_A {
			pathAentryA = pathA
		}
		if array[pathA.min].dataP[row_i] == ENTRY_B {
			pathAentryB = pathA
		}
	} else {
		pathAentryA = pathA.entryA
		pathAentryB = pathA.entryB
	}

	if pathB.min == pathB.max {
		if array[pathB.min].dataP[row_i] == ENTRY_A {
			pathBentryA = pathB
		}
		if array[pathB.min].dataP[row_i] == ENTRY_B {
			pathBentryB = pathB
		}
	} else {
		pathBentryA = pathB.entryA
		pathBentryB = pathB.entryB
	}

	csm.pathDAChecker(array, pathAentryA, pathBentryA, row_i+1, k, score, settingToResample, rowContributions)
	csm.pathDAChecker(array, pathAentryB, pathBentryB, row_i+1, k, score, settingToResample, rowContributions)
	csm.pathDAChecker(array, pathAentryB, pathBentryA, row_i+1, k, score, settingToResample, rowContributions)
	csm.pathDAChecker(array, pathAentryA, pathBentryB, row_i+1, k-1, score, settingToResample, rowContributions)

	// add row contributions
	if rowContributions != nil && pathAentryA != nil && pathBentryB != nil {
		rowContributions[row_i] += (pathAentryA.max - pathAentryA.min + 1) * (pathBentryB.max - pathBentryB.min + 1)
	}
}

func (csm *CSMatrix) quickSort(array []*CSCol, min int, max int, row_top int, row_len int) {
	if min == max {
		return
	}

	tempMin := min
	tempMax := max
	midCol := array[(min+max)/2]

	for {
		for tempMin <= tempMax && tempMin < max && csm.compare(array[tempMin], midCol, row_top, row_len) < 0 {
			tempMin++
		}

		for tempMin <= tempMax && tempMax > min && csm.compare(array[tempMax], midCol, row_top, row_len) > 0 {
			tempMax--
		}

		if tempMin == tempMax {
			if tempMin == min {
				tempMin++
			} else {
				tempMax--
			}
			break
		} else if tempMin <= tempMax {
			csm.swapColumns(array, tempMin, tempMax)
			tempMin++
			tempMax--
		} else {
			break
		}
	}

	csm.quickSort(array, min, tempMax, row_top, row_len)
	csm.quickSort(array, tempMin, max, row_top, row_len)
}

func (csm *CSMatrix) compare(csCol1 *CSCol, csCol2 *CSCol, row_top int, row_len int) int {
	for i := row_top; i < row_top+row_len; i++ {
		if csCol1.dataP[row_top] == csCol2.dataP[row_top] {
			continue
		} else if csCol1.dataP[row_top] < csCol2.dataP[row_top] {
			return -1
		} else if csCol1.dataP[row_top] > csCol2.dataP[row_top] {
			return 1
		}
	}
}

// LEGACY (to be removed later)
//cout << "Checking advanced for " << k << " differences" << endl;
//score = checkAdvanced(array, k, 0, getCols() - 1, 0, rows, settingToResample);
//cout << "Advanced Elapsed: " << elapsedTime << endl;
func (csm *CSMatrix) checkAdvanced(array []*CSCol, k int, min int, max int, row_top int, row_len int,
	settingToResample **FactorSetting) int64 {
	// no more differences to find or only 1 column
	if k <= 0 || min >= max {
		return 0
	} else if k > row_len {

		//		cout << "Issue " << getColName(array[min]) << " vs " << getColName(array[max]) << endl;
		//		cout << "Looking for " << k << " differences, but there are " << row_len << " rows left " << min << " vs " << max << endl;

		csCol1 := array[min]
		csCol2 := array[max]

		if settingToResample == nil {
			// choose random factor+level for resampling
			var csColToResample *CSCol

			if csCol1.factors != 0 && csCol2.factors != 0 {
				// neither has 0 factors, choose randomly
				if rand.Intn(2) {
					csColToResample = csCol1
				} else {
					csColToResample = csCol2
				}
			} else if csCol1.factors != 0 {
				// csCol1 has factors so choose it (csCol2 has no factors)
				csColToResample = csCol1
			} else if csCol2.factors != 0 {
				// csCol2 has factors so choose it (csCol1 has no factors)
				csColToResample = csCol2
			}

			// chose a random factor setting to resample
			if csColToResample != nil {
				settingToResample = &csColToResample.setting[rand.Intn(csColToResample.factors)]
			}
		}

		cols := max - min + 1
		return int64(math.Pow((cols * cols), k))
	} else {
		var score int64

		csm.rowSort(array, min, max, row_top, 1)

		var divide int
		for divide = min; divide < max; divide++ {
			if csm.compare(array[divide], array[divide+1], row_top, 1) != 0 {
				break
			}
		}
		for col_i := divide + 1; col_i < max; col_i++ {
			if csm.compare(array[col_i], array[col_i+1], row_top, 1) != 0 {
				fmt.Println("Integrity Issue")
			}
		}

		// check the left side
		score += checkAdvanced(array, k, min, divide, row_top+1, row_len-1, settingToResample)
		// check the right side
		score += checkAdvanced(array, k, divide+1, max, row_top+1, row_len-1, settingToResample)
		// check both
		score += checkAdvanced(array, k-1, min, max, row_top+1, row_len-1, settingToResample)

		return score
	}
}

// FUNCTIONS FOR 'fixla'

func (csm *CSMatrix) addRow(array []*CSCol, levelRow []byte) {
	csm.rows++

	csm.locatingArray.addLevelRow(levelRow)

	// add a row to each column of the CS matrix and populate
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		csm.addArray(array[col_i])
		populateColumnData(array[col_i], csm.locatingArray.getLevelMatrix(), rows-1, 1)
	}
}

func (csm *CSMatrix) remRow(array []*CSCol) {
	csm.rows--

	levelRow := csm.locatingArray.remLevelRow()
	// delete levelRow;

	// remove a row from each column of the CS matrix
	for col_i := 0; col_i < csm.getCols(); col_i++ {
		csm.remArray(array[col_i])
	}
}

func (csm *CSMatrix) resizeArray(array []*CSCol, newRows int) {
	factors := csm.locatingArray.getFactors()
	nConGroups := csm.locatingArray.getNConGroups()
	conGroups := csm.locatingArray.getConGroups()

	// increase size as needed
	for newRows > csm.rows {
		// allocate memory for new row of locating array
		levelRow := make([]byte, factors)

		// generate random row
		for factor_i := 0; factor_i < factors; factor_i++ {
			levelRow[factor_i] = rand.Intn(groupingInfo[factor_i].levels)
		}

		// resample constraint groups
		for conGroup_i := 0; conGroup_i < nConGroups; conGroup_i++ {
			conGroups[conGroup_i].randPopulateLevelRow(levelRow)
		}

		csm.addRow(array, levelRow)
	}

	// decrease size as needed
	for newRows < csm.rows && csm.rows > 1 {
		csm.remRow(array)
	}
}

func (csm *CSMatrix) checkColumnCoverability(csCol *CSCol) bool {
	factors := csm.locatingArray.getFactors()
	nConGroups := csm.locatingArray.getNConGroups()
	conGroups := csm.locatingArray.getConGroups()

	requireLevelRow := make([]byte, factors)
	avoidLevelRow := make([]byte, factors)
	for factor_i := 0; factor_i < factors; factor_i++ {
		requireLevelRow[factor_i] = -1
		avoidLevelRow[factor_i] = -1
	}

	// set to require csCol1 settings
	for setting_i := 0; setting_i < csCol.factors; setting_i++ {
		requireLevelRow[csCol.setting[setting_i].factor_i] = csCol.setting[setting_i].index
	}

	// check if possible for every constraint group
	satisfiableInGroups := true
	for conGroup_i := 0; conGroup_i < nConGroups; conGroup_i++ {
		if !conGroups[conGroup_i].satisfiableInGroupLA(requireLevelRow, avoidLevelRow) {
			satisfiableInGroups = false
			break
		}
	}

	/* for factors not in a constraint group, simply set their settings in csCol1 to be ENTRY_A,
	and the check passes because the interaction cannot be only */

	return satisfiableInGroups
}

func (csm *CSMatrix) checkDistinguishable(csCol1, csCol2 *CSCol) bool {
	return csm.checkOneWayDistinguishable(csCol1, csCol2) || csm.checkOneWayDistinguishable(csCol2, csCol1)
}

func (csm *CSMatrix) checkOneWayDistinguishable(csCol1, csCol2 *CSCol) bool {
	factors := csm.locatingArray.getFactors()
	nConGroups = csm.locatingArray.getNConGroups()
	conGroups := csm.locatingArray.getConGroups()

	requireLevelRow := make([]byte, factors)
	avoidLevelRow := make([]byte, factors)
	for factor_i := 0; factor_i < factors; factor_i++ {
		requireLevelRow[factor_i] = -1
		avoidLevelRow[factor_i] = -1
	}

	// set to require csCol1 settings
	for setting_i := 0; setting_i < csCol1.factors; setting_i++ {
		requireLevelRow[csCol1.setting[setting_i].factor_i] = csCol1.setting[setting_i].index
	}

	// at least one satisfying setting in csCol2 must be avoidable for every constraint group
	settingExists := false

	// try to avoid at least one csCol2 setting
	for setting_i = 0; setting_i < csCol2.factors; setting_i++ {
		// we are trying to avoid "csCol2->setting[setting_i].factor_i != csCol2->setting[setting_i].index" in this iteration
		// we cannot require the setting we are trying to avoid
		if requireLevelRow[csCol2.setting[setting_i].factor_i] != csCol.setting[setting_i].index {
			avoidLevelRow[csCol2.setting[setting_i].factor_i] = csCol2.setting[setting_i].index

			// check if possible for every constraint group
			satisfiableInGroups := true
			for conGroup_i := 0; conGroup_i < nConGroups; conGroup_i++ {
				if !conGroups[conGroup_i].satisfiableInGroupLA(requireLevelRow, avoidLevelRow) {
					satisfiableInGroups = false
					break
				}
			}

			if satisfiableInGroups {
				settingExists = true
				break
			}

			avoidLevelRow[csCol2.setting[setting_i].factor_i] = -1
		}
	}

	/* for factors not in a constraint group, simply set their settings in csCol1 to be ENTRY_A,
	and the check passes because the interaction cannot be only */

	return settingExists
}

func (csm *CSMatrix) randomizeArray(array []*CSCol) {
	factors := csm.locatingArray.getFactors()
	cols := csm.getCols()

	// resample entire locating array
	levelMatrix := csm.locatingArray.getLevelMatrix()
	for row_i := 0; row_i < csm.rows; row_i++ {
		for factor_i := 0; factor_i < factors; factor_i++ {
			levelMatrix[row_i][factor_i] = rand.Intn(csm.groupingInfo[factor_i].levels)
		}
	}

	// populate each column of the CS matrix
	for col_i := 0; col_i < cols; col_i++ {
		csm.populateColumnData(array[col_i], csm.locatingArray.getLevelMatrix(), 0, rows)
	}
}

func (csm *CSMatrix) addRowFix(array []*CSCol, csScore *int64) {
	// the maximum seconds without finding a better score
	secsMax := 2

	// total factors in locating array
	factors := csm.locatingArray.getFactors()

	var bestScore int64
	var newScore int64
	var csCol, bestCol *CSCol

	cols := csm.getCols()

	// for backing up the order of array
	backupArray := make([]*CSCol, cols)

	// allocate memory for new row of locating array
	levelRow := make([]byte, factors)
	oldLevelRow := make([]byte, factors)
	newLevelRow := make([]byte, factors)
	bestLevelRow := make([]byte, factors)

	// track if factors of new row are finalized
	finalized = make([]bool, factors)

	// generate random row non-finalized
	for factor_i := 0; factor_i < factors; factor_i++ {
		finalized[factor_i] = false
		levelRow[factor_i] = rand.Intn(csm.groupingInfo[factor_i].levels)
	}

	// add the row to locating array
	csm.addRow(array, levelRow)
	fmt.Println("The matrix now has ", rows, " rows")

	// smartly sort the array and score it
	csm.smartSort(array, csm.rows-1)
	*csScore = csm.getArrayScore(array)

	//	cout << "Score after random non-finalized row: " << csScore << endl;

	// used to track duplicate columns in CS matrix
	duplicate := false
	lastPairMatched := false

	// continue adding finalizing factors while CS score decreases (improves)
	for {

		// find the best column index (with the best score)
		bestScore = &csScore // the original best is the original score
		bestCol = nil        // no best column yet

		// backup the sorted order of the array
		copy(backupArray, array)

		// grab initial time
		start := time.Now()

		// go through all columns of CS matrix
		for col_i = 0; col_i < csm.getCols()-1; col_i++ {
			elapsedTime := time.Since(start).Seconds()
			if elapsedTime > secsMax {
				break
			}

			// check for duplicate column
			if csm.compare(array[col_i], array[col_i+1], 0, csm.rows) >= 0 {
				duplicate = true
				lastPairMatched = true
			} else if lastPairMatched {
				duplicate = true
				lastPairMatched = false
			} else {
				duplicate = false
			}

			// retrieve CSCol (and duplicate) from sorted array
			csCol = array[col_i]

			// if duplicate and the column has factors to change (column is not the INTERCEPT)
			if duplicate && csCol.dataP[rows-1] != 1 && csCol.factors > 0 {

				// grab duplicate column
				csDup := array[col_i+1]

				// the goal is to make the last row of csCol be a 1
				// all factors must be changed, check if this change conflicts with finalized factors
				changeAllowed := true
				var factor_i int

				// go through all relevant factors for this column
				for setting_i = 0; setting_i < csCol.factors; setting_i++ {
					factor_i = csCol.setting[setting_i].factor_i

					// backup the factor level
					oldLevelRow[factor_i] = levelRow[factor_i]

					// check if finalized yet
					if !finalized[factor_i] {
						// update factor to match the current column
						newLevelRow[factor_i] = csCol.setting[setting_i].index +
							rand.Intn(csCol.setting[setting_i].levelsInGroup)
					} else {
						// ensure this will actually make this column a 1
						changeAllowed &= (newLevelRow[factor_i] >= csCol.setting[setting_i].index &&
							newLevelRow[factor_i] < csCol.setting[setting_i].index+csCol.setting[setting_i].levelsInGroup)
					}

				}

				// check if this column change is allowed
				if changeAllowed {
					// CAREFUL!!!
					// col_i cannot be used after smartSort and before the rollback
					// because the array is reordered. Use csCol

					// try making the changes and see if the CS score gets smaller
					newScore = &csScore
					for setting_i := 0; setting_i < csCol.factors; setting_i++ {
						factor_i = csCol.setting[setting_i].factor_i

						if levelRow[factor_i] != newLevelRow[factor_i] {
							// update columns
							levelRow[factor_i] = newLevelRow[factor_i]
							csm.repopulateColumns(factor_i, oldLevelRow[factor_i], rows-1, 1)
							csm.repopulateColumns(factor_i, newLevelRow[factor_i], rows-1, 1)
						}
					}

					// smart sort and get new score
					csm.smartSort(array, rows-1)
					newScore = csm.getArrayScore(array)

					// check if we have a better score
					if newScore < bestScore {
						// store the better score as the best, and save the column
						bestScore = newScore
						bestCol = csCol

						for factor_i := 0; factor_i < factors; factor_i++ {
							bestLevelRow[factor_i] = newLevelRow[factor_i]
						}
					}

					// rollback the change
					for setting_i := 0; setting_i < csCol.factors; setting_i++ {
						factor_i = csCol.setting[setting_i].factor_i

						if oldLevelRow[factor_i] != levelRow[factor_i] {
							// rollback columns and get new score
							levelRow[factor_i] = oldLevelRow[factor_i]
							csm.repopulateColumns(factor_i, oldLevelRow[factor_i], rows-1, 1)
							csm.repopulateColumns(factor_i, newLevelRow[factor_i], rows-1, 1)
						}
					}

					// restore old order
					copy(array, backupArray)

				}
			}
		}

		// check if we found a better score
		if bestScore < csScore {
			//			cout << "Found better score! " << bestScore << endl;

			// make the change to improve the score
			var factor_i int
			for setting_i := 0; setting_i < bestCol.factors; setting_i++ {
				factor_i = bestCol.setting[setting_i].factor_i

				if bestLevelRow[factor_i] != levelRow[factor_i] {
					// update columns and get new score
					oldLevelRow[factor_i] = levelRow[factor_i]
					levelRow[factor_i] = bestLevelRow[factor_i]
					csm.repopulateColumns(factor_i, oldLevelRow[factor_i], rows-1, 1)
					csm.repopulateColumns(factor_i, bestLevelRow[factor_i], rows-1, 1)
				}

				// finalize the changes
				finalized[factor_i] = true
			}

			// do a final smart sort and finalize the changed factors
			csm.smartSort(array, rows-1)
			csScore = csm.getArrayScore(array)

		} else {
			break
		}
	}

	fmt.Println("Score after finalized row: ", &csScore)
}

func (csm *CSMatrix) getArrayScore(array []*CSCol) int64 {
	var (
		streak     int64 = 0
		squaredSum int64 = 0
	)

	for col_i := 0; col_i < csm.getCols()-1; col_i++ {

		// increment streak
		streak++

		// check if the streak ended
		if csm.compare(array[col_i], array[col_i+1], 0, csm.rows) < 0 {
			squaredSum += streak * streak
			streak = 0
		} else if csm.compare(array[col_i], array[col_i+1], 0, csm.rows) > 0 {
			fmt.Println("Mistake in CS matrix at column: ", col_i)
		} else {
			if csm.checkDistinguishable(array[col_i], array[col_i+1]) {
				// cout << "Duplicates found: " << getColName(array[col_i]) << " vs " << getColName(array[col_i + 1]) << endl;
			}
		}
	}

	// add the final streak
	streak++
	squaredSum += streak * streak

	return (squaredSum - getCols())
}

func (csm *CSMatrix) getBruteForceArrayScore(array []*CSCol, k int) int64 {
	var (
		score             int64 = 0
		differences       int
		indistinguishable int = 0
	)

	for col_i1 := 0; col_i1 < csm.getCols()-1; col_i1++ {
		for col_i2 := col_i1 + 1; col_i2 < csm.getCols(); col_i2++ {
			if array[col_i1].coverable && array[col_i2].coverable {
				if !csm.checkDistinguishable(array[col_i1], array[col_i2]) {
					fmt.Println("Indistinguishable pair: ", csm.getColName(array[col_i1]), " vs ", csm.getColName(array[col_i2]))
					indistinguishable++
				}

				differences = 0

				for row_i := 0; row_i < csm.rows; row_i++ {
					if array[col_i1].dataP[row_i] != array[col_i2].dataP[row_i] {
						differences++
					}
				}

				if differences < k {
					score += (k - differences)
				}
			}
		}
	}

	fmt.Println("Indistinguishable pairs: ", indistinguishable)
	fmt.Println("Minimum score: ", (indistinguishable * k))

	return score
}

func (csm *CSMatrix) writeResponse(responseDir string, responseCol string, terms int, coefficients []float64, columns []int) {
	fmt.Println(responseDir)

	// initialize all responses to 0
	responses := make([]float64, csm.rows)
	for row_i := 0; row_i < csm.rows; row_i++ {
		responses[row_i] = 0
	}

	// add contributions from all terms
	for term_i := 0; term_i < terms; term_i++ {

		csCol := csm.data[columns[term_i]]
		coefficient := coefficients[term_i]

		fmt.Println(coefficient, " * ", csm.getColName(csCol))

		for row_i := 0; row_i < csm.rows; row_i++ {
			responses[row_i] += coefficient * csCol.dataP[row_i]
		}
	}

	// open a response file for writing
	file, err := os.Create(responseDir + "/Response.tsv")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	w := bufio.NewWriter(file)

	// write response headers
	fmt.Fprintln(w, csm.rows)
	fmt.Fprintln(w, responseCol)

	// write response values
	for row_i := 0; row_i < csm.rows; row_i++ {
		fmt.Fprintln(w, responses[row_i])
	}

	w.Flush()
}
