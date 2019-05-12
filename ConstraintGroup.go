package main

import (
	"bufio"
	"fmt"
	"log"
	"math/rand"
	"strconv"
	"strings"
)

func newResult(array *LocatingArray, words []string, index *int) *Result {
	operand := words[*index]
	*index = *index + 1

	r := &Result{
		array:   array,
		operand: operand,
	}

	if operand == "F" {
		factor_i, err := strconv.Atoi(words[*index])
		if err != nil {
			log.Fatal(err)
		}
		*index = *index + 1
		r.factor_i = factor_i
	} else if operand == "C" {
		value, err := strconv.ParseFloat(words[*index], 32)
		if err != nil {
			log.Fatal(err)
		}
		*index = *index + 1
		r.value = value
	} else if operand == "==" || operand == "<=" || operand == ">" || operand == "IF" ||
		operand == "+" || operand == "*" || operand == "/" {
		r.result1 = newResult(array, words, index)
		r.result2 = newResult(array, words, index)
	} else {
		log.Fatal("Invalid operand")
	}

	return r
}

type Result struct {
	array    *LocatingArray
	result1  *Result // LHS of Eq equation
	result2  *Result // RHS of Eq equation
	operand  string  //
	value    float64 //
	factor_i int     //
}

func (op *Result) getResult(test int) float64 {
	if op.operand == "==" {
		if op.result1.getResult(test) == op.result2.getResult(test) {
			return 1
		} else {
			return 0
		}
	} else if op.operand == "<=" {
		if op.result1.getResult(test) <= op.result2.getResult(test) {
			return 1
		} else {
			return 0
		}
	} else if op.operand == ">" {
		if op.result1.getResult(test) > op.result2.getResult(test) {
			return 1
		} else {
			return 0
		}
	} else if op.operand == "IF" {
		if !(op.result1.getResult(test) == 1 || op.result2.getResult(test) == 1) {
			return 1
		} else {
			return 0
		}
	} else if op.operand == "+" {
		return op.result1.getResult(test) + op.result2.getResult(test)
	} else if op.operand == "*" {
		return op.result1.getResult(test) * op.result2.getResult(test)
	} else if op.operand == "/" {
		return op.result1.getResult(test) / op.result2.getResult(test)
	} else if op.operand == "C" {
		return op.value
	} else if op.operand == "F" {
		return op.array.getFactorData().getNumericFactorLevel(op.factor_i, op.array.getLevelMatrix()[test][op.factor_i])
	} else {
		log.Fatal("Invalid operand")
	}

	return 0
}

func (op *Result) string() string {
	if op.operand == "==" {
		return fmt.Sprintf("\t==%s%s", op.result1.string(), op.result2.string())
	} else if op.operand == "<=" {
		return fmt.Sprintf("\t<=%s%s", op.result1.string(), op.result2.string())
	} else if op.operand == ">" {
		return fmt.Sprintf("\t>%s%s", op.result1.string(), op.result2.string())
	} else if op.operand == "IF" {
		return fmt.Sprintf("\tIF%s%s", op.result1.string(), op.result2.string())
	} else if op.operand == "+" {
		return fmt.Sprintf("\t+%s%s", op.result1.string(), op.result2.string())
	} else if op.operand == "*" {
		return fmt.Sprintf("\t*%s%s", op.result1.string(), op.result2.string())
	} else if op.operand == "/" {
		return fmt.Sprintf("\t/%s%s", op.result1.string(), op.result2.string())
	} else if op.operand == "C" {
		return fmt.Sprintf("\tC\t%v", op.value)
	} else if op.operand == "F" {
		return fmt.Sprintf("\tF\t%v", op.factor_i)
	} else {
		log.Fatal("Invalid operand")
	}

	return ""
}

type ConstraintGroup struct {
	groupLA         *LocatingArray
	constraints     int
	boolConstraints []*Result

	// there exists a weight window for each row (larger weight window means more likely to be chosen)
	weightMin     []int // bottom of weight window for each row
	weightMax     []int // top of weight window for each row
	weightRandMax int
	factors       int
	factorIndeces []int
}

func newConstraintGroup(array *LocatingArray, lines []string) *ConstraintGroup {
	cg := &ConstraintGroup{}

	// load constraint group factors
	words := strings.Fields(lines[0])
	factors, err := strconv.Atoi(words[0])
	if err != nil {
		log.Fatal(err)
	}
	cg.factors = factors

	factorIndeces := make([]int, factors)
	for factor_i := 0; factor_i < factors; factor_i++ {
		val, err := strconv.Atoi(words[factor_i+1])
		if err != nil {
			log.Fatal(err)
		}
		factorIndeces[factor_i] = val
	}
	cg.factorIndeces = factorIndeces

	// allocate the level counts for the factors
	levelCounts := make([]int, factors)
	for factor_i := 0; factor_i < factors; factor_i++ {
		levelCounts[factor_i] = array.getGroupingInfo()[factorIndeces[factor_i]].levels
	}

	// create constraint group LA
	cg.groupLA = NewLocatingArrayFromFactors(factors, levelCounts)

	words = strings.Fields(lines[1])
	constraints, err := strconv.Atoi(words[0])
	if err != nil {
		log.Fatal(err)
	}
	cg.constraints = constraints

	// allocate memory for BoolResult objects
	boolConstraints := make([]*Result, constraints) // Always 1?

	for constraint_i := 0; constraint_i < constraints; constraint_i++ {
		index := 0
		boolConstraints[constraint_i] = newResult(array, words[1:], &index)
	}

	cg.boolConstraints = boolConstraints

	// link to factors in array
	for factor_i := 0; factor_i < factors; factor_i++ {
		array.getGroupingInfo()[factorIndeces[factor_i]].conGroup = cg
		array.getGroupingInfo()[factorIndeces[factor_i]].conGroupIndex = factor_i
	}

	// populate groupLA
	levelRow := make([]int, array.getFactors())
	factorLevels := make([]int, factors)

	// loop through all locating array factors and set all to 0
	for factor_i := 0; factor_i < array.getFactors(); factor_i++ {
		levelRow[factor_i] = 0
	}

	// add row of 0s to the array for use in recursive function
	array.addLevelRow(levelRow)

	// count of each level occurrence (used for weighting)
	settingCount := make([][]int, factors)
	for factor_i := 0; factor_i < factors; factor_i++ {
		settingCount[factor_i] = make([]int, levelCounts[factor_i])
		for level_i := 0; level_i < levelCounts[factor_i]; level_i++ {
			settingCount[factor_i][level_i] = 0
		}
	}

	// count of full-factorial combinations (used for weighting)
	fullFactorialCount := 0

	// call recursive populator
	cg.populateGroupLA(array, levelRow, factorLevels, 0, settingCount, &fullFactorialCount)

	// assign row weights based on groupLA entries
	groupLevelMatrix := cg.groupLA.getLevelMatrix()
	rowWeights := make([]float64, cg.groupLA.getTests())
	weightMin := make([]int, cg.groupLA.getTests())
	weightMax := make([]int, cg.groupLA.getTests())
	var prevWeight float64 = 0

	for row_i := 0; row_i < cg.groupLA.getTests(); row_i++ {
		rowWeights[row_i] = prevWeight
		for col_i := 0; col_i < factors; col_i++ {
			rowWeights[row_i] += float64(fullFactorialCount / settingCount[col_i][groupLevelMatrix[row_i][col_i]])
			//			cout << (int)groupLevelMatrix[row_i][col_i] << "\t";
		}

		weightMin[row_i] = int(prevWeight)
		weightMax[row_i] = int(rowWeights[row_i] - 1)

		//		cout << rowWeights[row_i] << "\t" << rowWeights[row_i] - prevWeight << "\t" << weightMin[row_i] << "\t" << weightMax[row_i] << endl;

		prevWeight = rowWeights[row_i]
	}

	cg.weightMin = weightMin
	cg.weightMax = weightMax
	cg.weightRandMax = int(prevWeight)

	levelRow = array.remLevelRow()

	return cg
}

func (cg *ConstraintGroup) populateGroupLA(array *LocatingArray, levelRow []int,
	factorLevels []int, fixedFactors int, settingCount [][]int, fullFactorialCount *int) {
	if fixedFactors == cg.factors {
		//		for (int factor_i = 0; factor_i < factors; factor_i++) {
		//			cout << (int)factorLevels[factor_i] << "\t";
		//		}
		//		cout << endl;

		// increment full-factorial count
		*fullFactorialCount++

		// loop through all factors in constraint group
		for factor_i := 0; factor_i < cg.factors; factor_i++ {
			levelRow[cg.factorIndeces[factor_i]] = factorLevels[factor_i]
		}

		// add factorLevels to groupLA if group constraints are satisfied
		if cg.getResult(0) {
			groupLevelRow := make([]int, cg.factors)
			for factor_i := 0; factor_i < cg.factors; factor_i++ {
				groupLevelRow[factor_i] = factorLevels[factor_i]

				// increment the setting count
				settingCount[factor_i][factorLevels[factor_i]]++
			}
			cg.groupLA.addLevelRow(groupLevelRow)
		}
	} else {
		// call recursively for full-factorial design
		for level_i := 0; level_i < array.getGroupingInfo()[cg.factorIndeces[fixedFactors]].levels; level_i++ {
			factorLevels[fixedFactors] = level_i
			cg.populateGroupLA(array, levelRow, factorLevels, fixedFactors+1, settingCount, fullFactorialCount)
		}
	}
}

func (cg *ConstraintGroup) satisfiableInGroupLA(requireLevelRow []int, avoidLevelRow []int) bool {
	levelMatrix := cg.groupLA.getLevelMatrix()

	for test_i := 0; test_i < cg.groupLA.getTests(); test_i++ {
		match := true

		// check requireLevelRow
		for factor_i := 0; factor_i < cg.groupLA.getFactors() && match; factor_i++ {
			// check if the required level disagrees with the LA row
			if requireLevelRow[cg.factorIndeces[factor_i]] != -1 &&
				requireLevelRow[cg.factorIndeces[factor_i]] != levelMatrix[test_i][factor_i] {
				match = false
			}
		}

		// check avoidLevelRow
		for factor_i := 0; factor_i < cg.groupLA.getFactors() && match; factor_i++ {
			// check if the avoided level disagrees with the LA row
			if avoidLevelRow[cg.factorIndeces[factor_i]] != -1 &&
				avoidLevelRow[cg.factorIndeces[factor_i]] == levelMatrix[test_i][factor_i] {
				match = false
			}
		}

		// check if the row satisfies
		if match {
			return true
		}
	}

	return false
}

func (cg *ConstraintGroup) getResult(test int) bool {
	for constraint_i := 0; constraint_i < cg.constraints; constraint_i++ {
		// TODO float equality check careful
		if cg.boolConstraints[constraint_i].getResult(test) == 0 {
			return false
		}
	}
	return true
}

func (cg *ConstraintGroup) randPopulateLevelRow(levelRow []int) {
	// groupingInfo := cg.groupLA.getGroupingInfo()
	levelMatrix := cg.groupLA.getLevelMatrix()

	weightRand := rand.Intn(cg.weightRandMax)

	// use binary search to find weight window with weightRand
	botRow := 0
	topRow := cg.groupLA.getTests()

	for topRow != botRow {
		row_i := (topRow + botRow) / 2

		if weightRand <= cg.weightMax[row_i] && weightRand >= cg.weightMin[row_i] {
			topRow = row_i
			botRow = row_i
		} else if weightRand < cg.weightMin[row_i] {
			topRow = row_i - 1
		} else if weightRand > cg.weightMax[row_i] {
			botRow = row_i + 1
		} else {
			// should never reach here
			log.Fatal("Mistake in binary search")
		}
	}

	for factor_i := 0; factor_i < cg.factors; factor_i++ {
		levelRow[cg.factorIndeces[factor_i]] = levelMatrix[topRow][factor_i]
	}
}

func (cg *ConstraintGroup) writeToStream(w *bufio.Writer) {
	// write constraint group factors
	factorsAndIndeces := make([]string, cg.factors+1)
	factorsAndIndeces[0] = string(cg.factors)
	for factor_i := 0; factor_i < cg.factors; factor_i++ {
		factorsAndIndeces[factor_i+1] = string(cg.factorIndeces[factor_i])
	}
	fmt.Fprintln(w, strings.Join(factorsAndIndeces, "\t"))

	// write constraints
	constraints := make([]string, cg.constraints+1)
	constraints[0] = strconv.Itoa(cg.constraints)
	for constraint_i := 0; constraint_i < cg.constraints; constraint_i++ {
		constraints[constraint_i+1] = cg.boolConstraints[constraint_i].string()
	}
	fmt.Fprintln(w, strings.Join(constraints, "\t"))
}
