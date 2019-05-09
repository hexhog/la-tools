package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type Factor struct {
	name        string    // name of the factor
	levels      int       // level count
	numeric     bool      // is numeric
	levelNames  []string  // names of the levels
	levelValues []float64 // values of the levels (NULL if numeric is false)
}

type FactorData struct {
	factorCount int       // total factors (75-100ish for large dataset)
	factors     []*Factor // array of factor pointers
}

func NewFactor() *Factor {
	return &Factor{}
}

func NewFactorDataFromArray(arrayInfo []*GroupingInfo, factorCount int) *FactorData {
	// allocate necessary memory
	factors := make([]*Factor, factorCount)

	// read in line by line
	for factor_i := 0; factor_i < factorCount; factor_i++ {
		// allocate the factor
		factor := NewFactor()

		// read the factor name
		factor.name = fmt.Sprintf("F%d", factor_i)

		// read the number of levels
		factor.levels = arrayInfo[factor_i].levels

		// read whether the factor has numeric levels
		factor.numeric = false

		// allocate memory for each level
		factor.levelNames = make([]string, factor.levels)

		// read individual level names
		for level_i := 0; level_i < factor.levels; level_i++ {
			factor.levelNames[level_i] = fmt.Sprintf("L%d", level_i)
		}

		factors[factor_i] = factor
	}

	return &FactorData{factorCount, factors}
}

func NewFactorDataFromFile(filePath string) *FactorData {
	file, err := os.Open(filePath)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	var (
		factorCount int
		factors     []*Factor
		e           error
	)

	scanner := bufio.NewScanner(file)
	i := -1
	for scanner.Scan() {
		if i == -1 {
			// grab the total number of factors
			factorCount, e = strconv.Atoi(scanner.Text())
			if e != nil {
				log.Fatal(e)
			}

			// allocate necessary memory
			factors = make([]*Factor, factorCount)
		} else {
			// read in line by line
			words := strings.Fields(scanner.Text())

			// allocate the factor
			factor := NewFactor()

			// read the factor name
			factor.name = words[0]

			// read the number of levels
			levels, err := strconv.Atoi(words[1])
			if err != nil {
				log.Fatal(err)
			}
			factor.levels = levels

			// read whether the factor has numeric levels
			numeric, err := strconv.ParseBool(words[2])
			if err != nil {
				log.Fatal(err)
			}
			factor.numeric = numeric

			// allocate memory for each level
			factor.levelNames = make([]string, factor.levels)

			// read individual level names
			factor.levelNames = words[3 : 3+factor.levels]

			// read the level values if they are numeric
			if factor.numeric {
				// allocate memory for each level value
				factor.levelValues = make([]float64, factor.levels)

				// read individual level values
				factor.levelNames = words[4+factor.levels : 4+2*factor.levels]
			}

			factors[i] = factor
		}
		i++
	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	return &FactorData{factorCount, factors}
}

func (f *FactorData) getFactor(factor_i int) *Factor {
	return f.factors[factor_i]
}

func (f *FactorData) getFactorName(factor_i int) string {
	return f.factors[factor_i].name
}

func (f *FactorData) getFactorLevelName(factor_i, level_i int) string {
	return f.factors[factor_i].levelNames[level_i]
}

func (f *FactorData) getNumericFactorLevel(factor_i, level_i int) float64 {
	if f.factors[factor_i].numeric {
		return f.factors[factor_i].levelValues[level_i]
	} else {
		return 0
	}
}
