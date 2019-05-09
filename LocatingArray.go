package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
)

const laVersion = "v2.0"

type GroupingInfo struct {
	levels        int              // total levels for this factor
	grouped       bool             //
	levelGroups   []int            //
	conGroup      *ConstraintGroup // constraint group to which this factor belongs
	conGroupIndex int              // this factor's index in the constraint group above
}

type LocatingArray struct {
	factorGrouping []*GroupingInfo
	levels         [][]int            // pointer to array of test levels (the main locating array)
	tests          int                // count of tests in locating array
	factors        int                // count of factors in locating array
	t              int                // covers t-way interactions
	nConGroups     int                //
	conGroups      []*ConstraintGroup //
	factorData     *FactorData        //
}

func NewLocatingArrayFromFactors(factors int, levelCounts []int) *LocatingArray {
	l := &LocatingArray{
		tests:      0,
		factors:    factors,
		t:          1,
		nConGroups: 0,
		conGroups:  make([]*ConstraintGroup, nConGroups),
	}

	factorGrouping := make([]*GroupingInfo, factors)

	// load the level counts and grouping for the factors
	for factor_i := 0; factor_i < factors; factor_i++ {
		factorGrouping[factor_i] = &GroupingInfo{
			levels:        levelCounts[factor_i],
			grouped:       false,
			conGroupIndex: -1,
		}
	}

	l.factorGrouping = factorGrouping

	return l
}

func NewLocatingArrayFromFile(filePath, factorDataFile string) *LocatingArray {
	l := &LocatingArray{}

	file, err := os.Open(filePath)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	var (
		tempData   int
		tempString string
	)

	scanner := bufio.NewScanner(file)
	i := 0
	for scanner.Scan() {
		// verify version
		if i == 0 {
			tempString = scanner.Text()
			if tempString != laVersion {
				log.Fatalf("LA version must be %s for this software but instead found %s", laVersion, tempString)
			}
		} else if i == 1 {
			// read the first 2 lines of locating array
			words := strings.Fields(scanner.Text())

			tests, err := strconv.Atoi(words[0])
			if err != nil {
				log.Fatal(err)
			}
			l.tests = tests

			factors, err := strconv.Atoi(words[1])
			if err != nil {
				log.Fatal(err)
			}
			l.factors = factors

			l.t = 2

			fmt.Printf("Tests %d Factors %d", l.tests, l.factors)

			l.factorGrouping = make([]*GroupingInfo, l.factors)
		} else if i == 2 {
			words := strings.Fields(scanner.Text())

			// load the level counts for the factors
			for factor_i := 0; factor_i < factors; factor_i++ {
				l.factorGrouping[factor_i] = &GroupingInfo{}

				levels, err := strconv.Atoi(words[factor_i])
				if err != nil {
					log.Fatal(err)
				}

				l.factorGrouping[factor_i].levels = levels
			}
		} else if 3 <= i <= 2+factors {
			// load the grouping for the factors
			words := strings.Fields(scanner.Text())

			// load grouped bool
			grouped, err := strconv.ParseBool(words[0])
			if err != nil {
				log.Fatal(err)
			}
			l.factorGrouping[factor_i].grouped = grouped

			if l.factorGrouping[factor_i].grouped {
				// allocate array for level grouping
				l.factorGrouping[factor_i].levelGroups = make([]byte, l.factorGrouping[factor_i].levels)

				for level_i := 0; level_i < l.factorGrouping[factor_i].levels; level_i++ {
					// load level group
					levelGroup, err := strconv.Atoi(words[level_i+1])
					if err != nil {
						log.Fatal(err)
					}
					l.factorGrouping[factor_i].levelGroups[level_i] = byte(levelGroup)
				}
			}

			l.factorGrouping[factor_i].conGroupIndex = -1

			// load factor data (must be done before constraint groups)
			if factorDataFile == "" {
				l.factorData = NewFactorDataFromArray(l.getGroupingInfo(), l.getFactors())
			} else {
				l.factorData = NewFactorDataFromFile(factorDataFile)
			}
		} else if i == 3+factors {
			// load constraint groups
			nConGroups, err := strconv.Atoi(scanner.Text())
			if err != nil {
				log.Fatal(err)
			}
			l.conGroups = make([]*ConstraintGroup, nConGroups)

			// load each constraint group
			for iConGroup := 0; iConGroup < nConGroups; iConGroup++ {
				words := []string{scanner.Text(), scanner.Text()}
				l.conGroups[iConGroup] = newConstraintGroup(l, words)
			}
		} else {
			// load the tests now
			for test_i := 0; test_i < l.tests; test_i++ {
				levelRow := make([]byte, l.factors)

				words := strings.Fields(scanner.Text())

				// now load each factor level for this specific test
				for factor_i := 0; factor_i < factors; factor_i++ {
					factor, err := strconv.Atoi(words[factor_i])
					if err != nil {
						log.Fatal(err)
					}
					levelRow[factor_i] = factor
				}

				l.addLevelRow(levelRow)
				l.tests--
			}
		}
		i++
	}
}

func (l *LocatingArray) addLevelRow(levelRow []byte) {
	l.levels = append(l.levels, levelRow)
	l.tests++
}

func (l *LocatingArray) remLevelRow() []byte {
	levelRow := l.levels[len(l.levels)-1]
	l.levels = l.levels[:len(l.levels)-1]
	l.tests--
	return levelRow
}

func (l *LocatingArray) getGroupingInfo() []*GroupingInfo {
	return l.factorGrouping
}

func (l *LocatingArray) getLevelMatrix() [][]int {
	return l.levels
}

func (l *LocatingArray) getFactors() int {
	return l.factors
}

func (l *LocatingArray) getTests() int {
	return l.tests
}

func (l *LocatingArray) getT() int {
	return l.t
}

func (l *LocatingArray) getNConGroups() int {
	return l.nConGroups
}

func (l *LocatingArray) getConGroups() []*ConstraintGroup {
	return l.conGroups
}

func (l *LocatingArray) getFactorData() *FactorData {
	return l.factorData
}

func (l *LocatingArray) writeToFile(filePath string) {
	file, err := os.Create(filePath)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	w := bufio.NewWriter(file)

	// write version
	fmt.Fprintln(w, laVersion)

	// initial tests and factors
	fmt.Fprintln(w, fmt.Sprintf("%d \t %d", l.tests, l.factors))

	// write the level counts for the factors
	levelCounts := make([]int, l.factors)
	for factor_i := 0; factor_i < factors; factor_i++ {
		levelCounts[factor_i] = l.factorGrouping[factor_i].levels
	}
	fmt.Fprintln(w, strings.Join(levelCounts, "\t"))

	// write the grouping for the factors
	for factor_i = 0; factor_i < factors; factor_i++ {
		levelGroups := make([]int, l.factorGrouping[factor_i].levels+1)

		// load grouped bool
		levelGroups[0] = l.factorGrouping[factor_i].grouped

		// check grouping
		if l.factorGrouping[factor_i].grouped {
			for level_i := 0; level_i < factorGrouping[factor_i].levels; level_i++ {
				// write level group
				levelGroups[level_i+1] = l.factorGrouping[factor_i].levelGroups[level_i]
			}
		}
		fmt.Fprintln(w, strings.Join(levelGroups, "\t"))
	}

	// write constraint groups
	fmt.Fprintln(w, l.nConGroups)
	for iConGroup := 0; iConGroup < l.nConGroups; iConGroup++ {
		// l.conGroups[iConGroup]->writeToStream(ofs);
	}

	// write the tests now
	levelMatrix := getLevelMatrix()
	for test_i := 0; test_i < tests; test_i++ {
		levelFactors := make([]byte, l.factors)
		// now write each factor level for this specific test
		for factor_i = 0; factor_i < factors; factor_i++ {
			levelFactors[factor_i] = levelMatrix[test_i][factor_i]
		}
		fmt.Fprintln(w, strings.Join(levelFactors, "\t"))
	}

	w.Flush()
}

// func main() {
// 	l = NewLocatingArrayFromFile("la.tsv", "fd.tsv")
// 	_ = l
// }
