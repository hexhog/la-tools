package main

import (
	"fmt"
	"log"
	"math"
)

// QR and least squares workspace
type WorkSpace struct {
	dataQ   [][]float64
	dataR   [][]float64
	workVec []float64
}

type TermIndex struct {
	termIndex int
	next      *TermIndex
}

var workSpace *WorkSpace // static

// static: setup the global workspace
func setupWorkSpace(rows, cols int) {
	// allocate initial memory
	workSpace = &WorkSpace{}

	// allocate memory for Q
	workSpace.dataQ = make([][]float64, rows)
	for row_i := 0; row_i < rows; row_i++ {
		workSpace.dataQ[row_i] = make([]float64, cols)
	}

	// allocate memory for R
	workSpace.dataR = make([][]float64, cols)
	for row_i := 0; row_i < cols; row_i++ {
		workSpace.dataR[row_i] = make([]float64, cols)
	}

	// allocate memory for work vector
	workSpace.workVec = make([]float64, rows)
}

type Model struct {
	csMatrix      *CSMatrix  // CSMatrix to work with
	maxTerms      int        // the maximum number of terms his model can accomodate
	tests         int        // the number of rows or tests
	terms         int        // number of terms in the model
	hTermIndex    *TermIndex // list with term indices
	response      *VectorXf  // response vector (this is actual data from the response file)
	modelResponse []float64  // actual response of the model
	coefVec       []float64  // n by 1
	resiVec       []float64  // m by 1
	rSquared      float64
}

// constructer - initialize the model
func newModel(response *VectorXf, maxTerms int, csMatrix *CSMatrix) *Model {
	m := &Model{
		csMatrix: csMatrix,
		maxTerms: maxTerms,
		response: response,
		tests:    response.getLength(),
	}

	// allocate memory for the coefficients vector
	m.coefVec = make([]float64, maxTerms)

	// allocate memory for the residuals vector
	m.resiVec = make([]float64, m.tests)

	// allocate memory for the model response
	m.modelResponse = make([]float64, m.tests)

	// add the intercept
	m.terms = 1
	m.hTermIndex = &TermIndex{}
	m.hTermIndex.termIndex = 0
	m.hTermIndex.next = nil

	// this is a one line way to get the same intercept above
	m.leastSquares()

	return m
}

// constructor - duplicate a model
func newDuplicateModel(model *Model) *Model {
	m := &Model{
		csMatrix: model.csMatrix,
		maxTerms: model.maxTerms,
		response: model.response,
		rSquared: model.rSquared,
		terms:    model.terms,
		tests:    model.response.getLength(),
	}

	// allocate memory for the coefficients vector and copy
	coefVec := make([]float64, m.maxTerms)
	for coef_i := 0; coef_i < m.maxTerms; coef_i++ {
		coefVec[coef_i] = model.coefVec[coef_i]
	}
	m.coefVec = coefVec

	// allocate memory for the residuals vector
	resiVec := make([]float64, m.tests)
	for resi_i := 0; resi_i < m.tests; resi_i++ {
		resiVec[resi_i] = model.resiVec[resi_i]
	}
	m.resiVec = resiVec

	// allocate memory for the model response
	modelResponse := make([]float64, m.tests)
	for resp_i := 0; resp_i < m.tests; resp_i++ {
		modelResponse[resp_i] = model.modelResponse[resp_i]
	}
	m.modelResponse = modelResponse

	// copy term index list
	destTermIndex := &m.hTermIndex
	for pTermIndex := model.hTermIndex; pTermIndex != nil; pTermIndex = pTermIndex.next {
		(*destTermIndex) = &TermIndex{}
		(*destTermIndex).termIndex = pTermIndex.termIndex
		(*destTermIndex).next = nil
		destTermIndex = &(*destTermIndex).next
	}

	return m
}

// print the model factors' names
func (m *Model) printModelFactors() {
	// var factor1i, factor2i, level1, level2 int

	// print out the number of terms
	fmt.Println("Model with ", m.terms, " terms")

	// print out the headers
	fmt.Println("Coefficient | Term")
	fmt.Println("------- | -------")

	// print out each term
	term_i := 0
	for pTermIndex := m.hTermIndex; pTermIndex != nil; pTermIndex = pTermIndex.next {
		// print out the coefficient
		termIndex := pTermIndex.termIndex

		// print out the factor names and level names
		fmt.Println(m.coefVec[term_i], " | ", m.csMatrix.getColName(m.csMatrix.getCol(termIndex)))
		term_i++
	}

	// calculate adjusted r-squared (terms - 2 to not include the intercept)
	adjustedRSquared := 1 - (1-m.rSquared)*(float64(m.tests)-1)/float64(m.tests-m.terms-2)

	// print model r-squared
	if m.rSquared == 1 {
		fmt.Println("Perfect Model!!!")
	}
	fmt.Println("R-Squared: ", m.rSquared)
	fmt.Println("Adjusted R-Squared: ", adjustedRSquared)
}

// perform least squares on this model
func (m *Model) leastSquares() {
	// used when accessing CS Matrix
	var csCol *CSCol

	// used when looping through term indices
	var pTermIndex *TermIndex

	// check dimensions
	if m.terms > m.tests {
		log.Fatal("Cols (terms) cannot be more than rows (tests) to perform QR")
	}

	// do QR here (column by column)
	pTermIndex = m.hTermIndex
	for col_i := 0; col_i < m.terms; col_i++ {

		csCol = m.csMatrix.getCol(pTermIndex.termIndex)
		// assign initial column A[col_i] to work vector
		for row_i := 0; row_i < m.tests; row_i++ {
			workSpace.workVec[row_i] = float64(csCol.dataP[row_i])
		}

		// subtract appropriate other vectors
		var dotProd float64
		for row_i := 0; row_i < col_i; row_i++ {
			// find the dot product of A[:][col_i] and Q[:][row_i]
			dotProd = 0
			for dotrow_i := 0; dotrow_i < m.tests; dotrow_i++ {
				dotProd += float64(csCol.dataP[dotrow_i]) * workSpace.dataQ[dotrow_i][row_i]
			}

			// assign the dot product to the R matrix
			workSpace.dataR[row_i][col_i] = dotProd

			// perform the necessary subtraction
			for subrow_i := 0; subrow_i < m.tests; subrow_i++ {
				workSpace.workVec[subrow_i] -= dotProd * workSpace.dataQ[subrow_i][row_i]
			}

		}

		// find the dot product of workVec and workVec
		dotProd = 0
		for dotrow_i := 0; dotrow_i < m.tests; dotrow_i++ {
			dotProd += workSpace.workVec[dotrow_i] * workSpace.workVec[dotrow_i]
		}
		dotProd = math.Sqrt(dotProd)

		// assign the dot product to the R matrix
		workSpace.dataR[col_i][col_i] = dotProd

		// assign the normalized column to Q
		if dotProd == 0 {
			for row_i := 0; row_i < m.tests; row_i++ {
				workSpace.dataQ[row_i][col_i] = 0
			}
		} else {
			for row_i := 0; row_i < m.tests; row_i++ {
				workSpace.dataQ[row_i][col_i] = workSpace.workVec[row_i] / dotProd
			}
		}

		pTermIndex = pTermIndex.next
	}

	/*
		cout << "Q matrix" << endl;
		for (int row_i = 0; row_i < tests; row_i++) {
			for (int col_i = 0; col_i < terms; col_i++) {
				cout << workSpace->dataQ[row_i][col_i] << "\t";
			}
			cout << endl;
		}

		cout << "R matrix" << endl;
		for (int row_i = 0; row_i < terms; row_i++) {
			for (int col_i = 0; col_i < terms; col_i++) {
				cout << workSpace->dataR[row_i][col_i] << "\t";
			}
			cout << endl;
		}
	*/

	/* m >> n
	** A is (m by n)
	** Q is (m by n)
	** R is (n by n)
	** b is (m by 1)
	** Q-transpose is (n by m)
	**
	 */
	// We now have QRx = b

	// perform workVec = Q-transpose * b

	// first zero out the work vector up until n (or cols)
	for row_i := 0; row_i < m.terms; row_i++ {
		workSpace.workVec[row_i] = 0
	}

	// we go column 1st because Q is transposed
	responseData := m.response.getData()
	for col_i := 0; col_i < m.tests; col_i++ {
		for row_i := 0; row_i < m.terms; row_i++ {
			// row and col are swapped because of the transpose
			workSpace.workVec[row_i] += workSpace.dataQ[col_i][row_i] * responseData[col_i]
		}
	}

	// We now have Rx = workVec but R is upper triangular (n by n)
	var rowSolution float64
	for row_i := m.terms - 1; row_i >= 0; row_i-- {

		// initialize row solution
		rowSolution = workSpace.workVec[row_i]

		// subtract other parts of LHS of equation
		for col_i := m.terms - 1; col_i > row_i; col_i-- {
			rowSolution -= workSpace.dataR[row_i][col_i] * m.coefVec[col_i]
		}

		// divide for final row solution
		if workSpace.dataR[row_i][row_i] == 0 {
			m.coefVec[row_i] = 0
		} else {
			m.coefVec[row_i] = rowSolution / workSpace.dataR[row_i][row_i]
		}

	}

	// coefVec holds the least squares coefficients
	//cout << "coefficient vector" << endl;
	//for (int term_i = 0; term_i < terms; term_i++) {
	//	cout << coefVec[term_i] << " ";
	//}
	//cout << endl;

	// calculate residuals and r-squared
	var SSres float64 = 0
	for row_i := 0; row_i < m.tests; row_i++ {
		m.modelResponse[row_i] = 0
	}

	// find model response
	pTermIndex = m.hTermIndex
	for term_i := 0; term_i < m.terms; term_i++ {
		csCol = m.csMatrix.getCol(pTermIndex.termIndex)
		for row_i := 0; row_i < m.tests; row_i++ {
			m.modelResponse[row_i] += float64(csCol.dataP[row_i]) * m.coefVec[term_i]
		}
		pTermIndex = pTermIndex.next
	}

	// find residuals and SSres
	for row_i := 0; row_i < m.tests; row_i++ {
		m.resiVec[row_i] = responseData[row_i] - m.modelResponse[row_i]

		// add the square of residual to SSres
		SSres += m.resiVec[row_i] * m.resiVec[row_i]
	}

	// find r-squared
	m.rSquared = 1 - (SSres / m.response.getSStot())

}

// get the residuals vector for this model
func (m *Model) getResiVec() []float64 {
	return m.resiVec
}

// get r-squared
func (m *Model) getRSquared() float64 {
	return m.rSquared
}

// check if term is part of the model
func (m *Model) termExists(col_i int) bool {
	for pTermIndex := m.hTermIndex; pTermIndex != nil; pTermIndex = pTermIndex.next {
		if pTermIndex.termIndex == col_i {
			// term was found!
			return true
		}
	}
	// term never found
	return false
}

// add a term (col_i of csMatrix) to the model. Returns true if added successfully
func (m *Model) addTerm(col_i int) bool {
	var pTermIndex **TermIndex
	for pTermIndex = &m.hTermIndex; *pTermIndex != nil; pTermIndex = &(*pTermIndex).next {

		if (*pTermIndex).termIndex == col_i {

			// already exists so cannot insert
			return false

		} else if (*pTermIndex).termIndex > col_i {

			// insert within the list
			termIndex := &TermIndex{}
			termIndex.termIndex = col_i
			termIndex.next = (*pTermIndex)
			(*pTermIndex) = termIndex
			m.terms++
			return true

		}

	}

	// create a new term index at the end
	termIndex := &TermIndex{}
	termIndex.termIndex = col_i
	termIndex.next = nil
	(*pTermIndex) = termIndex
	m.terms++

	return true
}

// remove a term from the model
func (m *Model) removeTerm(col_i int) bool {
	// loop through the term index list
	for pTermIndex := &m.hTermIndex; *pTermIndex != nil; pTermIndex = &(*pTermIndex).next {
		// check if it equals the term to be removed (col_i)
		if (*pTermIndex).termIndex == col_i {
			// take the current term index out of the list
			removed := (*pTermIndex)
			(*pTermIndex) = removed.next

			// decrease the number of terms
			m.terms--

			return true
		}
	}

	// term to remove was not found
	return false
}

// get the number of terms
func (m *Model) getTerms() int {
	return m.terms
}

// count occurrences in model
func (m *Model) countOccurrences(occurrence *Occurrence) {
	// count occurrences for each term
	term_i := 0
	for pTermIndex := m.hTermIndex; pTermIndex != nil; pTermIndex = pTermIndex.next {
		termIndex := pTermIndex.termIndex
		magnitude := m.coefVec[term_i+1]
		m.csMatrix.countOccurrences(m.csMatrix.getCol(termIndex), occurrence, 0, magnitude)
	}
}

// check if the model is a duplicate
func (m *Model) isDuplicate(n *Model) bool {
	// check if all terms match
	p1TermIndex := m.hTermIndex
	p2TermIndex := n.hTermIndex

	for !(p1TermIndex == nil && p2TermIndex == nil) {
		if p1TermIndex.termIndex != p2TermIndex.termIndex {
			return false
		}

		p1TermIndex = p1TermIndex.next
		p2TermIndex = p2TermIndex.next
	}

	return true
}
