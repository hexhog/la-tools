package main

import (
	"fmt"
)

type VectorXf struct {
	length    int
	data      []float64
	usesSStot bool    // indicates whether SStot has been calculated
	SStot     float64 // SStot for calculating r-squared
}

func NewVectorXf(length int) *VectorXf {
	return &VectorXf{
		length:    length,
		data:      make([]float, length),
		usesSStot: false,
	}
}

func (v *VectorXf) getData() []float64 {
	return v.data
}

func (v *VectorXf) getLength() int {
	return v.length
}

func (v *VectorXf) calculateSStot() {
	// calculate SStot for r-squared

	// 1st get yBar
	var yBar, ySum float64
	for row_i := 0; row_i < v.length; row_i++ {
		ySum += v.data[row_i]
	}
	yBar = ySum / length

	// 2nd get SStot
	v.SStot = 0
	for row_i = 0; row_i < length; row_i++ {
		v.SStot += (v.data[row_i] - yBar) * (v.data[row_i] - yBar)
	}

	v.usesSStot = true
}

func (v *VectorXf) getSStot() float64 {
	// verify SStot has been calculated
	if !v.usesSStot {
		fmt.Println("You have not calculated SStot yet")
		return 0
	}

	return v.SStot
}
