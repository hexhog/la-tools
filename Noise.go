package main

import (
	"fmt"
	"math/rand"
)

type Noise struct {
	ratio float64
}

func (n *Noise) addNoise(mean, _range float64) float64 {
	r := rand.New(rand.NewSource(99))
	noise := r.Float64()

	noiseRange := ratio * _range
	noisy := mean + (noise * noiseRange) - (noiseRange / 2)

	fmt.Printf("Range %v NoiseRange: %v value: %v noisy: %v \n", _range, noiseRange, mean, noisy)

	return noisy
}

func NewNoise(ratio float64) *Noise {
	return &Noise{
		ratio: ratio,
	}
}
