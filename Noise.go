package main

import (
	"fmt"
	"math/rand"
)

type Noise struct {
	ratio float32
}

func (n *Noise) addNoise(mean, _range float32) float32 {
	r := rand.New(rand.NewSource(99))
	noise = r.Float32()

	noiseRange := ratio * _range
	noisy := mean + (noise * noiseRange) - (noiseRange / 2)

	fmt.Printf("Range %v NoiseRange: %v value: %v noisy: %v \n", _range, noiseRange, mean, noisy)

	return noisy
}

func NewNoise(ratio float32) *Noise {
	return &Noise{
		ratio: ratio,
	}
}
