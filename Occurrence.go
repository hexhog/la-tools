package main

type Occurrence struct {
	factorList   []int
	factorList_n int
	count        int
	magnitude    float64
	list         []Occurrence
}

type ByCount []Occurrence

func (bc ByCount) Len() int           { return len(bc) }
func (bc ByCount) Less(i, j int) bool { return bc[i].count < bc[j].count }
func (bc ByCount) Swap(i, j int)      { bc[i], bc[j] = bc[j], bc[i] }
