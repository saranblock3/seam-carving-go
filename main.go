package main

import (
	"fmt"
	"image"
	"image/color"
	"image/jpeg"
	"log"
	"math"
	"os"
	"time"
)

// func fillEnergyMatrix(m image.Image) [][]float64 {
// 	start := time.Now()
// 	energyMatrix := make([][]float64, m.Bounds().Max.Y)
// 	for y := 0; y < m.Bounds().Max.Y; y++ {
// 		energyRow := make([]float64, m.Bounds().Max.X)
// 		for x := 0; x < m.Bounds().Max.X; x++ {
// 			energyRow[x] = calculatePixelEnergy(&m, x, y)
// 		}
// 		energyMatrix[y] = energyRow
// 	}
// 	fmt.Println("energy calculations:", time.Since(start))
// 	return energyMatrix
// }

// func updateEnergyMatrix(m image.Image, energyMatrix [][]float64, seam []point) [][]float64 {
// 	start := time.Now()
// 	newEnergyMatrix := make([][]float64, m.Bounds().Max.Y)
// 	for y := 0; y < m.Bounds().Max.Y; y++ {
// 		newRow := make([]float64, m.Bounds().Max.X)
// 		seamFound := false
// 		seamBehind := false
// 		for x := 0; x < len(energyMatrix[0]); x++ {
// 			if x == seam[y].x {
// 				seamFound = true
// 				continue
// 			}
// 			if seamFound {
// 				newRow[x-2] = calculatePixelEnergy(&m, x-2, y)
// 				newRow[x-1] = calculatePixelEnergy(&m, x, y)
// 				seamFound = false
// 				seamBehind = true
// 			} else if seamBehind {
// 				newRow[x-1] = energyMatrix[y][x]
// 			} else {
// 				newRow[x] = energyMatrix[y][x]
// 			}
// 		}
// 		newEnergyMatrix[y] = newRow
// 	}
// 	fmt.Println("energy calculations:", time.Since(start))
// 	return newEnergyMatrix
// }

// func fillEnergyMatrixAsync2(m image.Image) [][]float64 {
// 	start := time.Now()
// 	energyMatrix := make([][]float64, m.Bounds().Max.Y)
// 	ch := make(chan energyRowEntry)
// 	for y := 0; y < m.Bounds().Max.Y; y++ {
// 		go func(y int, m image.Image, ch chan energyRowEntry) {
// 			energyRow := make([]float64, m.Bounds().Max.X)
// 			for x := 0; x < m.Bounds().Max.X; x++ {
// 				energyRow[x] = calculatePixelEnergy(&m, x, y)
// 			}
// 			ch <- energyRowEntry{y, energyRow}
// 		}(y, m, ch)
// 	}
// 	for y := 0; y < m.Bounds().Max.Y; y++ {
// 		energyRow := <-ch
// 		energyMatrix[energyRow.y] = energyRow.energyRow
// 	}
// 	fmt.Println("energy calculations:", time.Since(start))
// 	// fmt.Println(energyMatrix)
// 	return energyMatrix
// }

// type seamEntry struct {
// 	index int
// 	seam  []point
// 	score float64
// }

// func findSeams(energyMatrix [][]float64) pkg.PriorityQueue {
// 	start := time.Now()
// 	seams := make(pkg.PriorityQueue, len(energyMatrix[0])-2)
// 	// seamScores := make(map[int]float64, len(energyMatrix[0])-2)
// 	ch := make(chan seamEntry)
// 	for x := 1; x < len(energyMatrix[0])-1; x++ {
// 		go func(x int, energyMatrix [][]float64, ch chan seamEntry) {
// 			seam := make([]point, len(energyMatrix))
// 			seam[0] = point{x, 0}
// 			seamScore := 0.0
// 			xSeam := x
// 			for y := 1; y < len(energyMatrix); y++ {
// 				// fmt.Println(xSeam)
// 				if xSeam-1 < 1 {
// 					energyMock := (energyMatrix[y][xSeam] + energyMatrix[y][xSeam+1]) / 2
// 					seamScore += varianceEnergies(energyMock, energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 					xSeam += minEnergy(math.Inf(1), energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 				} else if xSeam+1 >= len(energyMatrix[0])-1 {
// 					energyMock := (energyMatrix[y][xSeam-1] + energyMatrix[y][xSeam]) / 2
// 					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMock)
// 					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], math.Inf(1))
// 				} else {
// 					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 				}
// 				seam[y] = point{xSeam, y}
// 			}
// 			ch <- seamEntry{x, seam, seamScore}
// 		}(x, energyMatrix, ch)
// 		// seams[x-1] = seam
// 	}
// 	for x := 1; x < len(energyMatrix[0])-1; x++ {
// 		currentSeamEntry := <-ch
// 		item := pkg.NewItem(&currentSeamEntry.seam, currentSeamEntry.score, currentSeamEntry.index)
// 		seams[x-1] = item
// 		// seams[currentSeamEntry.index-1] = currentSeamEntry.seam
// 		// seamScores[currentSeamEntry.index-1] = currentSeamEntry.score
// 	}
// 	heap.Init(&seams)
// 	// minScore := math.Inf(1)
// 	// minSeamIdx := 0
// 	// for seamIdx, score := range seamScores {
// 	// 	if score < minScore {
// 	// 		minSeamIdx = seamIdx
// 	// 		minScore = score
// 	// 	}
// 	// }
// 	fmt.Println("seam score calulations:", time.Since(start))
// 	return seams
// }

// func scoreSeams(energyMatrix [][]float64) []point {
// 	start := time.Now()
// 	seams := make([][]point, len(energyMatrix[0])-2)
// 	seamScores := make(map[int]float64, len(energyMatrix[0])-2)
// 	ch := make(chan seamEntry)
// 	for x := 1; x < len(energyMatrix[0])-1; x++ {
// 		go func(x int, energyMatrix [][]float64, ch chan seamEntry) {
// 			seam := make([]point, len(energyMatrix))
// 			seam[0] = point{x, 0}
// 			seamScore := 0.0
// 			xSeam := x
// 			for y := 1; y < len(energyMatrix); y++ {
// 				// fmt.Println(xSeam)
// 				if xSeam-1 < 1 {
// 					energyMock := (energyMatrix[y][xSeam] + energyMatrix[y][xSeam+1]) / 2
// 					seamScore += varianceEnergies(energyMock, energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 					xSeam += minEnergy(math.Inf(1), energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 				} else if xSeam+1 >= len(energyMatrix[0])-1 {
// 					energyMock := (energyMatrix[y][xSeam-1] + energyMatrix[y][xSeam]) / 2
// 					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMock)
// 					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], math.Inf(1))
// 				} else {
// 					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
// 				}
// 				seam[y] = point{xSeam, y}
// 			}
// 			ch <- seamEntry{x, seam, seamScore}
// 		}(x, energyMatrix, ch)
// 		// seams[x-1] = seam
// 	}
// 	for x := 1; x < len(energyMatrix[0])-1; x++ {
// 		currentSeamEntry := <-ch
// 		seams[currentSeamEntry.index-1] = currentSeamEntry.seam
// 		seamScores[currentSeamEntry.index-1] = currentSeamEntry.score
// 	}
// 	minScore := math.Inf(1)
// 	minSeamIdx := 0
// 	for seamIdx, score := range seamScores {
// 		if score < minScore {
// 			minSeamIdx = seamIdx
// 			minScore = score
// 		}
// 	}
// 	fmt.Println("score calulations:", time.Since(start))
// 	return seams[minSeamIdx]
// }

// func varianceEnergies(energies ...float64) float64 {
// 	return meanSqaureEnergies(energies...) - math.Pow(meanEnergies(energies...), 2)
// }

// func meanEnergies(energies ...float64) float64 {
// 	sum := 0.0
// 	for _, v := range energies {
// 		sum += v
// 	}
// 	return sum / float64(len(energies))
// }

// func meanSqaureEnergies(energies ...float64) float64 {
// 	sum := 0.0
// 	for _, v := range energies {
// 		sum += math.Pow(v, 2)
// 	}
// 	return sum / float64(len(energies))
// }

// func minEnergy(energies ...float64) int {
// 	minEnergy := energies[0]
// 	minEnergyIdx := 0
// 	for i, v := range energies {
// 		if v < minEnergy {
// 			minEnergyIdx = i
// 		}
// 	}
// 	return minEnergyIdx - 1
// }

// func fillNewImage(m MyImage, seam []point) *MyImage {
// 	start := time.Now()
// 	newImage := image.NewRGBA(image.Rect(0, 0, m.Bounds().Max.X-1, m.Bounds().Max.Y))
// 	for y := 0; y < m.Bounds().Max.Y; y++ {
// 		seamFound := false
// 		for x := 0; x < m.Bounds().Max.X; x++ {
// 			if seam[y].x == x && seam[y].y == y {
// 				seamFound = true
// 				continue
// 			}
// 			if seamFound {
// 				newImage.Set(x-1, y, m.At(x, y))
// 			} else {
// 				newImage.Set(x, y, m.At(x, y))
// 			}
// 		}
// 	}
// 	fmt.Println("image calulations:", time.Since(start))
// 	return &MyImage{newImage, map[image.Point]color.Color{}}
// }

// func fillNewImage3(m MyImage, seam []point) *MyImage {
// 	start := time.Now()
// 	newImage := image.NewRGBA(image.Rect(0, 0, m.Bounds().Max.X-1, m.Bounds().Max.Y))
// 	for y := 0; y < m.Bounds().Max.Y; y++ {
// 		seamFound := false
// 		for x := 0; x < m.Bounds().Max.X; x++ {
// 			if seam[y].x == x && seam[y].y == y {
// 				seamFound = true
// 				continue
// 			}
// 			if seamFound {
// 				newImage.Set(x-1, y, m.At(x, y))
// 			} else {
// 				newImage.Set(x, y, m.At(x, y))
// 			}
// 		}
// 	}
// 	fmt.Println("image calulations:", time.Since(start))
// 	return &MyImage{newImage, map[image.Point]color.Color{}}
// }

// func getSeamsAsSlice(seams pkg.PriorityQueue, reduce int) [][]point {
// 	seamsSlice := make([][]point, reduce)
// 	for i := 0; i < reduce; i++ {
// 		seamsSlice[i] = *seams.Pop().(*pkg.Item).GetValue().(*[]point)
// 	}
// 	seamsSliceReflected := make([][]point, len(seamsSlice[0]))
// 	for y := 0; y < len(seamsSlice[0]); y++ {
// 		seamRow := make([]point, len(seamsSlice))
// 		for x := 0; x < len(seamsSlice); x++ {
// 			seamRow[x] = seamsSlice[x][y]
// 		}
// 		seamsSliceReflected[y] = seamRow
// 	}
// 	return seamsSliceReflected
// }

// func getXsFromPointSlice(pointSlice []point) []int {
// 	xSlice := make([]int, len(pointSlice))
// 	for i, p := range pointSlice {
// 		xSlice[i] = p.x
// 	}
// 	return xSlice
// }

// func fillNewImage2(m image.Image, seams pkg.PriorityQueue, reduce int) *MyImage {
// 	start := time.Now()
// 	newImage := image.NewRGBA(image.Rect(0, 0, m.Bounds().Max.X-reduce, m.Bounds().Max.Y))
// 	imageMatrix := getImage(m)
// 	seamSlice := getSeamsAsSlice(seams, reduce)
// 	// fmt.Println(seamSlice[0][0].x)
// 	for i, seamRow := range seamSlice {
// 		xRow := getXsFromPointSlice(seamRow)
// 		sort.Ints(xRow)
// 		// fmt.Println(xRow)
// 		imageRow := imageMatrix[i]
// 		bias := 0
// 		seamRowIdx := 0
// 		fmt.Println(len(imageRow))
// 		for j := 0; j < m.Bounds().Max.X; j++ {
// 			if seamRowIdx > 99 {
// 				break
// 			}
// 			// fmt.Println("j", j)
// 			// fmt.Println("x", xRow[seamRowIdx])
// 			// fmt.Println(xRow[seamRowIdx] == j)
// 			if xRow[seamRowIdx] == j {
// 				imageRow = slices.Delete(imageRow, j-bias, j+1-bias)
// 				// imageRow = append(imageRow[:xRow[seamRowIdx]-bias], imageRow[xRow[seamRowIdx]+1-bias:]...)
// 				bias += 1
// 				seamRowIdx += 1
// 				// j += 1
// 			}
// 		}
// 		fmt.Println(len(imageRow))
// 		fmt.Println()
// 		imageMatrix[i] = imageRow
// 	}
// 	for y := 0; y < m.Bounds().Max.Y; y++ {
// 		for x := 0; x < m.Bounds().Max.X-reduce; x++ {
// 			newImage.Set(x, y, imageMatrix[y][x])
// 		}
// 	}
// 	fmt.Println("image calulations:", time.Since(start))
// 	return &MyImage{newImage, map[image.Point]color.Color{}}
// }

// func findCrossingSeams(p [][]int, seam []point) []point {
// 	var crossingSeams []point
// 	for i, seamPoint := range seam {
// 		for j := -1; j < 2; j++ {
// 			x := seamPoint.x + j
// 			y := seamPoint.y + 1
// 			if y > len(seam)-1 {
// 				break
// 			}
// 			if x > len(p[0])-1 {
// 				continue
// 			}
// 			if x < 0 {
// 				continue
// 			}
// 			if seam[y].x == x {
// 				continue
// 			}
// 			if p[y][x] == -j {
// 				if x > seam[i+1].x {
// 					crossingSeams = append(crossingSeams, point{x - 1, y})
// 				} else {
// 					crossingSeams = append(crossingSeams, point{x, y})
// 				}
// 			}
// 		}
// 	}
// 	return crossingSeams
// }

// type rect struct {
// 	upper []point
// 	lower []point
// }

// func checkExclusion(x, y, exclusion rect) bool {
// 	return true
// }

// func findBestSeamAsync(xBound, yBound int, energyMatrix [][]float64) []point {
// 	start := time.Now()
// 	e := energyMatrix
// 	opt := make([][]float64, yBound)
// 	p := make([][]int, yBound)
// 	for i := 0; i < yBound; i++ {
// 		opt[i] = make([]float64, xBound)
// 		p[i] = make([]int, xBound)
// 	}
// 	for j := 0; j < xBound; j++ {
// 		opt[0][j] = e[0][j]
// 		p[0][j] = 0
// 	}
// 	for i := 1; i < yBound; i++ {
// 		for j := 0; j < xBound; j++ {
// 			opt[i][j] = opt[i-1][j]
// 			p[i][j] = 0
// 			if j > 1 {
// 				if opt[i-1][j-1] < opt[i][j] {
// 					opt[i][j] = opt[i-1][j-1]
// 					p[i][j] = -1
// 				}
// 			}
// 			if j < xBound-2 {
// 				if opt[i-1][j+1] < opt[i][j] {
// 					opt[i][j] = opt[i-1][j+1]
// 					p[i][j] = 1
// 				}
// 			}
// 			// fmt.Println(i, j)
// 			// fmt.Println(len(e))
// 			opt[i][j] = opt[i][j] + e[i][j]
// 		}
// 	}
// 	k := 1
// 	for j := 2; j < xBound-1; j++ {
// 		if opt[yBound-1][j] < opt[yBound-1][k] {
// 			k = j
// 		}
// 	}
// 	// return p, k
// 	seam := retraceSeam(p, k)
// 	// fmt.Println(seam)
// 	fmt.Println("seam calulations:", time.Since(start))
// 	return seam
// }

// func findSeam(xBound, yBound int, p [][]int, opt [][]float64) []point {
// 	k := 1
// 	for j := 2; j < xBound-1; j++ {
// 		if opt[yBound-1][j] < opt[yBound-1][k] {
// 			k = j
// 		}
// 	}
// 	seam := retraceSeam(p, k)
// 	return seam
// }

// func findSeamPaths(m image.Image, energyMatrix [][]float64) ([][]int, [][]float64) {
// 	e := energyMatrix
// 	opt := make([][]float64, m.Bounds().Max.Y)
// 	p := make([][]int, m.Bounds().Max.Y)
// 	for i := 0; i < m.Bounds().Max.Y; i++ {
// 		opt[i] = make([]float64, m.Bounds().Max.X)
// 		p[i] = make([]int, m.Bounds().Max.X)
// 	}
// 	for j := 0; j < m.Bounds().Max.X; j++ {
// 		opt[0][j] = e[0][j]
// 		p[0][j] = 0
// 	}
// 	for i := 1; i < m.Bounds().Max.Y; i++ {
// 		for j := 0; j < m.Bounds().Max.X; j++ {
// 			opt[i][j] = opt[i-1][j]
// 			p[i][j] = 0
// 			if j > 1 {
// 				if opt[i-1][j-1] < opt[i][j] {
// 					opt[i][j] = opt[i-1][j-1]
// 					p[i][j] = -1
// 				}
// 			}
// 			if j < m.Bounds().Max.X-2 {
// 				if opt[i-1][j+1] < opt[i][j] {
// 					opt[i][j] = opt[i-1][j+1]
// 					p[i][j] = 1
// 				}
// 			}
// 			opt[i][j] = opt[i][j] + e[i][j]
// 		}
// 	}
// 	return p, opt
// }

// func updateSeamPaths(xBound, yBound int, p [][]int, opt [][]float64, seam []point, energyMatrix [][]float64) ([][]int, [][]float64) {
// 	e := energyMatrix
// 	start := time.Now()
// 	crossingSeamPoints := findCrossingSeams(p, seam)
// 	fmt.Println("cross calculations:", time.Since(start))
// 	p = deleteSeamFromMatrix[int](p, seam)
// 	fmt.Println("p calculations:", time.Since(start))
// 	opt = deleteSeamFromMatrix[float64](opt, seam)
// 	fmt.Println("opt calculations:", time.Since(start))
// 	for _, crossingSeamPoint := range crossingSeamPoints {
// 		x := crossingSeamPoint.x
// 		y := crossingSeamPoint.y
// 		maxEnergy := math.Inf(-1)
// 		pIdx := 0
// 		for i := -1; i < 2; i++ {
// 			if e[y-1][x+i] >= maxEnergy && x+i >= 0 && x+i < len(e[0]) {
// 				maxEnergy = e[y-1][x+i]
// 				pIdx = i
// 			}
// 		}
// 		p[y][x] = pIdx
// 		// for y < len(e) {
// 		// 	opt[y][x] = opt[y-1][x+p[y][x]] + e[y][x]
// 		// 	x += p[]
// 		// 	y += 1
// 		// }
// 	}
// 	fmt.Println("p and opt calculations: ======", time.Since(start))
// 	// opt := make
// 	for j := 0; j < xBound; j++ {
// 		opt[0][j] = e[0][j]
// 	}
// 	for i := 1; i < yBound; i++ {
// 		for j := 0; j < xBound; j++ {
// 			opt[i][j] = opt[i-1][j]
// 			if j > 1 {
// 				if opt[i-1][j-1] < opt[i][j] {
// 					opt[i][j] = opt[i-1][j-1]
// 				}
// 			}
// 			if j < xBound-2 {
// 				if opt[i-1][j+1] < opt[i][j] {
// 					opt[i][j] = opt[i-1][j+1]
// 				}
// 			}
// 			opt[i][j] = opt[i][j] + e[i][j]
// 		}
// 	}

// 	fmt.Println("p and opt calculations:", time.Since(start))
// 	return p, opt
// }

// func deleteSeamFromMatrix[T any](matrix [][]T, seam []point) [][]T {
// 	for i := 0; i < len(matrix); i++ {
// 		var newRow []T
// 		newRow = matrix[i][:seam[i].x]
// 		newRow = append(newRow, matrix[i][seam[i].x+1:]...)
// 		matrix[i] = newRow
// 	}
// 	return matrix
// }

type Changeable interface {
	Set(x, y int, c color.Color)
}

type MyImage struct {
	image.Image
	custom map[image.Point]color.Color
}

type MyCarving struct {
	image.Image
	seam []point
}

type energyRowEntry struct {
	y         int
	energyRow []float64
}

func fillEnergyMatrixAsync(m image.Image) [][]float64 {
	// start := time.Now()
	energyMatrix := make([][]float64, m.Bounds().Max.Y)
	ch := make(chan energyRowEntry)
	for y := 0; y < m.Bounds().Max.Y; y++ {
		go func(y int, m *image.Image, ch chan energyRowEntry) {
			energyRow := make([]float64, (*m).Bounds().Max.X)
			for x := 0; x < (*m).Bounds().Max.X; x++ {
				energyRow[x] = calculatePixelEnergy(m, x, y)
			}
			ch <- energyRowEntry{y, energyRow}
		}(y, &m, ch)
	}
	for y := 0; y < m.Bounds().Max.Y; y++ {
		energyRow := <-ch
		energyMatrix[energyRow.y] = energyRow.energyRow
	}
	// fmt.Println("energy calculations:", time.Since(start))
	// fmt.Println(energyMatrix)
	return energyMatrix
}

func updateEnergyMatrix2(xBound, yBound int, energyMatrix [][]float64, seam []point) [][]float64 {
	newEnergyMatrix := make([][]float64, yBound)
	for y := 0; y < yBound; y++ {
		var newRow []float64
		// newRow = energyMatrix[y][:seam[y].x-1]
		// newRow = append(newRow, calculatePixelEnergy(&m, seam[y].x-1, y))
		// newRow = append(newRow, calculatePixelEnergy(&m, seam[y].x+1, y))
		// newRow = append(newRow, energyMatrix[y][seam[y].x+2:]...)
		newRow = energyMatrix[y][:seam[y].x]
		// newRow = append(newRow, calculatePixelEnergy(&m, seam[y].x-1, y))
		// newRow = append(newRow, calculatePixelEnergy(&m, seam[y].x+1, y))
		newRow = append(newRow, energyMatrix[y][seam[y].x+1:]...)
		newEnergyMatrix[y] = newRow
	}
	return newEnergyMatrix
}

type point struct {
	x int
	y int
}

func calculatePixelEnergy(m *image.Image, x, y int) float64 {
	squareError := 0.0
	for j := -1; j < 2; j++ {
		if y+j < (*m).Bounds().Min.Y {
			continue
		}
		if y+j >= (*m).Bounds().Max.Y {
			continue
		}
		for i := -1; i < 2; i++ {
			if x+i < (*m).Bounds().Min.X {
				continue
			}
			if x+i >= (*m).Bounds().Max.X {
				continue
			}
			if i == 0 && j == 0 {
				continue
			}
			neighbourX := x + i
			neighbourY := y + j
			c1 := (*m).At(x, y)
			c2 := (*m).At(neighbourX, neighbourY)
			squareError += squareErrorRGBA(c1, c2)
		}
	}
	meanSquareError := squareError / 8
	return meanSquareError
}

func squareErrorRGBA(c1, c2 color.Color) float64 {
	r1, g1, b1, a1 := c1.RGBA()
	r2, g2, b2, a2 := c2.RGBA()

	total := 0.0

	total += math.Pow(float64(r1-r2), 2)
	total += math.Pow(float64(g1-g2), 2)
	total += math.Pow(float64(b1-b2), 2)
	total += math.Pow(float64(a1-a2), 2)

	return total
}

func NewMyImg(img image.Image) *MyImage {
	return &MyImage{img, map[image.Point]color.Color{}}
}

func (m *MyImage) Set(x, y int, c color.Color) {
	m.custom[image.Point{x, y}] = c
}

func (m *MyImage) At(x, y int) color.Color {
	if c := m.custom[image.Point{x, y}]; c != nil {
		return c
	}
	return m.Image.At(x, y)
}

func fillNewImageMatrix(m [][]color.Color, seam []point) [][]color.Color {
	newImageMatrix := make([][]color.Color, len(m))
	for _, p := range seam {
		newImageMatrix[p.y] = m[p.y][:p.x]
		newImageMatrix[p.y] = append(newImageMatrix[p.y], m[p.y][p.x+1:]...)
	}
	return newImageMatrix
}

func fillImage(xBound, yBound int, imageMatrix [][]color.Color) *MyImage {
	newImage := image.NewRGBA(image.Rect(0, 0, xBound, yBound))
	for y := 0; y < yBound; y++ {
		for x := 0; x < xBound; x++ {
			newImage.Set(x, y, imageMatrix[y][x])
		}
	}
	return &MyImage{newImage, map[image.Point]color.Color{}}

}

func getImage(m image.Image) [][]color.Color {
	imageMatrix := make([][]color.Color, m.Bounds().Max.Y)
	for y := 0; y < m.Bounds().Max.Y; y++ {
		imageRow := make([]color.Color, m.Bounds().Max.X)
		for x := 0; x < m.Bounds().Max.X; x++ {
			imageRow[x] = m.At(x, y)
		}
		imageMatrix[y] = imageRow
	}
	return imageMatrix
}

func findBestSeam(xBound, yBound int, energyMatrix [][]float64) []point {
	e := energyMatrix
	opt := make([][]float64, yBound)
	p := make([][]int, yBound)
	for i := 0; i < yBound; i++ {
		opt[i] = make([]float64, xBound)
		p[i] = make([]int, xBound)
	}
	for j := 0; j < xBound; j++ {
		opt[0][j] = e[0][j]
		p[0][j] = 0
	}
	start := time.Now()
	for i := 1; i < yBound; i++ {
		for j := 0; j < xBound; j++ {
			opt[i][j] = opt[i-1][j]
			p[i][j] = 0
			if j > 1 {
				if opt[i-1][j-1] < opt[i][j] {
					opt[i][j] = opt[i-1][j-1]
					p[i][j] = -1
				}
			}
			if j < xBound-2 {
				if opt[i-1][j+1] < opt[i][j] {
					opt[i][j] = opt[i-1][j+1]
					p[i][j] = 1
				}
			}
			// fmt.Println(i, j)
			// fmt.Println(len(e))
			opt[i][j] = opt[i][j] + e[i][j]
		}
	}
	fmt.Println("fill opt:", time.Since(start))
	k := 1
	for j := 2; j < xBound-1; j++ {
		if opt[yBound-1][j] < opt[yBound-1][k] {
			k = j
		}
	}
	// return p, k
	seam := retraceSeam(p, k)
	// fmt.Println(seam)
	return seam
}

func retraceSeam(p [][]int, k int) []point {
	seam := make([]point, len(p))
	for i := len(p) - 1; i >= 0; i-- {
		seam[i] = point{k, i}
		k += p[i][k]
	}
	return seam
}

func main() {
	reader, err := os.Open("./Broadway_tower_edit.jpg")
	if err != nil {
		log.Fatal(err)
	}
	defer reader.Close()
	m, err := jpeg.Decode(reader)
	if err != nil {
		log.Fatal(err)
	}
	// bounds := m.Bounds()
	myImage := NewMyImg(m)
	newImageMatrix := getImage(m)
	energyMatrix := fillEnergyMatrixAsync(myImage)
	// p, opt := findSeamPaths(myImage, energyMatrix)
	xBound := len(newImageMatrix[0])
	yBound := len(newImageMatrix)
	for i := 0; i < 300; i++ {
		start := time.Now()
		seam := findBestSeam(xBound, yBound, energyMatrix)
		fmt.Println("seam:", time.Since(start))
		newImageMatrix = fillNewImageMatrix(newImageMatrix, seam)
		fmt.Println("image:", time.Since(start))
		energyMatrix = updateEnergyMatrix2(xBound, yBound, energyMatrix, seam)
		fmt.Println("energy:", time.Since(start))
		xBound = len(newImageMatrix[0])
		yBound = len(newImageMatrix)
		fmt.Println("total:", time.Since(start))
	}
	myImage = fillImage(xBound, yBound, newImageMatrix)
	writer, err := os.Create("./test2.jpg")
	if err != nil {
		log.Fatal(err)
	}
	jpeg.Encode(writer, myImage, nil)
}
