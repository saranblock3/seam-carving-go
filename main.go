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

func fillEnergyMatrix(m image.Image) [][]float64 {
	start := time.Now()
	energyMatrix := make([][]float64, m.Bounds().Max.Y)
	for y := 0; y < m.Bounds().Max.Y; y++ {
		energyRow := make([]float64, m.Bounds().Max.X)
		for x := 0; x < m.Bounds().Max.X; x++ {
			energyRow[x] = calculatePixelEnergy(m, x, y)
		}
		energyMatrix[y] = energyRow
	}
	fmt.Println("energy calculations:", time.Since(start))
	return energyMatrix
}

func fillEnergyMatrixAsync(m image.Image) [][]float64 {
	start := time.Now()
	energyMatrix := make([][]float64, m.Bounds().Max.Y)
	ch := make(chan energyRowEntry)
	for y := 0; y < m.Bounds().Max.Y; y++ {
		go func(y int, m image.Image, ch chan energyRowEntry) {
			energyRow := make([]float64, m.Bounds().Max.X)
			for x := 0; x < m.Bounds().Max.X; x++ {
				energyRow[x] = calculatePixelEnergy(m, x, y)
			}
			ch <- energyRowEntry{y, energyRow}
		}(y, m, ch)
	}
	for y := 0; y < m.Bounds().Max.Y; y++ {
		energyRow := <-ch
		energyMatrix[energyRow.y] = energyRow.energyRow
	}
	fmt.Println("energy calculations:", time.Since(start))
	// fmt.Println(energyMatrix)
	return energyMatrix
}

func updateEnergyMatrix(m image.Image, energyMatrix [][]float64, seam []point) [][]float64 {
	start := time.Now()
	newEnergyMatrix := make([][]float64, m.Bounds().Max.Y)
	for y := 0; y < m.Bounds().Max.Y; y++ {
		newRow := make([]float64, m.Bounds().Max.X)
		seamFound := false
		seamBehind := false
		for x := 0; x < len(energyMatrix[0]); x++ {
			if x == seam[y].x {
				seamFound = true
				continue
			}
			if seamFound {
				newRow[x-2] = calculatePixelEnergy(m, x-2, y)
				newRow[x-1] = calculatePixelEnergy(m, x, y)
				seamFound = false
				seamBehind = true
			} else if seamBehind {
				newRow[x-1] = energyMatrix[y][x]
			} else {
				newRow[x] = energyMatrix[y][x]
			}
		}
		newEnergyMatrix[y] = newRow
	}
	fmt.Println("energy calculations:", time.Since(start))
	return newEnergyMatrix
}

func fillEnergyMatrixAsync2(m image.Image) [][]float64 {
	start := time.Now()
	energyMatrix := make([][]float64, m.Bounds().Max.Y)
	ch := make(chan energyRowEntry)
	for y := 0; y < m.Bounds().Max.Y; y++ {
		go func(y int, m image.Image, ch chan energyRowEntry) {
			energyRow := make([]float64, m.Bounds().Max.X)
			for x := 0; x < m.Bounds().Max.X; x++ {
				energyRow[x] = calculatePixelEnergy(m, x, y)
			}
			ch <- energyRowEntry{y, energyRow}
		}(y, m, ch)
	}
	for y := 0; y < m.Bounds().Max.Y; y++ {
		energyRow := <-ch
		energyMatrix[energyRow.y] = energyRow.energyRow
	}
	fmt.Println("energy calculations:", time.Since(start))
	// fmt.Println(energyMatrix)
	return energyMatrix
}

type point struct {
	x int
	y int
}

type seamEntry struct {
	index int
	seam  []point
	score float64
}

func findSeams(energyMatrix [][]float64) [][]point {
	start := time.Now()
	seams := make([][]point, len(energyMatrix[0])-2)
	seamScores := make(map[int]float64, len(energyMatrix[0])-2)
	ch := make(chan seamEntry)
	for x := 1; x < len(energyMatrix[0])-1; x++ {
		go func(x int, energyMatrix [][]float64, ch chan seamEntry) {
			seam := make([]point, len(energyMatrix))
			seam[0] = point{x, 0}
			seamScore := 0.0
			xSeam := x
			for y := 1; y < len(energyMatrix); y++ {
				// fmt.Println(xSeam)
				if xSeam-1 < 1 {
					energyMock := (energyMatrix[y][xSeam] + energyMatrix[y][xSeam+1]) / 2
					seamScore += varianceEnergies(energyMock, energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
					xSeam += minEnergy(math.Inf(1), energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
				} else if xSeam+1 >= len(energyMatrix[0])-1 {
					energyMock := (energyMatrix[y][xSeam-1] + energyMatrix[y][xSeam]) / 2
					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMock)
					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], math.Inf(1))
				} else {
					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
				}
				seam[y] = point{xSeam, y}
			}
			ch <- seamEntry{x, seam, seamScore}
		}(x, energyMatrix, ch)
		// seams[x-1] = seam
	}
	for x := 1; x < len(energyMatrix[0])-1; x++ {
		currentSeamEntry := <-ch
		seams[currentSeamEntry.index-1] = currentSeamEntry.seam
		seamScores[currentSeamEntry.index-1] = currentSeamEntry.score
	}
	// minScore := math.Inf(1)
	// minSeamIdx := 0
	// for seamIdx, score := range seamScores {
	// 	if score < minScore {
	// 		minSeamIdx = seamIdx
	// 		minScore = score
	// 	}
	// }
	fmt.Println("score calulations:", time.Since(start))
	return seams
}

func scoreSeams(energyMatrix [][]float64) []point {
	start := time.Now()
	seams := make([][]point, len(energyMatrix[0])-2)
	seamScores := make(map[int]float64, len(energyMatrix[0])-2)
	ch := make(chan seamEntry)
	for x := 1; x < len(energyMatrix[0])-1; x++ {
		go func(x int, energyMatrix [][]float64, ch chan seamEntry) {
			seam := make([]point, len(energyMatrix))
			seam[0] = point{x, 0}
			seamScore := 0.0
			xSeam := x
			for y := 1; y < len(energyMatrix); y++ {
				// fmt.Println(xSeam)
				if xSeam-1 < 1 {
					energyMock := (energyMatrix[y][xSeam] + energyMatrix[y][xSeam+1]) / 2
					seamScore += varianceEnergies(energyMock, energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
					xSeam += minEnergy(math.Inf(1), energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
				} else if xSeam+1 >= len(energyMatrix[0])-1 {
					energyMock := (energyMatrix[y][xSeam-1] + energyMatrix[y][xSeam]) / 2
					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMock)
					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], math.Inf(1))
				} else {
					seamScore += varianceEnergies(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
					xSeam += minEnergy(energyMatrix[y][xSeam-1], energyMatrix[y][xSeam], energyMatrix[y][xSeam+1])
				}
				seam[y] = point{xSeam, y}
			}
			ch <- seamEntry{x, seam, seamScore}
		}(x, energyMatrix, ch)
		// seams[x-1] = seam
	}
	for x := 1; x < len(energyMatrix[0])-1; x++ {
		currentSeamEntry := <-ch
		seams[currentSeamEntry.index-1] = currentSeamEntry.seam
		seamScores[currentSeamEntry.index-1] = currentSeamEntry.score
	}
	minScore := math.Inf(1)
	minSeamIdx := 0
	for seamIdx, score := range seamScores {
		if score < minScore {
			minSeamIdx = seamIdx
			minScore = score
		}
	}
	fmt.Println("score calulations:", time.Since(start))
	return seams[minSeamIdx]
}

func varianceEnergies(energies ...float64) float64 {
	return meanSqaureEnergies(energies...) - math.Pow(meanEnergies(energies...), 2)
}

func meanEnergies(energies ...float64) float64 {
	sum := 0.0
	for _, v := range energies {
		sum += v
	}
	return sum / float64(len(energies))
}

func meanSqaureEnergies(energies ...float64) float64 {
	sum := 0.0
	for _, v := range energies {
		sum += math.Pow(v, 2)
	}
	return sum / float64(len(energies))
}

func minEnergy(energies ...float64) int {
	minEnergy := energies[0]
	minEnergyIdx := 0
	for i, v := range energies {
		if v < minEnergy {
			minEnergyIdx = i
		}
	}
	return minEnergyIdx - 1
}

func calculatePixelEnergy(m image.Image, x, y int) float64 {
	squareError := 0.0
	for j := -1; j < 2; j++ {
		if y+j < m.Bounds().Min.Y {
			continue
		}
		if y+j >= m.Bounds().Max.Y {
			continue
		}
		for i := -1; i < 2; i++ {
			if x+i < m.Bounds().Min.X {
				continue
			}
			if x+i >= m.Bounds().Max.X {
				continue
			}
			if i == 0 && j == 0 {
				continue
			}
			neighbourX := x + i
			neighbourY := y + i
			c1 := m.At(x, y)
			c2 := m.At(neighbourX, neighbourY)
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

func fillNewImage(m MyImage, seam []point) *MyImage {
	newImage := image.NewRGBA(image.Rect(0, 0, m.Bounds().Max.X-1, m.Bounds().Max.Y))
	for y := 0; y < m.Bounds().Max.Y; y++ {
		seamFound := false
		for x := 0; x < m.Bounds().Max.X; x++ {
			if seam[y].x == x && seam[y].y == y {
				seamFound = true
				continue
			}
			if seamFound {
				newImage.Set(x-1, y, m.At(x, y))
			} else {
				newImage.Set(x, y, m.At(x, y))
			}
		}
	}
	return &MyImage{newImage, map[image.Point]color.Color{}}
}

func main() {
	reader, err := os.Open("./955422-artwork-Simon-Stlenhag-science-fiction-futuristic.jpg")
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
	energyMatrix := fillEnergyMatrixAsync(myImage)
	for i := 0; i < 1000; i++ {
		seam := scoreSeams(energyMatrix)
		// fmt.Println(seam)
		// fmt.Println(i)
		newImage := fillNewImage(*myImage, seam)
		energyMatrix = updateEnergyMatrix(newImage, energyMatrix, seam)
		// for _, p := range seam {
		// 	color := color.RGBA{85, 165, 34, 255}
		// 	// r, g, b, a := m.At(90, y).RGBA()
		// 	myImage.Set(p.x, p.y, color)
		// }
		myImage = newImage
	}
	writer, err := os.Create("./test.jpg")
	if err != nil {
		log.Fatal(err)
	}
	jpeg.Encode(writer, myImage, nil)
}
