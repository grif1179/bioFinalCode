package main

import (
	"fmt"
	// "math"
	"sync"
	"flag"
	"strings"
	"os"
	"bufio"
	"log"
	"errors"
	"strconv"
)

// Creating type directions whose base type is an integer
type directions int

// enum for each direction the table could go.
const (
	up directions = iota
	left
	diagonal
	nothing
)

// A struct that represent each element in the
//matrix.
type score struct {
	value int
	direction directions
	topBase string
	leftBase string
}

/**
* Variables that store the values given
* by the user using the go flags library.
*/
var (
	s1 = flag.String("s1", "", "First sequence to align or the path to a fasta file.")
	s2 = flag.String("s2", "", "Second sequence to align or the path to a fasta file.")
	matchScore = flag.Int("m", 5, "The score for a match.")
	misMatchScore = flag.Int("mis", -3, "The score for a mismatch.")
	gap = flag.Int("gap", -4, "The Score for the a gap.")
	scoringMatrix = flag.String("scoringMatrix", "", "The scoring matrix used to grade mismatches. (default: PAM250)\nOptions: PAM(250|30)\n\t BLOSUM(62|45|80)")
	alignType = flag.String("type", "global", "Type of alignment that will be performed. (global or local)")
	isProtein = flag.Bool("protein", false, "Flag that tells the program the given sequence is a amino acid sequence. \nDefault assumes nucleotide sequence.")
	cores = flag.Int("cores", 12, "The number of cores to use.")
)


/**
* Creates a 2d matrix that uses the score struct for each element.
* x: number of Rows.
* y: number of Cols.
* returns the resulting matrix.
*/
func makeMatrix(x, y int) [][]score {
	matrix := make([][]score, y)
	for i := range matrix {
		row := make([]score, x)
		matrix[i] = row
	}
	return matrix
}

/**
* Initializes a given matrix for a global alignment.
* It initializes it by making the first row and column 
* contain the gap value multiplied by the current index.
*/
func initGlobalMatrix(matrix [][]score, gap int, s1 string, s2 string) {
	rows := len(matrix) + 1
	cols := len(matrix[0]) + 1
	matrix[0][0].value = 0

	for i := 1; i < cols-1; i++ {
		matrix[0][i].value = i * gap
		matrix[0][i].direction = left
		matrix[0][i].topBase = string(s1[i-1])
	}

	for j := 1; j < rows-1; j++ {
		matrix[j][0].value = j * gap
		matrix[j][0].direction = up
		matrix[j][0].leftBase = string(s2[j-1])
	}	
}

/**
* Initializes a given matrix for a local alignment.
* It initializes it by making the first row and column 
* contain zeroes.
*/
func initLocalMatrix(matrix [][]score, s1 string, s2 string) {
	rows := len(matrix) + 1
	cols := len(matrix[0]) + 1
	matrix[0][0].value = 0

	for i := 1; i < cols-1; i++ {
		matrix[0][i].value = 0
		matrix[0][i].direction = left
		matrix[0][i].topBase = string(s1[i-1])
	}

	for j := 1; j < rows-1; j++ {
		matrix[j][0].value = 0
		matrix[j][0].direction = up
		matrix[j][0].leftBase = string(s2[j-1])
	}	
}

/**
* Prints the score matrix for debugging.
*/
func printMatrix(matrix [][]score) {
	rows := len(matrix)
	cols := len(matrix[0])
	sep := " "
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			fmt.Printf("%v%s", matrix[i][j], sep)
		}
		fmt.Println()
	}
	fmt.Println()
}


/**
* Basic check to see whether the argument given is an
* actual sequence or the path to a fasta file containing the sequence.
*/
func isFasta(str string) (bool, error) {
	if strings.Contains(str, ".fa") || strings.Contains(str, ".fasta") {
		// file exists
		_, err := os.Stat(str);
		if  err == nil {
			return true, nil
		} 
		return false, err
	}
	return false, nil
}

/*
* Extract the sequence in a fasta file as a string.
* filename: The path to a fasta file.
* container: A string pointer that will contain the sequence.
*/
func extractFasta(filename string, container *string) {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var temp string
	for scanner.Scan() {
		if !strings.Contains(scanner.Text(), ">") {
			temp += scanner.Text()
		}
	}
	*container = temp
}

/*
* Turns a string slice into an integer slice.
*/
func sliceAtoi(strSlice []string) []int {
	tempSlice := make([]int, len(strSlice))
	for i, value := range strSlice {
		tempSlice[i], _ = strconv.Atoi(value)
	}
	return tempSlice
}

/**
* Based on the scoring matrix picked the corresponding 
* scoring matrix will be extracted.
* scoreType: Scoring matrix to use. (PAM250|PAM30|BLOSUM62|BLOSUM45|BLOSUM80)
* returns: A 2-D integer slice containing the scoring matrix,
*		   a map containing the index of each character in the matrix, and
*          an error message if anything goes wrong.
*/
func readScoringMatrix(scoreType string) ([][]int, map[string]int, error) {
	filename := "scoring_matrices/"
	switch scoreType {
	case "PAM250": 
		filename += "PAM250.csv"
	case "PAM30":
		filename += "PAM30.csv"
	case "BLOSUM62":
		filename += "BLOSUM62.csv"
	case "BLOSUM45":
		filename += "BLOSUM45.csv"
	case "BLOSUM80":
		filename += "BLOSUM80.csv"
	default:
		return nil, nil, errors.New("Invalid score matrix selected")
	}

	var referenceSlice []string
	misMatchMap := make(map[string]int)
	scoreMatrix := make([][]int, 25)
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var i int
	for scanner.Scan() {
		if i == 0 {
			referenceSlice = strings.Fields(scanner.Text())
			for i, v := range referenceSlice {
				misMatchMap[v] = i
			}
		} else {
			temp := strings.Fields(scanner.Text())
			scoreMatrix[i-1] = sliceAtoi(temp)
		}
		i++
	}
	return scoreMatrix, misMatchMap, nil
}

var jobsWg sync.WaitGroup
var wg sync.WaitGroup

type job struct {
	i int
	j int
}

func process(ch chan job, matrix [][]score, misMatchMap map[string]int, scoringMatrix [][]int, myID int) {
	match := *matchScore
	mis := *misMatchScore
	gap := *gap
	s1 := *s1
	s2 := *s2
	var upVal, leftVal, diaVal int
	for  {
		job, ok := <- ch
		if !ok {
			wg.Done()
			return
		} 

		i := job.i
		j := job.j
		upVal = matrix[i-1][j].value + gap
		leftVal = matrix[i][j-1].value + gap
		if misMatchMap != nil {
			firstCharLoc := misMatchMap[strings.ToUpper(string(s1[j-1]))]
			secondCharLoc := misMatchMap[strings.ToUpper(string(s2[i-1]))]
			matchScore := scoringMatrix[firstCharLoc][secondCharLoc]
			diaVal = matrix[i-1][j-1].value + matchScore
		} else {
			if s1[j-1] == s2[i-1] {
				diaVal = matrix[i-1][j-1].value + match
			} else {
				diaVal = matrix[i-1][j-1].value + mis
			}
		}

		if diaVal > upVal && leftVal < diaVal {
			matrix[i][j].value = diaVal
			matrix[i][j].direction = diagonal
			matrix[i][j].leftBase = string(s2[i-1])
			matrix[i][j].topBase = string(s1[j-1])

		} else if upVal > leftVal {
			matrix[i][j].value = upVal
			matrix[i][j].direction = up
			matrix[i][j].leftBase = string(s2[i-1])
			matrix[i][j].topBase = "-"
		} else {
			matrix[i][j].value = leftVal
			matrix[i][j].direction = left
			matrix[i][j].leftBase = "-"
			matrix[i][j].topBase = string(s1[j-1])
		}

		// jobsWg.Done()
	}
}

func diaGlobalAlignment(matrix [][]score, match, mis, gap int, s1, s2 string, scoringMatrix [][]int, misMatchMap map[string]int) {
	// Creating communication channels
	rows := len(matrix)
	cols := len(matrix[0])
	cores := *cores
	myChannel := make(chan job, 6)
	// myChannels := make([]chan []job, cores)
	// for i := 0; i < cores; i++ {
	// 	myChannels[i] = make(chan []job)
	// 	wg.Add(1)
	// 	go process(myChannels[i], matrix, misMatchMap, scoringMatrix, i)
	// }
	// batchSize := 12

	// Creating and assigning jobs for
	//diagnols that start from the first row.
	for j := 1; j < cols; j++ {
		// jobs := make([]job, 0)
		tempj := j
		// currCore := 0

		myChannel := make(chan job, 6)
		for i := 0; i < cores; i++ {
			wg.Add(1)
			go process(myChannel, matrix, misMatchMap, scoringMatrix, i)
		}
		for i := 1; i < rows; i++ {
			// jobsWg.Add(1)
			if tempj <= 1 {
				myChannel <- job{i, tempj}
				break
			} 

			myChannel <- job{i, tempj}
			tempj--
		}
		close(myChannel)

		wg.Wait()	
		// jobsWg.Wait()
	}
	// Creating and assigning jobs for
	//diagnols that start from the last col.
	for i := 1; i < rows; i++ {
		tempi := i
		tempj := cols - 1
		myChannel := make(chan job, 6)
		for i := 0; i < cores; i++ {
			wg.Add(1)
			go process(myChannel, matrix, misMatchMap, scoringMatrix, i)
		}
		for {
			// jobsWg.Add(1)
			if tempj == 1 || tempi == rows-1 {
				myChannel <- job{tempi, tempj}
				break
			} 
			myChannel <- job{tempi, tempj}

			tempi++
			tempj--		
		}
		close(myChannel)

		wg.Wait()	
		// jobsWg.Wait()
	}

	close(myChannel)

	wg.Wait()	
	
	// fmt.Println()
	// printMatrix(matrix)

	var topSeq, leftSeq, matcher []string
	currElement := matrix[rows-1][cols-1]
	currX := rows - 1
	currY := cols - 1
	totalScore := currElement.value

	for currElement.topBase != "" && currElement.leftBase != "" {
		topSeq = append([]string{currElement.topBase}, topSeq...)
		leftSeq = append([]string{currElement.leftBase}, leftSeq...)
		if currElement.direction == up {
			currX--	
			matcher = append([]string{" "}, matcher...)
		} else if currElement.direction == left {
			currY--
			matcher = append([]string{" "}, matcher...)
		} else {
			currX--
			currY--
			if currElement.leftBase == currElement.topBase {
				matcher = append([]string{"|"}, matcher...)
			} else {
				matcher = append([]string{":"}, matcher...)
			}
		}
		currElement = matrix[currX][currY]
	}

	// fmt.Println("Final Matrix:")
	// printMatrix(matrix)

	fmt.Printf("The Total Score is %d\n", totalScore)
	seqPerLine := 50
	for i := 0; i < len(topSeq); i += seqPerLine {
		if i+seqPerLine >= len(topSeq) {
			fmt.Printf("s1: %s\n", strings.Join(topSeq[i:], " "))
			fmt.Printf("    %s\n", strings.Join(matcher[i:], " "))
			fmt.Printf("s2: %s\n", strings.Join(leftSeq[i:], " "))
		} else {
			fmt.Printf("s1: %s\n", strings.Join(topSeq[i:i+seqPerLine], " "))
			fmt.Printf("    %s\n", strings.Join(matcher[i:i+seqPerLine], " "))
			fmt.Printf("s2: %s\n", strings.Join(leftSeq[i:i+seqPerLine], " "))

			fmt.Printf("\n")
		}
	}
}


func main() {
	// lengthX := 5
	// lengthY := 5
	// cores := 12

	flag.Parse() // Get user arguments

	// Is sequence 1 a fasta file or literal sequence
	if isS1Fasta, err := isFasta(*s1); err != nil {
		fmt.Fprintf(os.Stderr, "error: %v\n", err)
	} else if isS1Fasta == true {
		extractFasta(*s1, s1)
	}

	// Is sequence 2 a fasta file or literal sequence
	if isS2Fasta, err := isFasta(*s2); err != nil {
		fmt.Fprintf(os.Stderr, "error: %v\n", err)
	} else if isS2Fasta == true {
		extractFasta(*s2, s2)
	}

	var scoreMatrix [][]int
	var misMatchMap map[string]int
	var err error

	// Retrieve scoring matrix 
	if *scoringMatrix != "" {
		scoreMatrix, misMatchMap, err = readScoringMatrix(*scoringMatrix)
		if err != nil {
			fmt.Fprintf(os.Stderr, "error: %v\n", err)
		}
	} else if *isProtein {
		// Protein sequence given with no scoring matrix so default 
		//matrix is PAM250.
		fmt.Println("Using PAM250")
		scoreMatrix, misMatchMap, err = readScoringMatrix("PAM250")
	}

	matrix := makeMatrix(len(*s1)+1, len(*s2)+1)

	// Is alignment global or local
	if strings.ToLower(*alignType) == "global" {
		initGlobalMatrix(matrix, *gap, *s1, *s2)
		// printMatrix(matrix)
		diaGlobalAlignment(matrix, *matchScore, *misMatchScore, *gap, *s1, *s2, scoreMatrix, misMatchMap)
		// printMatrix(matrix)
	} else if strings.ToLower(*alignType) == "local" {
		initLocalMatrix(matrix, *s1, *s2)
		// localAlignment(matrix, *matchScore, *misMatchScore, *gap, *s1, *s2, scoreMatrix, misMatchMap)
	} else {
		fmt.Fprintf(os.Stderr, "error: Invalid alignment type selected.\n")
		os.Exit(1)
	}
}