# Environment
- Operating System: Linux Mint 19
- Programming Language: Golang
- Compiler: go compiler version 1.11.4

# Arguments
```
Usage of ./main:
  -h 
        Show the help menu.
  -gap int
        The Score for the a gap. (default -4)
  -m int
        The score for a match. (default 5)
  -mis int
        The score for a mismatch. (default -3)
  -protein
        Flag that tells the program the given sequence is an amino acid sequence. 
        Default assumes nucleotide sequence.
  -s1 string
        First sequence to align or the path to a fasta file.
  -s2 string
        Second sequence to align or the path to a fasta file.
  -scoringMatrix string
        The scoring matrix used to grade mismatches. (default: PAM250)
        Options: PAM(250|30)
                 BLOSUM(62|45|80)
  -type string
        Type of alignment that will be performed. (global or local) (default "global")
  -cores int
        The number of cores to use. (default 12)
```

# How to build the program
> ```go run main.go <ARGS>```
### or
> ```go build main.go```
> <br>```./main <ARGS>```

# How to run the other versions of the program
If you would like to run the previous versions of this program you can by downloading the code from github and reverting it to a previous version. <br>
Note that only the final version of the program can preform alignments using the local type.

1. Perform a git clone
> ```git clone https://github.com/grif1179/bioFinalCode.git```

2. Perform one of the following reverts based on the previous version you would like to run.

### Using dynamic batch after looping through all elements
>```git checkout b3ae6e626619999b884640f6b0f33b0fa19c4381```

### Using set batch size
> ```git checkout f8e75eac455b3a9703e54adceaa8e3e4d56f5cca```

### Using a buffered channel
> ```git checkout 3727daf7489d3ed3d6e062ec14accef73b8313d3```

3. Return to head
> ```git checkout HEAD```