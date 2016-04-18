#include "ReadMapping.hpp"

void usage()
{
    cout << "Usage:\n<read_mapper executable> <FASTA file containing reference genome sequence G> <FASTA file containing reads> <input alphabet file>" << endl;
}


//<read_mapper executable> <FASTA file containing reference genome sequence G> <FASTA file containing reads> <input alphabet file>
int main(int argc, char** argv)
{
    int minMatchLength;
    string inputFile, alphaFile, readFile;
    Sequence inputSequence;
    string alphabet;
    string input;
    SuffixTree st;
    Sequence seq;

    //For api for the hw
    if (argc != 4) {
        cout << "ERROR incorrect number of arguments." << endl;
        usage();
        return -1;
    }

    inputFile = argv[1];
    readFile = argv[2];
    alphaFile = argv[3];
    //minMatchLength = argv[4]

    inputFile= "Peach_Reference.fasta";
    readFile = "Peach_Reads.txt";
    alphaFile= "alphabet.txt";
    minMatchLength = 25;

    if (fileExists(inputFile)) {
        if (fileExists(readFile)) {
            if (fileExists(alphaFile)) {
                if (parseAlphabetFile(alphaFile, alphabet)) {
                    if (parseFastaFile(inputFile, inputSequence, alphabet)) {
                        ReadMapping* rm = new ReadMapping();
                        rm->MapReads(inputSequence, alphabet, readFile, minMatchLength);


                    }
                    else {
                        cout << "ERROR could not parse FASTA file" << endl;
                    }
                }
                else {
                    cout << "ERROR could not parse alphabet file" << endl;
                }
            }
            else {
                cout << "ERROR alphabet file not found: " << alphaFile << endl;
            }
        }
        else {
            cout << "ERROR file containing reads not found: " << readFile << endl;
        }
    }
    else {
        cout << "ERROR input file not found: " << inputFile << endl;
    }

    return 0;
}
