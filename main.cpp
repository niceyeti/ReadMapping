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
    string resultFile;
    string paramPath = "parameters.config"; //TODO: not part of the cli api for prg3, so hardcoded for now
    SuffixTree st;
    Sequence seq;

    /*
    //api for the hw
    if (argc != 4) {
        cout << "ERROR incorrect number of arguments." << endl;
        usage();
        return -1;
    }
    

    inputFile = argv[1];
    readFile = argv[2];
    alphaFile = argv[3];
    minMatchLength = 25; //hard coded per assignment spec
    */

    inputFile = "Peach_Reference.fasta";
    readFile = "Peach_Reads.txt";
    alphaFile= "alphabet.txt";
    minMatchLength = 15;

    //result location is generated based on input file name, where the input should be named as 'Peach_Reference.txt'
    size_t underscore = inputFile.find_first_of('_');
    if (underscore != inputFile.npos) {
        resultFile = inputFile.substr(0, underscore) + "_Results.txt";
    }
    else {
        resultFile = inputFile + "_Results.txt";
    }

    if (fileExists(inputFile)) {
        if (fileExists(readFile)) {
            if (fileExists(alphaFile)) {
                if (parseAlphabetFile(alphaFile, alphabet)) {
                    if (parseFastaFile(inputFile, inputSequence, alphabet)) {
                        ReadMapping* rm = new ReadMapping();
                        rm->MapReads(inputSequence, alphabet, readFile, paramPath, minMatchLength, resultFile);
                        cout << "Mapping complete." << flush;
                        cin >> readFile;
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
        cout << "ERROR input fasta file not found: " << inputFile << endl;
    }

    return 0;
}
