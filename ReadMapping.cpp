#include "ReadMapping.hpp"

ReadCollection::ReadCollection(const string& readsFile)
{
}

ReadCollection::~ReadCollection()
{}

void ReadCollection::_buildCollection(const string& readsFile)
{
    int seqnum;
    ifstream inputFile;
    string line;
    int progress, i;
    size_t fsize, bytesRead, thousandsRead;

    if (ReadVector.size() > 0) {
        ReadVector.clear();
    }

    if (fileExists(readsFile)) {
        inputFile.open(readsFile);
        if (inputFile.is_open()) {
            //init
            progress = i = 0;
            bytesRead = thousandsRead = 0;
            fsize = getFilesize(readsFile);
            //reserve a large amount of mem, to avoid frequent reallocs
            ReadVector.reserve(fsize / 1000); //read files are huge, so div by some constant factor less than the avg size of a read (desc+read) in bytes.

            cout << "Parsing reads from file: " << readsFile << ". This may require significant time, for read files > 200kb." << endl;
            while (getline(inputFile, line)) {
                //if line begins with '>' grab desription line
                if (seqnum == 0 && line.length() > 0 && line[0] == '>') {
                    ReadVector.resize(ReadVector.size()+1);
                    ReadVector.back().desc = line;
                }
                //else grab the read data itself
                else {
                    ReadVector.back().data = line;
                }

                //report progress
                bytesRead += line.length();
                if (bytesRead % 10000 > thousandsRead) {
                    thousandsRead = bytesRead % 1000;
                    cout << "\rProgress: " << (int)((float(bytesRead) / float(fsize)) * 100) << "% complete          " << endl;
                }
            }
            inputFile.close();
        }
        else {
            cout << "ERROR file could not be opened! " << readsFile << endl;
        }
    }
    else {
        cout << "ERROR, file not found: " << readsFile << endl;
    }
}

ReadMapping::ReadMapping()
{}

ReadMapping::~ReadMapping()
{}

bool ReadMapping::_outputReadMap(const string& outputPath)
{
    ofstream outputFile(outputPath, ofstream::trunc);
    bool result = false;

    if (fileExists(outputPath)) {
        if (outputFile.good()) {
            //iterate over the read mappings



            outputFile.close();
            result = true;
        }
        else {
            cout << "Failed to open output file: " << outputPath;
        }
    }
    else {
        cout << "ERROR output file not found, mapping aborted: " << outputPath << endl;
    }

    return result;
}

bool ReadMapping::RunMapping(Sequence& input, const string& alphabet, const string& outputPath)
{
    bool result = false;

    if (fileExists(outputPath)) {
        _suffixTree = new SuffixTree(input.seq, alphabet);
        _sequenceComparer = new SequenceComparer();





        result = true;
    }
    else {
        cout << "ERROR output file not found, mapping aborted: " << outputPath << endl;
    }

    return result;
}


