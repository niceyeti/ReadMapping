#include "ReadMapping.hpp"

ReadCollection::ReadCollection(const string& readsFile)
{
    _buildCollection(readsFile);
}

ReadCollection::~ReadCollection()
{}

void ReadCollection::_buildCollection(const string& readsFile)
{
    int seqnum = 0;
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
                    ReadVector.back().hitBegin = -1;
                    ReadVector.back().hitEnd = -1;
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

/*
Outputs the read mappings as described in the prg3Spec.
This function assumes any Read with negative map indices was a miss.
*/
void ReadCollection::Write(const string& outputPath)
{
    ofstream outputFile;

    if (fileExists(outputPath)) {
        outputFile.open(outputPath);
        if (outputFile.is_open()) {
            //output all the mapping results
            for (int i = 0; i < ReadVector.size(); i++) {
                //output formatted as: <Read_name> <Start index of hit> <End index of hit>.
                if (ReadVector[i].hitBegin >= 0 && ReadVector[i].hitEnd >= 0) {
                    outputFile << ReadVector[i].desc << " " << ReadVector[i].hitBegin << " " << ReadVector[i].hitEnd << "\n";
                }
                else {
                    outputFile << ReadVector[i].desc << "No hit found.\n";
                }
            }
        }
        else {
            cout << "ERROR output file could not be opened: " << outputPath << endl;
        }
    }
    else {
        cout << "ERROR output file for reads not found" << endl;
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

/*

*/
bool ReadMapping::MapReads(Sequence& input, const string& alphabet, const string& readsPath, const string& paramsPath, const int minMatchLength, const string& resultPath)
{
    float pctIdCoverage, pctLenCoverage;
    float maxLenCoverage;
    const float pctIdentityCoverage = 0.80, pctLengthCoverage = 0.90; //these params defined in prg3Spec
    Read maxRead;
    int maxReadIndex;
    Parameters params;
    vector<int> matchIndices;
    bool result = false;

    if (!fileExists(readsPath)) {
        cout << "ERROR file containing reads not found, mapping aborted: " << readsPath << endl;
        return false;
    }
    if (!fileExists(paramsPath)) {
        cout << "ERROR file containing alignment params not found, mapping aborted: " << paramsPath << endl;
        return false;
    }

    //get the sequence alignment params
    SequenceAlignment::ParseParamsFile(paramsPath, params);
    //make the tempAlignment
    Alignment alignment;
    //Make the aligner
    SequenceAlignment* sequenceAligner = new SequenceAlignment();
    //Construct ST
    SuffixTree* suffixTree = new SuffixTree(input.seq, alphabet);
    //Prepare ST, building A array
    suffixTree->PrepareST();

    //Build the read collection; this is very fuzzy, since it will take so long (TODO: read files can be huge, so this could be done concurrently (eg. read-k reads, process them, read another k reads, process, etc)
    ReadCollection* reads = new ReadCollection(readsPath);

    //Map the reads
    for (int i = 0; i < reads->ReadVector.size(); i++) {
        Read& curRead = reads->ReadVector[i]; //TODO: ugh, fix the ReadCollection api when its i/o req's are clear

        if (curRead.data.length() > minMatchLength) {
            //get the deepest node of the longest match(es)
            string& readStr = curRead.data;
            TreeNode* deepest = suffixTree->FindLoc(readStr, minMatchLength);

            cout << "Computing maximal match over " << (deepest->EndLeafIndex - deepest->StartLeafIndex) << " candidates" << endl;
            //start/endLeafIndices of deepest span the leaves representing sufficient matching strings; this iterates them and computes their local alignments
            maxLenCoverage = 0;
            for (int j = deepest->StartLeafIndex; j <= deepest->EndLeafIndex; j++) {
                //run SmithWaterman on each extracted substring
                //TODO: rather than generate new substrings using substr() could instead rewrite SmithWaterman to accept string index boundaries and a string ref
                //extract local string around match from genome input (big string)
                string genomeSubstring = input.seq.substr(std::max(0, j - (int)readStr.length()), 2 * (int)readStr.length());
                sequenceAligner->SmithWaterman(genomeSubstring, readStr, params, alignment, false);

                //get the coverage threshold vals
                pctIdCoverage = (float)alignment.matches / (float)(alignment.Length());
                pctLenCoverage = (float)alignment.Length() / (float)readStr.length();
                //firstly check if the alignment is of 'good' overall quality as defined by our prg3Spec thresholds
                if (pctIdCoverage >= pctIdentityCoverage && pctLenCoverage >= pctLengthCoverage) {
                    //track max read via its pctLenCoverage value
                    if (pctLenCoverage > maxLenCoverage) {
                        maxLenCoverage = pctLenCoverage;
                        curRead.hitBegin = j;
                        curRead.hitEnd = j + readStr.length();
                    }
                }
            }
        }
    }

    //output the reads to file
    reads->Write(resultPath);

    return result;
}


