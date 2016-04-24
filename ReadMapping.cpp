#include "ReadMapping.hpp"

ReadCollection::ReadCollection(const string& readsFile)
{
    _buildCollection(readsFile);
}

ReadCollection::~ReadCollection()
{}

void ReadCollection::_buildCollection(const string& readsFile)
{
    bool getData = false;
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
                if (line.length() > 2) {
                    //if line begins with '>' grab desription line of this read
                    if (line[0] == '>') {
                        ReadVector.resize(ReadVector.size() + 1);
                        ReadVector.back().desc = line;
                        ReadVector.back().hitBegin = -1;
                        ReadVector.back().hitEnd = -1;
                        getData = true;
                    }
                    //else grab the read data itself
                    else if(getData){
                        ReadVector.back().data = line;
                        getData = false;
                    }

                    if (ReadVector.size() % 10000 == 9999) {
                        cout << "\r" << ReadVector.size() << " reads read.      " << flush;
                    }
                }
            }
            ReadVector.shrink_to_fit();
            cout << "\rProgress: 100% complete             " << endl;
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

void ReadCollection::Clear()
{
    ReadVector.clear();
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
                    outputFile << ReadVector[i].desc << "  " << ReadVector[i].hitBegin << " " << ReadVector[i].hitEnd << "\n";
                }
                else {
                    outputFile << ReadVector[i].desc << "  No hit found.\n";
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
    int locBegin, numChars, rootMisses, similarityMisses, hitCount, numAlignments;
    int stBuildtime, stPreptime, mapTime, readsBuildtime, outputTime;
    float pctIdCoverage, pctLenCoverage;
    float maxLenCoverage;
    const float minIdentityCoverage = (float)0.90; //these params defined in prg3Spec: 0.8 0.9
    const float minLengthCoverage = (float)0.80;
    Read maxRead;
    Parameters params;
    Alignment alignment;
    vector<int> matchIndices;
    bool result = false;
    clock_t c_start;
    string dummy;

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
    //Make the aligner
    SequenceAlignment* sequenceAligner = new SequenceAlignment();
    //Construct ST
    c_start = clock();
    SuffixTree* suffixTree = new SuffixTree(input.seq, alphabet);
    stBuildtime = clock() - c_start;

    //Prepare ST, building A array
    c_start = clock();
    suffixTree->PrepareST();
    stPreptime = clock() - c_start;

    //Build the read collection; this is very fuzzy, since it will take so long (TODO: read files can be huge, so this could be done concurrently (eg. read-k reads, process them, read another k reads, process, etc)
    c_start = clock();
    ReadCollection* reads = new ReadCollection(readsPath);
    readsBuildtime = clock() - c_start;

    //Map the reads
    c_start = clock();
    rootMisses = similarityMisses = hitCount = numAlignments = 0;
    for (int i = 0; i < reads->ReadVector.size(); i++) {
        Read& curRead = reads->ReadVector[i]; //TODO: ugh, fix the ReadCollection api when its i/o req's are clear

        if (curRead.data.length() > minMatchLength) {
            //get the deepest node of the longest match(es)
            string& readStr = curRead.data;
            TreeNode* deepest = suffixTree->FindLoc(readStr, minMatchLength);

            if (deepest == deepest->Parent) {
                //most similar node is root; this is a miss, or else every suffix in the genome will be aligned with the read!
                //cout << "ERROR hit root" << endl;
                //cin >> del;
                rootMisses++;
            }
            else {
                //cout << "Start: " << deepest->StartLeafIndex << "  End: " << deepest->EndLeafIndex << ":  ";
                //for (int k = deepest->StartLeafIndex; k <= deepest->EndLeafIndex; k++) {
                //    cout << suffixTree->A[k] << ", ";
                //}
                //cout << endl;
                //suffixTree->PrintPrefix(deepest);
                //TODO: if deepest node found is root, skip the read; otherwise the read will be aligned with every suffix of the genome string!

                //cout << "Computing maximal match over " << (deepest->EndLeafIndex - deepest->StartLeafIndex + 1) << " candidate(s)" << " for read: " << curRead.desc << endl;
                //start/endLeafIndices of deepest span the leaves representing sufficient matching strings; this iterates them and computes their local alignments
                maxLenCoverage = 0;
                for (int j = deepest->StartLeafIndex; j <= deepest->EndLeafIndex; j++) {
                    //get absolute string index from A array, minus some padding, making sure we don't go below zero
                    locBegin = max(0, suffixTree->A[j] - (int)readStr.length());
                    //get the number of chars to extract (2*read.length), making sure we don't exceed remainder of string
                    numChars = min(2 * (int)readStr.length(), (int)input.seq.length() - locBegin);

                    //cout << j << " < " << deepest->EndLeafIndex << endl;
                    //only run SmithWaterman alignment if more than one candidate location
                    //if (deepest->StartLeafIndex < deepest->EndLeafIndex) {
                        //TODO: rather than generate new substrings using substr() could instead rewrite SmithWaterman to accept string index boundaries and a string ref
                        //extract local string around match from genome input (big string)
                    string genomeSubstring = input.seq.substr(locBegin, numChars);
                    //cout << locBegin << "(" << numChars << "):  " << genomeSubstring << endl;
                    sequenceAligner->SmithWaterman(genomeSubstring, readStr, params, alignment, false);

                    //track the number of alignments
                    numAlignments++;

                    //get the coverage threshold vals
                    pctIdCoverage = (float)alignment.matches / (float)(alignment.Length());
                    pctLenCoverage = (float)alignment.Length() / (float)readStr.length();
                    //cout << "id coverage: " << pctIdCoverage << "  len coverage: " << pctLenCoverage << endl;
                    //firstly check if the alignment is of 'good' overall quality as defined by our prg3Spec thresholds
                    if (pctIdCoverage >= minIdentityCoverage && pctLenCoverage >= minLengthCoverage) {
                        //track max read via its pctLenCoverage value
                        if (pctLenCoverage > maxLenCoverage) {
                            maxLenCoverage = pctLenCoverage;
                            curRead.hitBegin = locBegin;
                            curRead.hitEnd = locBegin + numChars;
                        }
                    }
                    else {
                        similarityMisses++;
                        //cout << "miss, insufficient similarity:  " << genomeSubstring << "\n                         " << readStr << endl;
                        //cout << j << " miss, insufficient similarity: " << pctIdCoverage << " id coverage,  " << pctLenCoverage << " len coverage   " << alignment.Length() << " alignment length" << endl;
                    }
                    //}
                    //only one possible alignment found, so just store it as the max
                    //else {
                    //    curRead.hitBegin = locBegin;
                    //    curRead.hitEnd = locBegin + numChars;
                    //}

                    //dbg 
                    if (curRead.hitBegin < -1) {
                        cout << "ERROR hitBegin preceded beginning of string: " << curRead.hitBegin << " for read: " << curRead.desc << endl;
                    }
                    if (curRead.hitEnd >= (int)input.seq.length()) {
                        cout << "ERROR hitEnd off end of string: " << curRead.hitEnd << " for read: " << curRead.desc << endl;
                    }
                } //end inner for

                //track hits by whether or not the read was assigned a non-negative index
                if (curRead.hitBegin >= 0) {
                    hitCount++;
                }
            }
        }

        if ((i & 0x00000FFF) == 0) {
            cout << "\rMapping " << (100.0 * (float)i / (float)reads->ReadVector.size()) << "% complete               " << flush;
        }
    } // end outer for
    mapTime = clock() - c_start;
    
    cout << "\rMapping 100% complete.    " << endl;

    //output the reads to file
    c_start = clock();
    reads->Write(resultPath);
    outputTime = clock() - c_start;

    //Report miss/hit stats
    cout << "Mapped reads: " << hitCount << endl;
    cout << "Number of alignments: " << numAlignments << endl;
    cout << "Min-match length misses: " << rootMisses << endl;
    cout << "Similarity misses: " << similarityMisses;
    if (numAlignments > 0) {
        cout << "    Miss rate (similarity misses / num alignments): " << ((float)similarityMisses / (float)numAlignments);
    }
    cout << endl;
    cout << "Hit ratio: " << ((float)hitCount / (float)reads->ReadVector.size()) << endl;
    cout << "Alignments per total num reads: " << ((float)numAlignments / (float)reads->ReadVector.size()) << endl;
    //Report timing stats
    cout << "Timing stats: " << endl;
    cout << "  SuffixTree Build time : " << ((1000.0 * (float)stBuildtime) / CLOCKS_PER_SEC) << "ms" << endl;
    cout << "  SuffixTree PrepareST() time: " << ((1000.0 * (float)stPreptime) / CLOCKS_PER_SEC) << "ms" << endl;
    cout << "  Read collection build time: " << ((1000.0 * (float)readsBuildtime) / CLOCKS_PER_SEC) << "ms" << endl;
    cout << "  MapReads() time: " << ((1000.0 * (float)mapTime) / CLOCKS_PER_SEC) << "ms" << endl;
    cout << "  Output time: " << ((1000.0 * (float)outputTime) / CLOCKS_PER_SEC) << "ms" << endl;

    cin >> dummy;

    //destroy the read collection
    reads->Clear();
    //destroy the tree
    suffixTree->Clear();

    return result;
}


