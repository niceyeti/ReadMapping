#include "ReadMapping.hpp"

ReadCollection::ReadCollection(const string& readsFile)
{
    _buildCollection(readsFile);
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

bool ReadMapping::MapReads(Sequence& input, const string& alphabet, const string& readsPath, const int minMatchLength)
{
    vector<int> matchIndices;
    bool result = false;

    if (fileExists(readsPath)) {
        SequenceAlignment* sequenceAligner = new SequenceAlignment();
        //Construct ST
        SuffixTree* suffixTree = new SuffixTree(input.seq, alphabet);
        //Prepare ST, building A array
        suffixTree->PrepareST();

        //Build the read collection; this is very fuzzy, since it will take so long (TODO: read files can be huge, so this could be done concurrently (eg. read-k reads, process them, read another k reads, process, etc)
        ReadCollection* reads = new ReadCollection(readsPath);

        //Map the reads: foreach read, get the 
        for (int i = 0; i < reads->ReadVector.size(); i++) {
            Read& thisRead = reads->ReadVector[i]; //TODO: ugh, fix the ReadCollection api when its i/o req's are clear
            if (reads->ReadVector[i].data.length() > minMatchLength) {
                //get the deepest node of the longest match(es)
                TreeNode* deepest = suffixTree->FindLoc(reads->ReadVector[i].data, minMatchLength);
                
                //start/endLeafIndices of deepest span the leaves representing sufficient matching strings; this iterates them and computes their local alignments
                for (int j = deepest->StartLeafIndex; j <= deepest->EndLeafIndex; j++) {
                    //run SmithWaterman on each extracted substring

                }
            }
        }
        
        result = true;
    }
    else {
        cout << "ERROR file containing reads not found, mapping aborted: " << readsPath << endl;
    }

    return result;
}


