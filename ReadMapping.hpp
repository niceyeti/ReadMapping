#pragma once

#include "../../StringAlignment/StringAlignment/Header.hpp"
#include "../../mccreight//mccreight/SuffixTree.hpp"
#include "../../mccreight/mccreight/Util.hpp"


#ifndef _IOSTREAM_
#include <iostream>
#endif

#ifndef _STRING_
#include <string>
#endif

#ifndef _VECTOR_
#include <vector>
#endif

#ifndef _LIST_
#include <list>
#endif

#ifndef _IOMANIP_
#include <iomanip>
#endif

#ifndef _CHRONO_
#include <chrono>
#endif

#ifndef _CTIME_
#include <ctime>
#endif

//This mirrors the primitive read structure of the read test input files:
// > some description
// atcgtcgagtgcg...[the read sequence]
struct Read {
    string desc;
    string data;
};

//Primitive source-container for a bunch of reads, to be mapped to some sequence
class ReadCollection
{
public:
    ReadCollection(const string& readsFile);
    ~ReadCollection();
    vector<Read> ReadVector;
private:
    void _buildCollection(const string& readsFile);

};



class ReadMapping {
public:
//ReadMapping(inputSequence, alphabet);
    ReadMapping();
    ~ReadMapping();
    bool RunMapping(Sequence& input, const string& alphabet, const string& outputPath);

private:
    string _alphabet;
    Sequence* _inputSequence;
    SuffixTree* _suffixTree;
    SequenceComparer* _sequenceComparer;
    vector<int> _findLoci();
    bool _outputReadMap(const string& outputPath);

};



