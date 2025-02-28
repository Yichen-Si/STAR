#ifndef CODE_ReadAlignChunk
#define CODE_ReadAlignChunk

#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "OutSJ.h"
#include "Transcriptome.h"
#include "BAMoutput.h"
#include "Quantifications.h"

class ReadAlignChunk {//chunk of reads and alignments
public:
    Parameters& P;
    ReadAlign* RA;
    std::shared_ptr<SbcdWL> cbWL;

    Transcriptome *chunkTr;

    char **chunkIn; //space for the chunk of input reads
    array<uint64, MAX_N_MATES> chunkInSizeBytesTotal;

    char *chunkOutBAM, *chunkOutBAM1;//space for the chunk of output SAM
    OutSJ *chunkOutSJ, *chunkOutSJ1;

    BAMoutput *chunkOutBAMcoord, *chunkOutBAMunsorted, *chunkOutBAMquant;
    Quantifications *chunkQuants;

    istringstream** readInStream;
    ostringstream*  chunkOutBAMstream;
    ofstream chunkOutBAMfile;
    string chunkOutBAMfileName;

    bool noReadsLeft;
    uint iChunkIn; //current chunk # as read from .fastq
    uint iChunkOutSAM; //current chunk # writtedn to Aligned.out.sam
    int iThread; //current thread
    uint chunkOutBAMtotal; //total number of bytes in the write buffer

    ReadAlignChunk(Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk, std::shared_ptr<SbcdWL> _cbWL = nullptr);
    ~ReadAlignChunk();
    void processChunks();
    void mapChunk();
    void chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut);
    void chunkFstreamCat (fstream &chunkOut, ofstream &allOut, bool mutexFlag, pthread_mutex_t &mutexVal);
    void chunkFilesCat(ostream *allOut, string filePrefix, uint &iC);

    Genome &mapGen;
private:
};
#endif
