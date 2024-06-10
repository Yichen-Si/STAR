#include "ReadAlign.h"
#include "readLoad.h"
#include "readBarcodeLoad.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"

int ReadAlign::oneRead() {//process one read: load, map, write

    //load read name, sequence, quality from the streams into internal arrays
    int readStatus[P.readNends];

    for (uint32 im=0; im<P.readNends; im++) { // load all mates for one read
        readStatus[im] = readLoad(*(readInStream[im]), P, readLength[im], readLengthOriginal[im], readNameMates[im], Read0[im], Read1[im], Qual0[im], clipMates[im], iReadAll, readFilesIndex, readFilter, readNameExtra[im]);
        if (readStatus[im] != readStatus[0]) {//check if the end of file was reached or not for all files
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: read files are not consistent, reached the end of the one before the other one\n";
            errOut << "SOLUTION: Check you your input files: they may be corrupted\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
    };

    if (readStatus[0]==-1) {//finished with the stream
        return -1;
    };
    readFileType=readStatus[0];

    if (P.outFilterBySJoutStage != 2) {
        for (uint32 im=0; im<P.readNmates; im++) {//not readNends: the barcode quality will be calculated separately
            for (uint64 ix=clipMates[im][0].clippedN; ix<readLengthOriginal[im]-clipMates[im][1].clippedN; ix++) {
                qualHist[im][(uint8)Qual0[im][ix]]++;
            };
        };
    };

    // Processed the part to be aligned
    if (P.readNmates==2) {//combine two mates together
        Lread=readLength[0]+readLength[1]+1;
        readLengthPairOriginal=readLengthOriginal[0]+readLengthOriginal[1]+1;
        if (Lread>DEF_readSeqLengthMax) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: Lread of the pair = " << Lread << "   while DEF_readSeqLengthMax=" << DEF_readSeqLengthMax <<endl;
            errOut << "Read Name="<<readNameMates[0]<<endl;
            errOut << "SOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile STAR"<<endl<<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        //marker for spacer base
        Read1[0][readLength[0]]=MARK_FRAG_SPACER_BASE;

        //copy 2nd mate into Read1[0] & reverse-complement
        complementSeqNumbers(Read1[1],Read1[0]+readLength[0]+1,readLength[1]);//complement. Here Read1[1] is still the 2nd mate's numeric-sequence. Later Read1[1] will be reverse complement of the combined read.
        for (uint ii=0;ii<readLength[1]/2;ii++) {
            swap(Read1[0][Lread-ii-1],Read1[0][ii+readLength[0]+1]); //reverse
        };

    } else {//1 mate

        if (readStatus[0]==-1) {//finished with the stream
            return -1;
        };

        Lread=readLength[0];
        readLengthPairOriginal=readLengthOriginal[0];
        readLength[1]=0;

    };
    complementSeqNumbers(Read1[0],Read1[1],Lread); //returns complement of Reads[ii]
    for (uint ii=0;ii<Lread;ii++) {//reverse
        Read1[2][Lread-ii-1]=Read1[1][ii];
    };

    // Process identifiers (for solo counting)
    // Could be more general // 2024UM
    if (seqScope) {
        cbMatch = -1;
        outputCB = 0;
        outputAnno = 0;
        // spatial barcode
        if ((int32) readLengthOriginal[1] < P.cbS+P.cbL) {
            cbMatch = -11;
            return 2; // Wrong length, should not happen
        }
        char cb[P.cbL+1];
        int32_t x, y;
        cbMatch = cbWL->query(Read0[1]+P.cbS, cb, x, y);
        if (cbMatch < 0 || cbMatch > 1) {
            if (P.skipAlignUnmatchWL) {
                soloRead->readBar->addStats(cbMatch);
                return 1; // Does not match white list, skip align
            }
        }
        if (cbMatch == 0 && !P.skipCBifExact) {
            cbCorrected = std::string(Read0[1]+P.cbS, P.cbL);
            outputCB = 1;
        } else if (cbMatch > 0) {
            cbCorrected = std::string(cb, P.cbL);
            outputCB = 1;
        }
        if (cbMatch >= 0) {
            cbInfo =  std::to_string(x) + "," + std::to_string(y);
            outputAnno = 1;
            sb = (uint64) x << 32 | y;
        }
        // UMI
        if ((int32) readLengthOriginal[1] < P.ubS+P.ubL) {
            cbMatch = -11;
            return 2; // Wrong length, should not happen
        }
        std::vector<uint8_t> nonACGTs;
        umint4 = seq2bits2(Read0[1]+P.ubS, P.ubL, nonACGTs);
        umiAmbig = nonACGTs.size(); // We could allow ambiguous bases in UMI, but STARsolo does not
        umiHomopoly = 0;
        if (umiAmbig > 0) {
            cbMatch = -23;
        } else {
            for (uint i = 0; i < 4; ++i) {
                if (umint4 == homopolymer[i]) {
                    umiHomopoly = 1;
                    cbMatch = -24;
                    break;
                }
            }
        }
        soloRead->readBar->addStats(cbMatch);
    }

    statsRA.readN++;
    statsRA.readBases += readLength[0]+readLength[1];

    //max number of mismatches allowed for this read
    outFilterMismatchNmaxTotal=min(P.outFilterMismatchNmax, (uint) (P.outFilterMismatchNoverReadLmax*(readLength[0]+readLength[1])));

    //map the read
    if (P.pGe.gType==101) {//SpliceGraph
        mapOneReadSpliceGraph();
    } else {//all other cases - standard alignment algorithm
        mapOneRead();
    };

    peOverlapMergeMap();

    multMapSelect();

    mappedFilter(); // assign unmapType and add to statsRA

    transformGenome();//for now genome transformation happens after multimapper selection, and mapping filter

    if (!peOv.yes) {//if the alignment was not mates merged - otherwise the chimeric detection was already done
        chimericDetection();
    };

    if (P.pCh.out.bam && chimRecord) {//chimeric alignment was recorded in main BAM files, and it contains the representative portion, so non-chimeric aligmnent is not output
        return 0;
    };

    waspMap();

    #ifdef OFF_BEFORE_OUTPUT
        #warning OFF_BEFORE_OUTPUT
        return 0;
    #endif

    //write out alignments
    outputAlignments();

    {
    #ifdef DEBUG_OutputLastRead
        lastReadStream.seekp(ios::beg);
        lastReadStream << iReadAll <<" "<< readName <<endl;
    #endif
    };

    return 0;

};
