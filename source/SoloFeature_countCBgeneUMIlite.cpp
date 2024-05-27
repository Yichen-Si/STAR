#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "systemFunctions.h"

void SoloFeature::countCBgeneUMIlite()
{
    time_t rawTime;

    rguStride=2;
    rGeneUMI = new uint32[rguStride*nReadsMapped]; //big array for all CBs - each element is gene and UMI
    rCBp = new uint32*[nCB+1]; // start position of each barcode in the storage
    uint32 **rCBpa = new uint32*[nCB+1];

    rCBp[0]=rGeneUMI;
    rCBpa[0]=rGeneUMI;
    for (uint32 ii=0; ii<nCB; ii++) { // some WL barcode may have zero read?
        rCBp[ii+1] = rCBp[ii] + rguStride*readFeatSum->cbReadCount[ii];
        rCBpa[ii+1]= rCBp[ii+1];
    };

    //read and store the CB/gene/UMI from files
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished allocating arrays for Solo " << nReadsMapped*rguStride*4.0/1024/1024/1024 <<" GiB" <<endl;

    ///////////////////////////////////////////////////////////////////////////
    ////////////// Input records
    readFlagCounts.flagCounts.reserve(nCB*3/2);
    nReadPerCB.resize(nCB);       // total reads
    nReadPerCBunique.resize(nCB); // uniquely mapped reads
    // temprorary lazy version that does not count multi-mapped reads
    for (int ii=0; ii<P.runThreadN; ii++) {//TODO: this can be parallelized
        readFeatAll[ii]->inputRecords(rCBpa, rguStride, readFlagCounts, nReadPerCBunique, nReadPerCB, cbWL);
        readFeatSum->addStats(*readFeatAll[ii]);//sum stats: has to be done after inputRecords, since the stats values are updated there
    }
    nReadPerCBmax=0;
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        nReadPerCBmax=max(nReadPerCBmax,nReadPerCB[iCB]);
    };

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax<<endl;
    // *P.inOut->logStdOut << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax<<endl;

    ///////////////////////////////////////////////////////////////////////////
    /////////////////////////// collapse each CB
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);

                     //dedup options        //gene ID
    countMatStride = pSolo.umiDedup.yes.N + 1;
    countCellGeneUMI.resize(nReadsMapped*countMatStride/5+16); //matrix
    countCellGeneUMIindex.resize(nCB+1, 0); //Starting position in the count matrix corresponding to each CB
    collapseUMIall();

    P.inOut->logMain << "RAM for solo feature "<< SoloFeatureTypes::Names[featureType] <<"\n"
                     <<  linuxProcMemory() << flush;
    delete[] rGeneUMI;
    delete[] rCBp;
    delete[] rCBpa;

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
    if (P.debug)
        *P.inOut->logStdOut << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
};
