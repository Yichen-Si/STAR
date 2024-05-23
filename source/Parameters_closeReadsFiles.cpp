#include "Parameters.h"
#include "ErrorWarning.h"
#include <fstream>
#include <sys/stat.h>
void Parameters::closeReadsFiles() {
    for (uint imate=0; imate<readFilesIn.size(); imate++) {//open readIn files
        if ( inOut->readIn[imate].is_open() ) {
std::cerr << "Closing readIn file " << readFilesIn[imate] << std::endl;
            inOut->readIn[imate].close();
        }
        if (readFilesCommandPID[imate]>0) {
std::cerr << "Killing  CommandPID " << imate << std::endl;
            kill(readFilesCommandPID[imate],SIGKILL);
        }
    };
};
