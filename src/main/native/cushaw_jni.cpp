/**
  * Copyright 2016 José Manuel Abuín Mosquera <josemanuel.abuin@usc.es>
  *
  * This file is part of SparkBWA.
  *
  * SparkBWA is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * SparkBWA is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  * GNU General Public License for more details.
  *
  * You should have received a copy of the GNU General Public License
  * along with SparkBWA. If not, see <http://www.gnu.org/licenses/>.
  */

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <string.h>

#include "com_github_sparkbwa_CushawJni.h"
#include "cushaw3-v3.0.3/core/Sequence.h"
#include "cushaw3-v3.0.3/core/Options.h"
#include "cushaw3-v3.0.3/core/BSOptions.h"
#include "cushaw3-v3.0.3/core/Genome.h"
#include "cushaw3-v3.0.3/core/SAMSpark.h"
#include "cushaw3-v3.0.3/core/PairedEndSpark.h"


PairedEndSpark *pairedSpark;
Genome *genome;
Options *options;
SAMSpark *samSpark;
char *outputSAM;

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(CushawInit(JNIEnv * env, jobject object)) {

    //Parte argumentos
   	//char **argv;
   	//char **argvTmp;

   	//int stringCount = (*env).GetArrayLength(arguments);//env->GetArrayLength(stringArray);

   	//argvTmp = (char **) malloc(stringCount*sizeof(char **));

    outputSAM = NULL;

    options = new BSOptions();
    genome = new Genome(options, true);
    samSpark = new SAMSpark(options, genome, outputSAM);
    pairedSpark = new PairedEndSpark(options, genome, samSpark);

   	return 1;

}

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(executeEstimateJNI(JNIEnv * env, jobject object, jstring seq1, jstring seq2)) {

    const char *buf1, *buf2;
    //const jbyte *str;
    buf1 = (*env).GetStringUTFChars(seq1, NULL);

    if (buf1 == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    char *charSeq1 = (char *) malloc(sizeof(char)*(strlen(buf1)+1));
    strcpy(charSeq1, buf1);


    buf2 = (*env).GetStringUTFChars(seq2, NULL);

    if (buf2 == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    char *charSeq2 = (char *) malloc(sizeof(char)*(strlen(buf2)+1));
    strcpy(charSeq2, buf2);

    pairedSpark->executeEstimate(charSeq1, charSeq2);

    free(charSeq1);
    free(charSeq2);

    (*env).ReleaseStringUTFChars(seq1, buf1);
    (*env).ReleaseStringUTFChars(seq2, buf2);

    return 1;
}

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(loadIndexJNI(JNIEnv * env, jobject object, jstring indexPath)) {


    std::string bwtFileName, saFileName, annFileName;
    std::string pacFileName, basePacFileName;
    std::string baseBitmapFileName;

    const char *buf1;
    //const jbyte *str;
    buf1 = (*env).GetStringUTFChars(indexPath, NULL);

    if (buf1 == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    char *charIndex = (char *) malloc(sizeof(char)*(strlen(buf1)+1));
    strcpy(charIndex, buf1);

    options->setBwtFileBase(charIndex);

    std::string& bwtBase = options->getBwtFileBase();
    /*using the reverse orientation of the genome*/
    bwtFileName = bwtBase + ".rbwt";
    saFileName = bwtBase + ".rsa";
    annFileName = bwtBase + ".ann";
    baseBitmapFileName = bwtBase + ".map";
    pacFileName = bwtBase + ".pac";
    basePacFileName = bwtBase + ".nt.pac"; /*for color-space*/

    /*read the data from files*/
    genome->_init(bwtFileName, saFileName, annFileName, baseBitmapFileName,
    	pacFileName, basePacFileName, options->isColorSpace(),
    	options->maskAmbiguous());

    return 1;

}

JNIEXPORT jstring JNICALL JNIFUNCTION_CUSHAW(alignJNI(JNIEnv * env, jobject object, jstring seq1, jstring seq2)) {

    const char *buf1, *buf2;
    //const jbyte *str;
    buf1 = (*env).GetStringUTFChars(seq1, NULL);

    if (buf1 == NULL) {
        return (*env).NewStringUTF(""); /* OutOfMemoryError already thrown */
    }

    char *charSeq1 = (char *) malloc(sizeof(char)*(strlen(buf1)+1));
    strcpy(charSeq1, buf1);


    buf2 = (*env).GetStringUTFChars(seq2, NULL);

    if (buf2 == NULL) {
        return (*env).NewStringUTF(""); /* OutOfMemoryError already thrown */
    }

    char *charSeq2 = (char *) malloc(sizeof(char)*(strlen(buf2)+1));
    strcpy(charSeq2, buf2);

    char *result = (char *)malloc(sizeof(char) * 8192);

    strcpy(result, pairedSpark->align(charSeq1, charSeq2));

    free(charSeq1);
    free(charSeq2);

    (*env).ReleaseStringUTFChars(seq1, buf1);
    (*env).ReleaseStringUTFChars(seq2, buf2);


    return (*env).NewStringUTF(result);

}