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

#define BIG_BUFFER_SIZE 25*1024*1024

PairedEndSpark *pairedSpark;
Genome *genome;
Options *options;
SAMSpark *samSpark;
char *outputSAM;

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(CushawInit(JNIEnv * env, jobject object, jobjectArray stringArray)) {

    //Parte argumentos
   	//char **argv;
   	char **argvTmp;

   	int myArgc = (*env).GetArrayLength(stringArray);//env->GetArrayLength(stringArray);

   	argvTmp = (char **) malloc(myArgc*sizeof(char *));

    //outputSAM = NULL;
    outputSAM = (char *) malloc (sizeof(char) * BIG_BUFFER_SIZE);

    fprintf(stderr,"[%s] Creating Cushaw C++ Objects ...\n", __func__);
    fprintf(stderr,"[%s] Creating Options ...\n", __func__);

    //char **myArgv;

    int i = 0;

    for (i=0; i<myArgc; i++) {
        jstring string = (jstring) (*env).GetObjectArrayElement(stringArray, i);
        argvTmp[i] = (char *)(*env).GetStringUTFChars(string, 0);


    }


    options = new BSOptions();
    options->parse(myArgc, argvTmp);

    fprintf(stderr,"[%s] Creating Genome ...\n", __func__);
    genome = new Genome(options, true);

    fprintf(stderr,"[%s] Creating SAM ...\n", __func__);
    samSpark = new SAMSpark(options, genome, outputSAM);

    fprintf(stderr,"[%s] Creating PairedEnd ...\n", __func__);
    pairedSpark = new PairedEndSpark(options, genome, samSpark);

    fprintf(stderr,"[%s] Objects created ...\n", __func__);

   	return 1;

}

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(executeEstimateJNI(JNIEnv * env, jobject object, jstring seq1, jstring seq2)) {

    const char *buf1, *buf2;
    //const jbyte *str;

    fprintf(stderr, "[%s] JMAbuin Inside executeEstimate\n",__func__);


    buf1 = (*env).GetStringUTFChars(seq1, NULL);

    if (buf1 == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    char *charSeq1 = (char *) malloc(sizeof(char)*(strlen(buf1)+1));
    strcpy(charSeq1, buf1);

    fprintf(stderr, "[%s] JMAbuin length of first split %d\n", __func__, strlen(charSeq1));

    buf2 = (*env).GetStringUTFChars(seq2, NULL);

    if (buf2 == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    char *charSeq2 = (char *) malloc(sizeof(char)*(strlen(buf2)+1));
    strcpy(charSeq2, buf2);

    fprintf(stderr, "[%s] JMAbuin length of second split %d\n", __func__, strlen(charSeq2));

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