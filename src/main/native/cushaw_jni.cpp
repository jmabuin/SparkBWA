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

#include "com_github_sparkbwa_CushawJni.h"
#include "cushaw3-v3.0.3/core/Sequence.h"



Sequence sequence1;
Sequence sequence2;

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(parseSequence1(JNIEnv * env, jobject object, jstring seq1)) {

    const char *buf;
    //const jbyte *str;
    buf = (*env).GetStringUTFChars(seq1, NULL);

    if (buf == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    printf("%s parsing sequence: %s",__func__, buf);

    sequence1 = Sequence(buf, FILE_FORMAT_FASTQ);

    (*env).ReleaseStringUTFChars(seq1, buf);

    return sequence1._length;

}

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(parseSequence2(JNIEnv * env, jobject object, jstring seq2)) {

    const char *buf;
    //const jbyte *str;
    buf = (*env).GetStringUTFChars(seq2, NULL);

    if (buf == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    printf("%s parsing sequence: %s",__func__, buf);

    sequence2 = Sequence(buf, FILE_FORMAT_FASTQ);

    (*env).ReleaseStringUTFChars(seq2, buf);

    return sequence2._length;
}

JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(CushawInit(JNIEnv * env, jobject object, jint argN, jobjectArray arguments, jintArray argumentsSizes)) {

    //Parte argumentos
   	char **argv;
   	char **argvTmp;

   	int stringCount = (*env).GetArrayLength(arguments);//env->GetArrayLength(stringArray);

   	argvTmp = (char **) malloc(stringCount*sizeof(char **));

   	return 1;

}