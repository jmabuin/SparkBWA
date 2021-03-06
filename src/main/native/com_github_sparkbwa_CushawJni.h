/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class com_github_sparkbwa_CushawJni */

#ifndef _Included_com_github_sparkbwa_CushawJni
#define _Included_com_github_sparkbwa_CushawJni
#ifdef __cplusplus
extern "C" {
#endif

// Utility preprocessor directive so only one change needed if Java class name changes
#define JNIFUNCTION_CUSHAW(sig) Java_com_github_sparkbwa_CushawJni_##sig

//JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(parseSequence1(JNIEnv * env, jobject object, jstring seq1));
//JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(parseSequence2(JNIEnv * env, jobject object, jstring seq2));
JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(CushawInit(JNIEnv * env, jobject object, jobjectArray stringArray));
JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(executeEstimateJNI(JNIEnv * env, jobject object, jstring seq1, jstring seq2));
JNIEXPORT jint JNICALL JNIFUNCTION_CUSHAW(loadIndexJNI(JNIEnv * env, jobject object, jstring indexPath));
JNIEXPORT jstring JNICALL JNIFUNCTION_CUSHAW(alignJNI(JNIEnv * env, jobject object, jstring seq1, jstring seq2));

#ifdef __cplusplus
}
#endif


#endif
