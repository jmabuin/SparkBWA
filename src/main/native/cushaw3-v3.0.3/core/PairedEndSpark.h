/*
 * PairedEndSpark.h
 *
 *  Created on: 25 nov, 2016
 *      Author: chema
 */

#ifndef CORE_PAIREDENDSPARK_H_
#define CORE_PAIREDENDSPARK_H_

#include "Macros.h"
#include "Options.h"
#include "MemEngine.h"
#include "Thread.h"
#include "SeqSparkParser.h"
#include "SAMSpark.h"

class PairedEndSpark {
public:
	PairedEndSpark(Options* options, Genome * genome, SAMSpark* sam);
	virtual ~PairedEndSpark();

	/*execute the paired-end alignment*/
	void executeEstimate(char *seq1, char *seq2);
	void estimateInsertSize(int minAlignedPairs, int maxReadBatchSize, int numReadBatchs, SeqSparkParser *parser1, SeqSparkParser *parser2);

	char *align(char *seq1, char*seq2);

	inline void lock() {
		if (_numThreads > 1) {
			pthread_mutex_lock(&_mutex);
		}
	}

	inline void unlock() {
		if (_numThreads > 1) {
			pthread_mutex_unlock(&_mutex);
		}
	}

	MemEngine *engine;

private:

    void imprimir();
	/*private member variables*/
	Options* _options;
	Genome* _genome;
	SAMSpark* _sam;

	int _numThreads;
	/*thread parameters*/
	vector<ThreadParams*> _threads;

	pthread_mutex_t _mutex;

	static void* _threadFunc(void*);
};

#endif /* CORE_PAIREDENDSPARK_H_ */
