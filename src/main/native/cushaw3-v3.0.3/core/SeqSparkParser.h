/*
 * SeqSparkParser.h
 *
 *  Created on: 25 nov, 2016
 *      Author: chema
 */

#ifndef CORE_SEQSPARKPARSER_H_
#define CORE_SEQSPARKPARSER_H_

#include <vector>
#include "Macros.h"
#include "Sequence.h"
#include "Utils.h"
#include "MyFile.h"
#include "Options.h"

class SeqSparkParser {
public:
	SeqSparkParser(Options* options, bool withLock, int type, size_t BUFFER_SIZE = 4095, char *sequencesRaw = NULL);
	virtual ~SeqSparkParser();

	size_t insertSequence(char *seq);
	size_t getFastqSeq(Sequence& seq, char *rawSeq);

	inline size_t getSeqLockFree(Sequence& seq) {
			size_t ret;
			int32_t numNs = 0, i;

			/*read the sequence from the file*/
			/*

			if (_format == FILE_FORMAT_FASTA) {
				ret = getFastaSeq(seq);
			} else if (_format == FILE_FORMAT_FASTQ) {
				ret = getFastqSeq(seq);
			} else { //BAM/SAM format
				ret = getBSamSeq(seq);
			}
			*/

			if(_format == FILE_FORMAT_FASTQ) {
				ret = getFastqSeq(seq);
			}
			else {
				ret = -1;
			}

			/*compute the reverse complement*/
			if (ret > 0) {
				/*trim all consecutive Ns at the 3 prime end for each read*/
				for (i = seq._length - 1; i >= 0; --i){
					if(seq._bases[i] != UNKNOWN_BASE){
						break;
					}
				}
				seq._length = max(1, i + 1);

				/*compute the number of Ns in the full length*/
				for (i = seq._length - 1; i >= 0; --i) {
					if (seq._bases[i] == UNKNOWN_BASE) {
						++numNs;
					}
				}
				seq._tlength = seq._length - numNs;

				//compute the reverse complement of the sequence
				reverseComp(seq._rbases, seq._bases, seq._length, _colorspace);
			}

			return ret;
		}

		inline size_t getSeqLockFree(Sequence* seqs, int maxSeqs) {
			int index;
			for (index = 0; index < maxSeqs; ++index) {
				if (!getSeqLockFree(seqs[index])) {
					break;
				}
			}

			return index;
		}

private:

	/*buffered file operations*/
	inline int myfgetc(char *character) {
		int m = *character;
		character++;

		return m;
	}

	inline void reverseComp(uint8_t* rbases, uint8_t* bases, size_t length,
				bool colorspace) {
			size_t off;
			size_t halfLength = length / 2;

			if (colorspace) {
				for (size_t i = 0; i < halfLength; i++) {
					off = length - i - 1;
					rbases[off] = bases[i];
					rbases[i] = bases[off];
				}
				if (length & 1) {
					rbases[halfLength] = bases[halfLength];
				}
			} else {
				for (size_t i = 0; i < halfLength; i++) {
					off = length - i - 1;
					rbases[off] = _complements[bases[i]];
					rbases[i] = _complements[bases[off]];
				}
				if (length & 1) {
					rbases[halfLength] = _complements[bases[halfLength]];
				}
			}
		}
		inline void reverseComp(uint8_t* bases, size_t length, bool colorspace) {
			uint8_t ch;
			size_t off;
			size_t halfLength = length / 2;

			if (colorspace) {
				for (size_t i = 0; i < halfLength; i++) {
					off = length - i - 1;
					ch = bases[off];
					bases[off] = bases[i];
					bases[i] = ch;
				}
			} else {
				for (size_t i = 0; i < halfLength; i++) {
					off = length - i - 1;
					ch = bases[off];
					bases[off] = _complements[bases[i]];
					bases[i] = _complements[ch];
				}
				if (length & 1) {
					bases[halfLength] = _complements[bases[halfLength]];
				}
			}
		}


	/*private member functions*/
	void resizeBuffer(size_t nsize);
	size_t getFastqSeq(Sequence& seq);

	std::vector<Sequence> sequences;


	int currentSequenceIndex;

	/*private member variables*/
	//buffer for file reading
	uint8_t* _buffer;
	size_t _length;
	size_t _size;


	int _format;
	bool _colorspace;
	bool _trimPrimer;

	/* Sequences raw string*/
	char *sequencesRaw;

	/*interal lock for concurrent accesses to a single file*/
	pthread_mutex_t _mutex;
	bool _withLock;

	static const uint8_t _codeTab[26];
	static const uint8_t _decodeTab[5];
	static const uint8_t _complements[5];
};

#endif /* CORE_SEQSPARKPARSER_H_ */
