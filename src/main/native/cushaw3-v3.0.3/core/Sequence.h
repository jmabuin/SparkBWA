/*
 * Sequence.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#include "Macros.h"

struct Sequence
{
	/*member functions*/
	Sequence();
	Sequence(const Sequence& s);
	Sequence(const char *newSeq, int fileFormat);
	~Sequence();
	void setNameSize(size_t size);
	void setSequenceSize(size_t size, bool quals);
	inline void releaseQuals() {
		if (_quals) {
			delete[] _quals;
			_quals = NULL;
		}
	}
	void clear();
	void print(FILE* file);

	size_t parseFastqSeq(char *seq);
	void resizeBuffer(size_t nsize);

	size_t BUFFER_SIZE = 4095;

	/*member variables*/
	uint8_t* _name;
	uint8_t* _bases;
	uint8_t* _rbases;
	uint8_t* _quals;
	uint32_t _length;
	uint32_t _tlength;
	uint32_t _nameSize;
	uint32_t _seqSize;

	uint8_t* _buffer;
	size_t _size;
	bool _trimPrimer;

	static const uint8_t _codeTab[26];
};

#endif /* SEQUENCE_H_ */
