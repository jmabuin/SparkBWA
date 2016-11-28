/*
 * SeqSparkParser.cpp
 *
 *  Created on: 25 nov, 2016
 *      Author: chema
 */

#include "SeqSparkParser.h"
#include "Utils.h"
#include <zlib.h>

const uint8_t SeqSparkParser::_codeTab[26] = { 0, 4, 1, 4, 4, 4, 2, //A -> G
		4, 4, 4, 4, 4, 4, 4, //H->N
		4, 4, 4, 4, 4, 3, 4, //O->U
		4, 4, 4, 4, 4 //V->Z
		};
const uint8_t SeqSparkParser::_decodeTab[5] = { 'A', 'C', 'G', 'T', 'N' };

const uint8_t SeqSparkParser::_complements[5] = { 3, 2, 1, 0, 4 };

SeqSparkParser::SeqSparkParser(Options* options, bool withLock, int type, size_t BUFFER_SIZE, char *sequencesRaw) {

	//char modes[5];
	/*colorspace or base space*/
	_colorspace = options->isColorSpace();
	_trimPrimer = options->trimPrimer();

	this->sequencesRaw = sequencesRaw;

	currentSequenceIndex = 0;

	/*create mutex*/
	_withLock = withLock;
	pthread_mutex_init(&_mutex, NULL);

	// For now only FASTQ format
	_format = FILE_FORMAT_FASTQ;
	_size = BUFFER_SIZE;
	_length = 0;

	_buffer = new uint8_t[_size + 1];
	if (_buffer == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
			__FUNCTION__, __LINE__);
	}

	Utils::log("FASTQ format in Spark parser\n");
}

SeqSparkParser::~SeqSparkParser() {

	/*destroy mutex*/
	pthread_mutex_destroy(&_mutex);

}

size_t SeqSparkParser::insertSequence(char *seq) {

	Sequence sequence = Sequence();

	size_t n = getFastqSeq(sequence);

	sequences.push_back(sequence);

	return n;

}

void SeqSparkParser::resizeBuffer(size_t nsize) {
	if (nsize <= _size) {
		return;
	}

	//allocate a new buffer
	_size = nsize * 2;
	uint8_t* nbuffer = new uint8_t[_size];
	if (!nbuffer) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	//copy the old data
	memcpy(nbuffer, _buffer, _length);

	//release the old buffer
	delete[] _buffer;
	_buffer = nbuffer;
}

size_t SeqSparkParser::getFastqSeq(Sequence& seq) {
	int ch;
	int nprimers;
	bool trimmable;

	char *currentSequence = this->sequencesRaw;//(char *)sequences[currentSequenceIndex].c_str();

	//find the header
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '@')
		;
	if (ch == -1)
		return 0; //reach the end of file

	//read the sequence name (only one line)
	_length = 0;
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (isspace(ch)) {
			ch = '\0';
		}
		_buffer[_length++] = ch;
	}
	if (ch == -1) {
		Utils::exit("Incomplete file\n");
	}
	_buffer[_length] = '\0';

	/*trim characters /[12]$ like BWA*/
	if (_length > 2 && _buffer[_length - 2] == '/'
			&& (_buffer[_length - 1] == '1' || _buffer[_length - 1] == '2')) {
		_length -= 2;
		_buffer[_length] = '\0';
	}

	//save the sequence name
	seq.setNameSize(_length + 1); /*adjust the name buffer size*/
	strcpy((char*) seq._name, (char*) _buffer);

	//read the sequence bases
	_length = 0;
	nprimers = 0;
	do {
		//filter out the blank lines
		while ((ch = myfgetc(currentSequence)) != -1 && (ch == '\r' || ch == '\n'))
			;
		if (ch == -1)
			Utils::exit("Incomplete FASTQ file\n");
		if (ch == '+')
			break; //the comment line

		//encode and save the base
		if (!_colorspace) {
			/*base space encoding*/
			if (ch >= 'A' && ch <= 'Z') {
				ch -= 'A';
			} else if (ch >= 'a' && ch <= 'z') {
				ch -= 'a';
			} else {
				Utils::exit("Unexpected character %c at line %d in file %s\n",
						ch, __LINE__, __FILE__);
			}
			ch = _codeTab[ch];
		} else {
			/*color space encoding: the leading nucleotide will be encoded in base-space*/
			if (ch >= 'A' && ch <= 'Z') {
				ch = _codeTab[ch - 'A'];
				++nprimers;
			} else if (ch >= 'a' && ch <= 'z') {
				ch = _codeTab[ch - 'a'];
				++nprimers;
			} else if (ch >= '0' && ch <= '3') {
				ch -= '0';
			} else {
				ch = UNKNOWN_BASE;
			}
		}
		//save the current encoded base
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = ch;

	} while (1);

	/*trim the primer base?*/
	trimmable = _trimPrimer && nprimers == 1 && nprimers < _length;

	/*save the sequence length as well as the bases*/
	seq._length = _length - (trimmable ? 1 : 0);
	seq.setSequenceSize(seq._length, true);
	memcpy(seq._bases, _buffer + (trimmable ? 1 : 0), seq._length);

	//read the comment line (only one line)
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '\n')
		;
	if (ch == -1)
		Utils::exit("Incomplete FASTQ file\n");

	//read the quality scores
	_length = 0;
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (ch >= 33 && ch <= 127) {
			_buffer[_length++] = ch;
		}

		if (_length > seq._length)
			break;
	}

	/*for base-space reads*/
	if (!_colorspace && seq._length != _length) {
		Utils::exit(
				"The number of bases is not equal to the number of quality scores\n");
	}
	if (_colorspace) {
		int32_t i, j;
		/*copy the base quality scores from the 3'end o the 5'end*/
		for(i = _length - 1, j = seq._length - 1; i >= 0 && j >= 0; --i, --j){
			seq._quals[j] = _buffer[i];
		}
		/*pseudo base quality score*/
		for(;j >= 0; --j){
			seq._quals[j] = '!';
		}
	}else{
		memcpy(seq._quals, _buffer, _length);
	}

	return seq._length;
}

size_t SeqSparkParser::getFastqSeq(Sequence& seq, char *rawSeq) {
	int ch;
	int nprimers;
	bool trimmable;

	char *currentSequence = rawSeq;

	//find the header
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '@')
		;
	if (ch == -1)
		return 0; //reach the end of file

	//read the sequence name (only one line)
	_length = 0;
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (isspace(ch)) {
			ch = '\0';
		}
		_buffer[_length++] = ch;
	}
	if (ch == -1) {
		Utils::exit("Incomplete file\n");
	}
	_buffer[_length] = '\0';

	/*trim characters /[12]$ like BWA*/
	if (_length > 2 && _buffer[_length - 2] == '/'
			&& (_buffer[_length - 1] == '1' || _buffer[_length - 1] == '2')) {
		_length -= 2;
		_buffer[_length] = '\0';
	}

	//save the sequence name
	seq.setNameSize(_length + 1); /*adjust the name buffer size*/
	strcpy((char*) seq._name, (char*) _buffer);

	//read the sequence bases
	_length = 0;
	nprimers = 0;
	do {
		//filter out the blank lines
		while ((ch = myfgetc(currentSequence)) != -1 && (ch == '\r' || ch == '\n'))
			;
		if (ch == -1)
			Utils::exit("Incomplete FASTQ file\n");
		if (ch == '+')
			break; //the comment line

		//encode and save the base
		if (!_colorspace) {
			/*base space encoding*/
			if (ch >= 'A' && ch <= 'Z') {
				ch -= 'A';
			} else if (ch >= 'a' && ch <= 'z') {
				ch -= 'a';
			} else {
				Utils::exit("Unexpected character %c at line %d in file %s\n",
						ch, __LINE__, __FILE__);
			}
			ch = _codeTab[ch];
		} else {
			/*color space encoding: the leading nucleotide will be encoded in base-space*/
			if (ch >= 'A' && ch <= 'Z') {
				ch = _codeTab[ch - 'A'];
				++nprimers;
			} else if (ch >= 'a' && ch <= 'z') {
				ch = _codeTab[ch - 'a'];
				++nprimers;
			} else if (ch >= '0' && ch <= '3') {
				ch -= '0';
			} else {
				ch = UNKNOWN_BASE;
			}
		}
		//save the current encoded base
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = ch;

	} while (1);

	/*trim the primer base?*/
	trimmable = _trimPrimer && nprimers == 1 && nprimers < _length;

	/*save the sequence length as well as the bases*/
	seq._length = _length - (trimmable ? 1 : 0);
	seq.setSequenceSize(seq._length, true);
	memcpy(seq._bases, _buffer + (trimmable ? 1 : 0), seq._length);

	//read the comment line (only one line)
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '\n')
		;
	if (ch == -1)
		Utils::exit("Incomplete FASTQ file\n");

	//read the quality scores
	_length = 0;
	while ((ch = myfgetc(currentSequence)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (ch >= 33 && ch <= 127) {
			_buffer[_length++] = ch;
		}

		if (_length > seq._length)
			break;
	}

	/*for base-space reads*/
	if (!_colorspace && seq._length != _length) {
		Utils::exit(
				"The number of bases is not equal to the number of quality scores\n");
	}
	if (_colorspace) {
		int32_t i, j;
		/*copy the base quality scores from the 3'end o the 5'end*/
		for(i = _length - 1, j = seq._length - 1; i >= 0 && j >= 0; --i, --j){
			seq._quals[j] = _buffer[i];
		}
		/*pseudo base quality score*/
		for(;j >= 0; --j){
			seq._quals[j] = '!';
		}
	}else{
		memcpy(seq._quals, _buffer, _length);
	}

	return seq._length;
}
