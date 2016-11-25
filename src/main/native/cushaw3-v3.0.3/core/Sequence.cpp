/*
 * Sequence.cpp
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#include "Sequence.h"
#include "Utils.h"

const uint8_t Sequence::_codeTab[26] = { 0, 4, 1, 4, 4, 4, 2, //A -> G
		4, 4, 4, 4, 4, 4, 4, //H->N
		4, 4, 4, 4, 4, 3, 4, //O->U
		4, 4, 4, 4, 4 //V->Z
		};

Sequence::Sequence() {
	_name = NULL;
	_bases = NULL;
	_rbases = NULL;
	_quals = NULL;
	_nameSize = _seqSize = 0;
	_length = 0;
	_tlength = 0;
	_size = 0;
	_buffer = NULL;
	_trimPrimer = false;
}

Sequence::Sequence(const Sequence & s) {
	_length = s._length;
	_tlength = s._tlength;
	_nameSize = s._nameSize;
	_seqSize = s._seqSize;
	_size = 0;
	_buffer = NULL;
	_trimPrimer = false;

	if (_length == 0) {
		_name = NULL;
		_bases = NULL;
		_rbases = NULL;
		_quals = NULL;
		_nameSize = 0;
		_seqSize = 0;
		return;
	}
	if (s._name) {
		_name = new uint8_t[_nameSize];
		if (_name == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		strcpy((char*) _name, (const char*) s._name);
	}
	if (s._bases) {
		_bases = new uint8_t[_seqSize];
		if (_bases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_bases, s._bases, _length);
	}
	if (s._rbases) {
		_rbases = new uint8_t[_seqSize];
		if (_rbases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_rbases, s._rbases, _length);
	}
	if (s._quals) {
		_quals = new uint8_t[_seqSize];
		if (_quals == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_quals, s._quals, _length);
	}
}

Sequence::Sequence(const char *newSeq, int fileFormat) {

    _name = NULL;
    _bases = NULL;
    _rbases = NULL;
    _quals = NULL;
    _nameSize = _seqSize = 0;
    _length = 0;
    _tlength = 0;
    _trimPrimer = false;

    if( fileFormat == FILE_FORMAT_FASTQ) {

    	this->_size = this->BUFFER_SIZE;
    	_buffer = new uint8_t[_size + 1];
    	if (_buffer == NULL) {
    		Utils::exit("Memory allocation failed in file %s in line %d\n",
    				__FUNCTION__, __LINE__);
    	}

    	// Call to parse fastq
    }

}


Sequence::~Sequence() {
	clear();
}
void Sequence::clear() {
	if (_name) {
		delete[] _name;
	}
	if (_bases) {
		delete[] _bases;
	}
	if (_rbases) {
		delete[] _rbases;
	}
	if (_quals) {
		delete[] _quals;
	}

	_name = NULL;
	_bases = NULL;
	_rbases = NULL;
	_quals = NULL;
	_length = 0;
	_tlength = 0;
	_nameSize = 0;
	_seqSize = 0;
}
void Sequence::setNameSize(size_t size) {
	if (size >= _nameSize) {
		_nameSize = size * 2;
		if (_name) {
			delete[] _name;
		}
		_name = new uint8_t[_nameSize];
		if (_name == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
}
void Sequence::setSequenceSize(size_t size, bool quals) {
	if (size >= _seqSize) {
		_seqSize = size * 2;
		/*forward strand*/
		if (_bases) {
			delete[] _bases;
		}
		_bases = new uint8_t[_seqSize];
		if (_bases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		/*reverse strand*/
		if (_rbases) {
			delete[] _rbases;
		}
		_rbases = new uint8_t[_seqSize];
		if (_rbases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}

		/*allocate space for quality scores*/
		if (quals) {
			if (_quals) {
				delete[] _quals;
			}
			_quals = new uint8_t[_seqSize];
			if (_quals == NULL) {
				Utils::exit("Memory allocation failed in function %s line %d\n",
						__FUNCTION__, __LINE__);
			}
		}
	}
}
void Sequence::print(FILE* file) {
	//print the sequence name
	if (_quals) {
		fputc('@', file);
	} else {
		fputc('>', file);
	}
	fprintf(file, "%s\n", _name);

	//print the query sequence
	for (uint32_t i = 0; i < _length; ++i) {
		fputc(decode(_bases[i]), file);
	}
	fputc('\n', file);

	//print the quality scores if available
	if (_quals) {
		/*printout comments*/
		fputc('+', file);
		fputc('\n', file);
		for (uint32_t i = 0; i < _length; ++i) {
			fputc(_quals[i], file);
		}
	}
	fputc('\n', file);
}

size_t Sequence::parseFastqSeq(char *seq) {

	int ch;
	int nprimers;
	bool trimmable;

	bool _colorspace = false; // For now we only focus in base space
	//find the header
	ch = *seq;

	// Check, In FASTQ format the read name must start with @
	while (ch != '\n' && ch != '@') {
		seq++;
		ch = *seq;
	}

	if ((ch == '\n') || (ch == EOF)) {
		return 0; //reach the end of file
	}

	//=============== 1st line. Read the sequence name (only one line)
	_length = 0;
	//while ((ch = myfgetc(_fp)) != -1 && ch != '\n') {

	// Pointer advance
	if(ch != '@') {
		seq++;
		ch = *seq;
	}

	while((ch != '\n') || (ch != EOF)) {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (isspace(ch)) {
			ch = '\0';
		}
		_buffer[_length++] = ch;

		seq++;
		ch = *seq;
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
	setNameSize(_length + 1); /*adjust the name buffer size*/
	strcpy((char*) _name, (char*) _buffer);

	//read the sequence bases
	_length = 0;
	nprimers = 0;
	do {


		//filter out the blank lines
		while (ch != -1 && (ch == '\r' || ch == '\n')) {
			seq++;
			ch = *seq;
		}

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
	_length = _length - (trimmable ? 1 : 0);
	setSequenceSize(_length, true);
	memcpy(_bases, _buffer + (trimmable ? 1 : 0), _length);

	//read the comment line (only one line, the 3rd, which is started by a '+')
	while (ch != -1 && ch != '\n') {
		seq++;
		ch = *seq;
	}

	if (ch == -1)
		Utils::exit("Incomplete FASTQ file\n");
	else if(ch == '\n') {
		seq++;
		ch = *seq;
	}

	// =============== 4th line read the quality scores
	_length = 0;



	while (ch != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (ch >= 33 && ch <= 127) {
			_buffer[_length++] = ch;
		}

		seq++;
		ch = *seq;
	}

	/*for base-space reads*/
	if (!_colorspace && _length != _length) {
		Utils::exit(
				"The number of bases is not equal to the number of quality scores\n");
	}
	if (_colorspace) {
		int32_t i, j;
		/*copy the base quality scores from the 3'end o the 5'end*/
		for(i = _length - 1, j = _length - 1; i >= 0 && j >= 0; --i, --j){
			_quals[j] = _buffer[i];
		}
		/*pseudo base quality score*/
		for(;j >= 0; --j){
			_quals[j] = '!';
		}
	}else{
		memcpy(_quals, _buffer, _length);
	}

	return _length;

}

void Sequence::resizeBuffer(size_t nsize) {
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
