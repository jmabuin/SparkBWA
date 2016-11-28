/*
 * SAMSpark.h
 *
 *  Created on: 28 nov, 2016
 *      Author: chema
 */

#ifndef CORE_SAMSPARK_H_
#define CORE_SAMSPARK_H_

#include "Genome.h"
#include "Options.h"
#include "Sequence.h"
#include "Mapping.h"

class SAMSpark {
public:
	SAMSpark(Options* options, Genome* genome, char *str);
	virtual ~SAMSpark();

	void print(Sequence& seq,char *str, int flags = 0);
	void print(Sequence& seq, Mapping& mapping, char *str, int flags = 0);
	void printPaired(Sequence& seq, Mapping& mate, char *str, int flags = 0);
	void printPaired(Sequence& seq, Mapping& self, Mapping& mate, char *str, bool properlyPaired, int _flags);

private:
	Genome* _genome;
	bool _colorspace;
	bool _pairingMode;
	char* _groupID;
	char* _groupLB;

	unsigned int BUFF_SIZE = 4096;
};

#endif /* CORE_SAMSPARK_H_ */
