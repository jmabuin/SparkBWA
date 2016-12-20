/*
 * SAMSpark.cpp
 *
 *  Created on: 28 nov, 2016
 *      Author: chema
 */

#include "SAMSpark.h"
#include "Utils.h"

SAMSpark::SAMSpark(Options* options, Genome* genome, char *str) {

	_genome = genome;
	_colorspace = options->isColorSpace();
	_pairingMode = options->getPairingMode();

	/*print out the header*/
	sprintf(str, "@HD\tVN:1.0\tSO:unsorted\n");
	char *buffer = (char *) malloc(sizeof(char) * BUFF_SIZE);

	/*print out the genome sequence information in SAM format*/
	_genome->genomeNamesOut(buffer);

	strcat(str,buffer);
    //fprintf(stderr,"[%s] Estamos aqui \n",__func__);
	/*print out the read group header information*/
	_groupID = NULL;
	_groupLB = NULL;


	if (options->getRgId().length() > 0) {
		_groupID = (char*) options->getRgId().c_str();
		/*tag ID*/
		sprintf(buffer, "@RG\tID:%s\tSM:%s", options->getRgId().c_str(),
				options->getRgSm().c_str());
		strcat(str,buffer);

		/*tag LB*/
		if (options->getRgLb().length() > 0) {
			sprintf(buffer, "\tLB:%s", options->getRgLb().c_str());
			strcat(str, buffer);
			_groupLB = (char*)options->getRgLb().c_str();
		}
		/*tag PL*/
		if (options->getRgPl().length() > 0) {
			sprintf(buffer, "\tPL:%s", options->getRgPl().c_str());
			strcat(str,buffer);
		}
		/*tag PU*/
		if (options->getRgPu().length() > 0) {
			sprintf(buffer, "\tPU:%s", options->getRgPu().c_str());
			strcat(str,buffer);
		}
		/*tag CN*/
		if (options->getRgCn().length() > 0) {
			sprintf(buffer, "\tCN:%s", options->getRgCn().c_str());
			strcat(str,buffer);
		}
		/*tag DS*/
		if (options->getRgDs().length() > 0) {
			sprintf(buffer, "\tDS:%s", options->getRgDs().c_str());
			strcat(str,buffer);
		}
		/*tag DT*/
		if (options->getRgDt().length() > 0) {
			sprintf(buffer, "\tDT:%s", options->getRgDt().c_str());
			strcat(str,buffer);
		}
		/*tag PI*/
		if (options->getRgPi().length() > 0) {
			sprintf(buffer, "\tPI:%s", options->getRgPi().c_str());
			strcat(str,buffer);
		}
		/*end of the line*/
		//fputc('\n', buffer);
		strcat(str,"\n");
	}

	/*print out the aligner information*/
	sprintf(buffer, "@PG\tID:%s\tVN:%s\n", PROGRAM_NAME, PROGRAM_VERSION);
	strcat(str,buffer);

	free(buffer);

}

SAMSpark::~SAMSpark() {
	// TODO Auto-generated destructor stub
}

/*print the unaligned read information*/
void SAMSpark::print(Sequence& seq, char *str, int _flags) {
	int flags = SAM_FSU | _flags;

	char buffer[BUFF_SIZE];
	/*print out query name, bitwise-flag, and reference sequence name*/
	sprintf(buffer, "%s\t%d\t*\t", seq._name, flags);
	strcat(str, buffer);
	strcat(str, buffer);

	//print 1-based leftmost mapping position, mapping quality (phred-scale)*/
	sprintf(buffer, "0\t0\t");
	strcat(str, buffer);

	//print extended CIGAR
	strcat(str, "*\t");

	//print paired-end information, INAVAILABLE
	sprintf(buffer, "*\t0\t0\t");
	strcat(str, buffer);

	//print the query sequence
	uint8_t* bases = seq._bases;
	for (uint32_t i = 0; i < seq._length; ++i) {
		//fputc(decode(bases[i]), buffer);
		strcat(str, (char *) decode(bases[i]));
	}
	//fputc('\t', buffer);
	strcat(str, "\t");

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals = seq._quals;
		for (uint32_t i = 0; i < seq._length; ++i) {
			//fputc(*quals, buffer);
			strcat(str, (char *)*quals);
			++quals;
		}
	} else {
		//fputc('*', buffer);
		strcat(str, "*");
	}
	/*print tags*/
	sprintf(buffer, "\tPG:Z:%s", PROGRAM_NAME);
	strcat(str, buffer);

	if (_groupID) {
		sprintf(buffer, "\tRG:Z:%s", _groupID);
		strcat(str, buffer);
	}
	if(_groupLB){
		sprintf(buffer, "\tLB:Z:%s", _groupLB);
		strcat(str, buffer);
	}

	//fputc('\n', buffer);
	strcat(str, "\n");
	free(buffer);
}
/*print the alignment information through the dynamic programming*/
void SAMSpark::print(Sequence& seq, Mapping& mapping, char *str, int _flags) {
	int32_t flags = _flags;
	uint32_t length;
	uint8_t* bases;
	CigarAlign* align = mapping._align;

	char buffer[BUFF_SIZE];

	/*check the sequence strand*/
	if (mapping._strand) {
		flags |= SAM_FSR;
	}

	/*get sequence length*/
	length = _colorspace ? mapping._seqLength : seq._length;

	/*print out query name, bitwise-flag, and reference sequence name*/
	sprintf(buffer, "%s\t%d\t%s\t", seq._name, flags,
			_genome->getGenomeName(mapping._genomeIndex));
	strcat(str, buffer);

	//print 1-based leftmost mapping position, mapping quality (phred-scale)*/
	sprintf(buffer, "%u\t%d\t", (uint32_t) mapping._position, mapping._mapQual);
	strcat(str, buffer);

	//print extended CIGAR
	if (align->getCigar()) {
		align->cigarOut(buffer);
	} else {
		sprintf(buffer, "%dM", length);
		strcat(str, buffer);
	}
	//fputc('\t', buffer);
	strcat(str, "\t");

	//print paired-end information, INAVAILABLE
	sprintf(buffer, "*\t0\t0\t");
	strcat(str, buffer);

	//print the query sequence
	if (_colorspace) {
		bases = (mapping._strand == 0) ? mapping._data : mapping._data + length;
	} else {
		bases = (mapping._strand == 0) ? seq._bases : seq._rbases;
	}

	for (uint32_t i = 0; i < length; ++i) {
		//fputc(decode(bases[i]), buffer);
		strcat(str, (char *)decode(bases[i]));
	}
	//fputc('\t', buffer);
	strcat(str, "\t");

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals;
		if (mapping._strand == 0) {
			quals = _colorspace ? mapping._data + 2 * length : seq._quals;
			for (uint32_t i = 0; i < length; ++i) {
				//fputc(*quals, buffer);
				strcat(str, (char *)*quals);
				++quals;
			}
		} else {
			quals = _colorspace ?
					mapping._data + 3 * length - 1 :
					seq._quals + seq._length - 1;
			for (uint32_t i = 0; i < length; ++i) {
				//fputc(*quals, buffer);
				strcat(str, (char *)*quals);
				--quals;
			}
		}
	} else {
		//fputc('*', buffer);
		strcat(str,"*");
	}

	/*print tags*/
	sprintf(buffer, "\tPG:Z:%s", PROGRAM_NAME);
	strcat(str, buffer);
	if (_groupID) {
		sprintf(buffer, "\tRG:Z:%s", _groupID);
		strcat(str, buffer);
	}
	if(_groupLB){
		sprintf(buffer, "\tLB:Z:%s", _groupLB);
		strcat(str, buffer);
	}
	sprintf(buffer, "\tNM:i:%d", align->getEditDistance());
	strcat(str, buffer);
	sprintf(buffer, "\tAS:i:%d", align->getAlignScore());
	strcat(str, buffer);

	/*end of line*/
	//fputc('\n', buffer);
	strcat(str,"\n");
	free(buffer);

}
/*print out the alignment with itself unaligned and its mate aligned*/
void SAMSpark::printPaired(Sequence& seq, Mapping& mate, char *str, int _flags) {
	int flags = _flags;

	/*set the unmapped flag*/
	flags |= SAM_FSU;

	char buffer[BUFF_SIZE];

	//check the strand of the mate
	if (mate._strand) {
		flags |= SAM_FMR;
	}

	/*print out query name, bitwise-flag, and reference sequence name*/
	sprintf(buffer, "%s\t%d\t*\t", seq._name, flags);
	strcat(str, buffer);
	//print 1-based leftmost mapping position, mapping quality (phred-scale)*/
	sprintf(buffer, "0\t0\t");
	strcat(str, buffer);

	//print extended CIGAR
	//fputs("*\t", buffer);
	strcat(str, "*\t");

	/****************************
	 print the mate information
	 1-base mate mapping position and estimated insert size)
	 *****************************/
	//print the reference sequence mapped by the mate
	sprintf(buffer, "%s\t", _genome->getGenomeName(mate._genomeIndex));
	strcat(str, buffer);
	//print the mapping position of the mate and the distance

	int64_t distance = 0; /*means not properly paired*/
	//Utils::log("self %d mate %d length %d strand %d-%d\n", self._position, mate._position, length, self._strand, mate._strand);

	sprintf(buffer, "%d\t%ld\t", (int) mate._position, distance);
	strcat(str, buffer);

	//print the query sequence
	uint8_t* bases = seq._bases;
	for (uint32_t i = 0; i < seq._length; ++i) {
		//fputc(decode(bases[i]), buffer);
		strcat(str, (char *)bases[i]);
	}
	//fputc('\t', buffer);
	strcat(str, "\t");

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals = seq._quals;
		for (uint32_t i = 0; i < seq._length; ++i) {
			//fputc(*quals, buffer);
			strcat(str, (char *)*quals);
			++quals;
		}
	} else {
		//fputc('*', buffer);
		strcat(str, "*");
	}
	/*print tags*/
	sprintf(buffer, "\tPG:Z:%s", PROGRAM_NAME);
	strcat(str, buffer);

	if (_groupID) {
		sprintf(buffer, "\tRG:Z:%s", _groupID);
		strcat(str, buffer);
	}
	if(_groupLB){
		sprintf(buffer, "\tLB:Z:%s", _groupLB);
		strcat(str, buffer);
	}
	//fputc('\n', buffer);
	strcat(str, "\n");
	free(buffer);
}


void SAMSpark::printPaired(Sequence& seq, Mapping& self, Mapping& mate, char *str,
		bool properlyPaired, int _flags) {
	int flags = _flags;
	uint32_t length;
	uint8_t* bases;
	CigarAlign* align = self._align;

	char buffer[BUFF_SIZE];

	/*the reads are paired*/
	if (properlyPaired) {
		flags |= SAM_FPP;
	}

	//check the strand of itself
	if (self._strand) {
		flags |= SAM_FSR;
	}

	//check the strand of the mate
	if (mate._strand) {
		flags |= SAM_FMR;
	}

	/*get sequence length*/
	length = _colorspace ? self._seqLength : seq._length;

	//print query-name, bitwise-flag, and reference-sequence-name
	sprintf(buffer, "%s\t%d\t%s\t", seq._name, flags,
			_genome->getGenomeName(self._genomeIndex));
	strcat(str, buffer);

	//print 1-based leftmost mapping position, mapping quality (phred-scaled)
	sprintf(buffer, "%d\t%d\t", (int) self._position, self._mapQual);
	strcat(str, buffer);

	//print extended CIGAR if applicable
	if (align->getCigar()) {
		align->cigarOut(buffer);
	} else {
		sprintf(buffer, "%dM", length);
		strcat(str, buffer);
	}
	//fputc('\t', buffer);
	strcat(str, "\t");

	/****************************
	 print the mate information
	 1-base mate mapping position and estimated insert size)
	 *****************************/
	//print the reference sequence mapped by the mate
	sprintf(buffer, "%s\t",
			(self._genomeIndex == mate._genomeIndex) ?
					"=" : _genome->getGenomeName(mate._genomeIndex));
	strcat(str, buffer);

	//print the mapping position of the mate and the distance

	int64_t distance = 0; /*means not properly paired*/
	if (properlyPaired) {
		if (_pairingMode == MATE_PAIRED_READS) {
			distance = (int64_t) self._position - mate._position;
		} else {
			if (self._strand == 0) {
				distance = (int64_t) self._position - mate._position - length;
			} else {
				distance = (int64_t) self._position + length - mate._position;
			}
		}
	}
	//Utils::log("self %d mate %d length %d strand %d-%d\n", self._position, mate._position, length, self._strand, mate._strand);

	sprintf(buffer, "%d\t%ld\t", (int) mate._position, distance);
	strcat(str, buffer);

	//print the query sequence on the same strand as the reference sequence
	if (_colorspace) {
		bases = (self._strand == 0) ? self._data : self._data + length;
	} else {
		bases = (self._strand == 0) ? seq._bases : seq._rbases;
	}
	for (uint32_t i = 0; i < length; ++i) {
		//fputc(decode(bases[i]), buffer);
		strcat(str, (char *)decode(bases[i]));
	}
	//fputc('\t', buffer);
	strcat(str, "\t");

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals;
		if (self._strand == 0) {
			quals = _colorspace ? self._data + 2 * length : seq._quals;
			for (uint32_t i = 0; i < length; ++i) {
				//fputc(*quals, buffer);
				strcat(str, (char*)*quals);
				++quals;
			}
		} else {
			quals = _colorspace ?
					self._data + 3 * length - 1 : seq._quals + seq._length - 1;
			for (uint32_t i = 0; i < length; ++i) {
				//fputc(*quals, buffer);
				strcat(str, (char*)*quals);
				--quals;
			}
		}
	} else {
		//fputc('*', buffer);
		strcat(str, "*");
	}

	/*print tags*/
	sprintf(buffer, "\tPG:Z:%s", PROGRAM_NAME);
	strcat(str, buffer);
	if (_groupID) {
		sprintf(buffer, "\tRG:Z:%s", _groupID);
		strcat(str, buffer);
	}
	if(_groupLB){
		sprintf(buffer, "\tLB:Z:%s", _groupLB);
		strcat(str, buffer);
	}
	sprintf(buffer, "\tNM:i:%d", align->getEditDistance());
	strcat(str, buffer);
	sprintf(buffer, "\tAS:i:%d", align->getAlignScore());
	strcat(str, buffer);
	strcat(str,"\n");
	free(buffer);
}
