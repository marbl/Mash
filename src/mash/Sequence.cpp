#include "kseq.h"
#include "Sequence.h"
#include <algorithm>
#include <array>
#include <string>
#include <vector>
#include <cctype>


void Sequence::stripNonCanonical()
{
	static const std::array<char,256> alphabet = []{
		std::array<char,256> alphabet = {};
		alphabet['A'] = 'A';
		alphabet['C'] = 'C';
		alphabet['G'] = 'G';
		alphabet['T'] = 'T';
		alphabet['a'] = 'A';
		alphabet['c'] = 'C';
		alphabet['g'] = 'G';
		alphabet['t'] = 'T';
		return alphabet;
	}();

	auto it = std::remove_if(sequence.begin(), sequence.end(), [&](char c){
		return alphabet[c] == 0;
	});
	sequence.erase(it, sequence.end());
}

void Sequence::toupper()
{
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), [](unsigned char c){
		return std::toupper(c);
	});
}
