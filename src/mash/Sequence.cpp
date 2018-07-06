#include "kseq.h"
#include "Sequence.h"
#include <algorithm>
#include <array>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <zlib.h>
#include <cctype>

KSEQ_INIT(gzFile, gzread)

Sequence Sequence::fromFile(std::string filename)
{
	auto ret = Sequence{};

	auto fd = filename == "-" ? 0 : open(filename.c_str(), O_RDONLY);

	gzFile fp = gzdopen(fd, "r");
	kseq_t *seq = kseq_init(fp);

	int l;

	while ((l = kseq_read(seq)) >= 0) {
		ret.seqs.emplace_back(seq->seq.s); // TODO: find out if seq->seq.l can be trusted
		ret.names.emplace_back(seq->name.s);
		if (seq->comment.s) {
			ret.comments.emplace_back(seq->comment.s);
		} else {
			ret.comments.push_back("");
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
	close(fd);

	return ret;
}

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

	for (auto& seq: seqs){
		auto it = std::remove_if(seq.begin(), seq.end(), [&](char c){
			return alphabet[c] == 0;
		});
		seq.erase(it, seq.end());
	};
}

void Sequence::toupper()
{
	for (auto &seq: seqs) {
		std::transform(seq.begin(), seq.end(), seq.begin(), [](unsigned char c){
			return std::toupper(c);
		});
	}
}
