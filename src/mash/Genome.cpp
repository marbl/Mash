#include "kseq.h"
#include "Sequence.h"
#include "Genome.h"
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

std::string extract_name(const std::string &file_name);


KSEQ_INIT(gzFile, gzread)

Genome Genome::fromFile(std::string filename)
{
	auto ret = Genome{};
	ret.name = extract_name(filename);

	auto fd = filename == "-" ? 0 : open(filename.c_str(), O_RDONLY);

	gzFile fp = gzdopen(fd, "r");
	kseq_t *kseq = kseq_init(fp);

	int l;

	while ((l = kseq_read(kseq)) >= 0) {
		// TODO: find out if seq->seq.l can be trusted
		auto name = std::string(kseq->seq.s);
		auto sequence = std::string(kseq->seq.s);

		auto comment = std::string();
		if (kseq->comment.s) {
			comment += kseq->comment.s;
		}

		ret.sequences.emplace_back(name, comment, sequence);
	}

	kseq_destroy(kseq);
	gzclose(fp);
	close(fd);

	return ret;
}

std::string extract_name(const std::string &file_name)
{
	// find the last path separator
	auto left = file_name.rfind('/');
	left = (left == std::string::npos) ? 0 : left + 1;
	// left is the position one of to the right of the path separator

	// find the extension
	auto right = file_name.find('.', left);
	right = (right == std::string::npos) ? file_name.size() : right;

	// copy only the file name, not its path or extension
	return file_name.substr(left, right - left);
}

// void Sequence::stripNonCanonical()
// {
// 	static const std::array<char,256> alphabet = []{
// 		std::array<char,256> alphabet = {};
// 		alphabet['A'] = 'A';
// 		alphabet['C'] = 'C';
// 		alphabet['G'] = 'G';
// 		alphabet['T'] = 'T';
// 		alphabet['a'] = 'A';
// 		alphabet['c'] = 'C';
// 		alphabet['g'] = 'G';
// 		alphabet['t'] = 'T';
// 		return alphabet;
// 	}();

// 	for (auto& seq: seqs){
// 		auto it = std::remove_if(seq.begin(), seq.end(), [&](char c){
// 			return alphabet[c] == 0;
// 		});
// 		seq.erase(it, seq.end());
// 	};
// }

// void Sequence::toupper()
// {
// 	for (auto &seq: seqs) {
// 		std::transform(seq.begin(), seq.end(), seq.begin(), [](unsigned char c){
// 			return std::toupper(c);
// 		});
// 	}
// }
