#pragma once

#include <vector>
#include <string>

class Sequence
{
public:
	std::vector<std::string> names = {};
	std::vector<std::string> comments = {};
	std::vector<std::string> seqs = {};

	Sequence() = default;
	
	static Sequence fromFile(std::string filename);
	void stripNonCanonical();
	void toupper();
};
