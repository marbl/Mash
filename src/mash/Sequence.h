#pragma once

#include <vector>
#include <string>

class Sequence
{
public:
	std::string name = {};
	std::string comment = {};
	std::string sequence = {};

	Sequence() = default;
	Sequence(std::string n, std::string c, std::string s):
		name(std::move(n)), comment(std::move(c)), sequence(std::move(s))
	{
	}
	
	static Sequence fromFile(std::string filename);
	void stripNonCanonical();
	void toupper();
};
