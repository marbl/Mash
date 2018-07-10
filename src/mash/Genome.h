#pragma once

#include <string>
#include <vector>
#include "Sequence.h"

struct Genome
{
	std::string name = {};
	std::vector<Sequence> sequences = {};

	static Genome fromFile(std::string filename);
};

