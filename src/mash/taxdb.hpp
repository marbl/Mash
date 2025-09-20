#ifndef TAXD_DB_H_
#define TAXD_DB_H_

// Florian Breitiweser

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <iomanip>
#include <sstream>
#include <stdexcept>


using TaxID = uint64_t;

using std::vector;
using std::string;

namespace mash {

class TaxEntry {
  public:
    TaxID taxID;
    TaxEntry* parent;
    vector<TaxEntry*> children;

    string rank;
    string name; // scientific name

    TaxEntry() : taxID(0), parent(NULL) {}

    TaxEntry(TaxID taxID, string rank) : taxID(taxID), rank(rank) {}
};

struct TaxCounts {
  uint64_t cladeCount = 0;
  uint64_t taxCount = 0;
  uint64_t taxHashCount = 0;
  uint64_t cladeHashCount = 0;
  vector<TaxID> children;
};

class TaxDB {
  public:
    TaxDB(const string namesDumpFileName, const string nodesDumpFileName);
    TaxDB(const string inFileName);
    TaxDB();

    void writeTaxIndex(std::ostream & outs) const;
    void readTaxIndex(const string inFileName);

    TaxID getLowestCommonAncestor(TaxID a, TaxID b) const;
    string getLineage(TaxID taxID) const;
    string getMetaPhlAnLineage(TaxID taxID) const;
    TaxEntry const * getEntry(TaxID taxID) const;

    unordered_map<TaxID, TaxEntry> entries;

    void writeReport(FILE* FP, const unordered_map<TaxID, TaxCounts> & counts, 
                     unsigned long totalCounts, 
                     unsigned long totalHashCounts, 
                     TaxID taxID = 0, int depth = 0);

  private:
    unordered_map<TaxEntry*, TaxID> parseNodesDump(const string nodesDumpFile);
    void parseNamesDump(const string namesDumpFile);
};

TaxEntry const * TaxDB::getEntry(TaxID taxID) const {
  auto it = entries.find(taxID);
  if (it == entries.end()) {
    cerr << "Couldn't find tax entry with taxID " << taxID << endl;
    return NULL;
  } else {
    return &it->second;
  }
}


TaxDB::TaxDB(const string namesDumpFileName, const string nodesDumpFileName) {
  unordered_map<TaxEntry*, TaxID> parentMap = parseNodesDump(nodesDumpFileName);

  // set parent links correctly
  for (auto const & c : parentMap) {
    if (c.first->taxID != c.second) {
      auto p = entries.find(c.second);
      if (p == entries.end()) {
         cerr << "Could not find parent with tax ID " << c.second << " for tax ID " << c.first->taxID << endl;
      } else {
        c.first->parent = &p->second;
      }
    } else {
      c.first->parent = NULL;
    }
  }
  parseNamesDump(namesDumpFileName);
  cerr << "   " << entries.size() << " distinct taxa\n";
}

unordered_map<TaxEntry*,TaxID> TaxDB::parseNodesDump(const string nodesDumpFileName) {
  std::ifstream nodesDumpFile(nodesDumpFileName);
  if (!nodesDumpFile.is_open())
    throw std::runtime_error("unable to open nodes file");

  string line;
  TaxID taxID;
  TaxID parentTaxID;
  string rank;
  char delim;
  unordered_map<TaxEntry*,TaxID> parentMap;

  while (nodesDumpFile >> taxID >> delim >> parentTaxID >> delim) {
    nodesDumpFile.ignore(1);
    getline(nodesDumpFile, rank, '\t');
    // TODO: Insert
    //auto res = entries.insert(taxID, TaxEntry(taxID, rank));
    auto res = entries.emplace(taxID, TaxEntry(taxID, rank));
    parentMap.emplace(&res.first->second, parentTaxID);
    nodesDumpFile.ignore(2560, '\n');
  }
  return parentMap;
}

void TaxDB::parseNamesDump(const string namesDumpFileName) {
  std::ifstream namesDumpFile(namesDumpFileName);
  if (!namesDumpFile.is_open())
    throw std::runtime_error("unable to open names file");
  string line;

  TaxID taxID;
  string name, type;
  char delim;
  while (namesDumpFile >> taxID) {
    namesDumpFile.ignore(3);
    getline(namesDumpFile, name, '\t');
    namesDumpFile.ignore(3);
    namesDumpFile.ignore(256, '|');
    namesDumpFile.ignore(1);
    getline(namesDumpFile, type, '\t');

    if (type == "scientific name") {
      auto entryIt = entries.find(taxID);
      if (entryIt == entries.end()) {
        cerr << "Entry for " << taxID << " does not exist - it should!" << '\n';
      } else {
        entryIt->second.name = name;
      }
    }
    namesDumpFile.ignore(2560, '\n');
  }
}

TaxID TaxDB::getLowestCommonAncestor(TaxID a, TaxID b) const {
  if (b == 0) { return a; }
  if (a == 0) { return b; } 

  // create a path from a to the root
  unordered_set<TaxEntry const *> a_path;
  std::unordered_map<TaxID, TaxEntry>::const_iterator ta = entries.find(a);
  if (ta == entries.end()) {
    cerr << "TaxID " << a << " not in database - ignoring it.\n";
    return 1;
  }

  std::unordered_map<TaxID, TaxEntry>::const_iterator tb = entries.find(b);
  if (tb == entries.end()) {
    cerr << "TaxID " << b << " not in database - ignoring it.\n";
    return 1;
  }
  TaxEntry const * pta = &(ta->second);
  while (pta != NULL && pta->taxID > 1 && pta->parent != NULL) {
    if (pta->taxID == b) { return b; }
    a_path.insert(pta);
    pta = pta->parent;
  }
  TaxEntry const * ptb = &(tb->second);
  // search for b in the path from a to the root
  while (ptb->taxID > 0 && ptb->parent != NULL) {
    if (a_path.count(ptb)) {
      return ptb->taxID;
    }
    ptb = ptb->parent;
  }
  return 1;
}

void TaxDB::writeReport(FILE* FP,
			const unordered_map<TaxID, TaxCounts> & counts,
			unsigned long totalCounts,
			unsigned long totalHashCounts,
			TaxID taxID, int depth) {

	unordered_map<TaxID, TaxCounts>::const_iterator it = counts.find(taxID);
	unsigned int cladeCount = it == counts.end()? 0 : it->second.cladeCount;
	unsigned int cladeHashCount = it == counts.end()? 0 : it->second.cladeHashCount;
	unsigned int taxCount = it == counts.end()? 0 : it->second.taxCount;
	unsigned int taxHashCount = it == counts.end()? 0 : it->second.taxHashCount;
	if (taxID == 0) {
    // TODO: Write header?
    // identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment
    fprintf(FP, "%%\thashes\ttaxHashes\thashesDB\ttaxHashesDB\ttaxID\trank\tname\n");
		if (cladeCount > 0) { // Should not happen
			fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
					100 * cladeCount / double(totalCounts),
					cladeCount, taxCount);
		}
	  writeReport(FP, counts, totalCounts, totalHashCounts, 1, 0);
	} else {
		if (cladeCount == 0) {
			return;
		}
		TaxEntry const * taxon = getEntry(taxID);
		fprintf(FP, "%.4f\t%i\t%i\t%i\t%i\t%s\t%lu\t%s%s\n",
				100*cladeCount/double(totalCounts), 
        cladeCount, 
        taxCount, 
        cladeHashCount,
        taxHashCount,
				taxon->rank.c_str(), taxID, std::string(2*depth, ' ').c_str(), taxon->name.c_str());

		std::vector<TaxID> children = it->second.children;
		std::sort(children.begin(), children.end(), [&](int a, int b) { return counts.at(a).cladeCount > counts.at(b).cladeCount; });
		for (TaxID childTaxId : children) {
			if (counts.count(childTaxId)) {
				writeReport(FP, counts, totalCounts, totalHashCounts, childTaxId, depth + 1);
			} else {
				break;
			}
		}
	}
}

/*
string TaxDB::getLineage(TaxEntry tax) const {
  string lineage;
  while (true) {
    // 131567 = Cellular organisms
    if (taxID != 131567) {
      if (lineage.size()) lineage.insert(0, "; ");
      lineage.insert(0, getScientificName(taxID));
      if (getRank(taxID) == "species") lineage.clear();
    }
    taxID = getParentTaxID(taxID);
    if (taxID == 0) {
      if (lineage.size()) lineage.append(".");
      break;
    }
  }
  return lineage;
}

string TaxDB::getMetaPhlAnLineage(TaxID taxID) const {
  string rank = getRank(taxID);
  if (rank == "superphylum") return string();
  string lineage;
  while (true) {
    // 131567 = Cellular organisms
    //if (taxID != 131567) {
      string rank = getRank(taxID);
      if (rank == "species") {
  lineage.insert(0, "|s__");
  lineage.insert(4, getScientificName(taxID));
      } else if (rank == "genus") {
  lineage.insert(0, "|g__");
  lineage.insert(4, getScientificName(taxID));
      } else if (rank == "family") {
  lineage.insert(0, "|f__");
  lineage.insert(4, getScientificName(taxID));
      } else if (rank == "order") {
  lineage.insert(0, "|o__");
  lineage.insert(4, getScientificName(taxID));
      } else if (rank == "class") {
  lineage.insert(0, "|c__");
  lineage.insert(4, getScientificName(taxID));
      } else if (rank == "phylum") {
  lineage.insert(0, "|p__");
  lineage.insert(4, getScientificName(taxID));
      } else if (rank == "superkingdom") {
  lineage.insert(0, "|k__");
  lineage.insert(4, getScientificName(taxID));
      } else {
  lineage.insert(0, "|-__");
  lineage.insert(4, getScientificName(taxID));

	 // }
    }
    taxID = getParentTaxID(taxID);
    if (taxID == 0) {
      break;
    }
  }
  std::replace(lineage.begin(), lineage.end(), ' ', '_');
  return lineage;
}
*/

}

#endif /* TAXD_DB_H_ */
