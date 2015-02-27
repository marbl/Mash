using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("capnp");

@0xc4c8b1ada05e7704;

struct MinHash
{
	struct ReferenceList
	{
		struct Reference
		{
			sequence @0 : Text;
			quality @1 : Text;
			length @ 2 : UInt32;
			name @3 : Text;
			comment @4 : Text;
		}
		
		references @0 : List(Reference);
	}
	
	struct HashTable
	{
		struct HashBin
		{
			struct Locus
			{
				sequence @0 : UInt32;
				position @1 : UInt32;
			}
			
			hash @0 : UInt32;
			loci @ 1 : List(Locus);
		}
		
		hashBins @0 : List(HashBin);
	}
	
	referenceList @0 : ReferenceList;
	hashTable @1 : HashTable;
	kmerSize @2 : UInt32;
	compressionFactor @3 : Float32;
}
