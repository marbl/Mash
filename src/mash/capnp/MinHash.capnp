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
	
	struct LocusList
	{
		struct Locus
		{
			sequence @0 : UInt32;
			position @1 : UInt32;
			hash @2 : UInt64;
		}
		
		loci @0 : List(Locus);
	}
	
	referenceList @0 : ReferenceList;
	locusList @1 : LocusList;
	kmerSize @2 : UInt32;
	windowSize @3 : UInt32;
	minHashesPerWindow @4 : UInt32;
}
