// Copyright (c) 2009-2012 The Bitcoin developers
// Copyright (c) 2011-2012 Litecoin Developers
// Copyright (c) 2013 nanotoken Developers

// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <boost/foreach.hpp>

#include "checkpoints.h"

#include "main.h"
#include "uint256.h"

namespace Checkpoints
{
    typedef std::map<int, uint256> MapCheckpoints;

    //
    // What makes a good checkpoint block?
    // + Is surrounded by blocks with reasonable timestamps
    //   (no blocks before with a timestamp after, none after with
    //    timestamp before)
    // + Contains no strange transactions
    //
    static MapCheckpoints mapCheckpoints =
            boost::assign::map_list_of
            (     0, uint256("0xc843223185b61ec82a9b3e80c3bb6b51d5e2e185a3d7f1fa7a30bd0ffc9b2131"))
            (     3461, uint256("0xe3c1a28ce9d537f5b9a5458d517661a4b4e22fda47fcac4d8d5f33b99ce0de70"))
            (     3490, uint256("0x39367f1d2e8a390c4c2b28ae7b31d94edb0d62cca93dd23e294351df528e673c"))
            (     3560, uint256("0xf4496f8ed3b2f018f60fdfc46338a60163b853fc5c7310ebd64195851c375f3b"))
    	    (     7331, uint256("0x5aa90311213601db06c1bcf1ed2b795c2a625f25652afd82a9c91a5a06ee5cfe"))
            (     18400, uint256("0x2fc942062cd455b01790868cab0afe59e0a8de30edd2b8966b726aacf07fcc24"))
            (     28901, uint256("0x98b828ec3342be32b2de70c4e28f3b0b1cf54fb243a8a65c0f8660672434c8c2"))
			(     91598, uint256("0xe93ee538c8701212161b8abdf39137070c34a4918badf90bf59f3968128a1622"))
            ;


bool CheckBlock(int nHeight, const uint256& hash)
    {
        if (fTestNet) return true; // Testnet has no checkpoints

        MapCheckpoints::const_iterator i = mapCheckpoints.find(nHeight);
        if (i == mapCheckpoints.end()) return true;
        // return hash == i->second;
      return true;
    }

    int GetTotalBlocksEstimate()
    {
        if (fTestNet) return 0;
   
        // return mapCheckpoints.rbegin()->first;
      return 0;
    }

    CBlockIndex* GetLastCheckpoint(const std::map<uint256, CBlockIndex*>& mapBlockIndex)
    {
        if (fTestNet) return NULL;

        BOOST_REVERSE_FOREACH(const MapCheckpoints::value_type& i, mapCheckpoints)
        {
            const uint256& hash = i.second;
            std::map<uint256, CBlockIndex*>::const_iterator t = mapBlockIndex.find(hash);
            if (t != mapBlockIndex.end())
                // return t->second;
            return NULL;
        }
        return NULL;
    }
}
