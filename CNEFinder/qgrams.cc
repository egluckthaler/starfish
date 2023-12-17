/* ================================================================= *
 *  qgrams.cc                                       *
 *                                                                   *
 *  E-MEM: An efficient (MUMmer-like) tool to retrieve Maximum Exact *
 *         Matches using hashing based algorithm                     *
 *                                                                   *
 *  Copyright (c) 2014, Nilesh Khiste                                *
 *  All rights reserved                                              *
 *                                                                   *
 *  This program is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public License as   *
 *  published by the Free Software Foundation, either version 3 of   *
 *  the License, or (at your option) any later version.              *
 *                                                                   *
 *  This program is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 *  GNU General Public License for more details.                     *
 *                                                                   *
 *  You should have received a copy of the GNU General Public        *
 *  License along with this program.                                 *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'LICENSE', which is part of this source code package.       *
 * ================================================================= */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <iterator>
#include <omp.h>
#include "boost/algorithm/string.hpp"
#include "boost/tokenizer.hpp"
#include <sys/stat.h>
#include "qgrams.h"
#include "cnef.h"
#include "file.h"
#include "qlist.h"

using namespace std;
using namespace boost;

/* 
 * Function builds a kmer hash for a reference sequence.
 * Input: empty refHash
 * Output: populated refHash
 */
void buildRefHash(Knode* &refHash, uint64_t totalBits, seqFileReadInfo &RefFile)
{
    uint64_t j=0;
    uint64_t currKmerPos=0, currKmer=0;
    int32_t offset=0; 
    int nextKmerPosition = commonData::minMemLen - commonData::kmerSize + 2;
    vector<mapObject>::iterator it;
    it = upper_bound(RefFile.blockOfNs.begin(), RefFile.blockOfNs.end(), currKmerPos, mapObject()); 
    while (currKmerPos<=totalBits)
    {
        if (currKmerPos + commonData::kmerSize - 2 > totalBits)
            break;

        if(RefFile.checkKmerForNs(currKmerPos, it)){
            currKmerPos+=nextKmerPosition; // Move L-K+2 bits = 50-28+1=23 char = 46 bits
            continue;
        }

        offset = currKmerPos%DATATYPE_WIDTH;
        j=currKmerPos/DATATYPE_WIDTH; // next loc in binReads 
        
        currKmer = RefFile.binReads[j];
        currKmer <<= offset;

        if (offset > DATATYPE_WIDTH-commonData::kmerSize) // Kmer split in two integers
            currKmer |= ((RefFile.binReads[j+1] & global_mask_left[(commonData::kmerSize-(DATATYPE_WIDTH-offset))/2 -1])>>(DATATYPE_WIDTH-offset));
        else
            currKmer &= global_mask_left[commonData::kmerSize/2 - 1];
        /* Add kmer to the hash table */
        refHash->addKmerNode(currKmer, currKmerPos);
        currKmerPos+=nextKmerPosition; // Move L-K+2 bits = 50-28+1=23 char = 46 bits
    }
}

/* 
 * Function extends the kmer match in left/right direction for 
 * possible MEMs.
 * Input: currRPos : current position of matching reference Kmer 
 * Input: currRPos : current position of matching query Kmer 
 * Input: totalRBits : total number of bits in reference  
 * Input: totalQBits : total number of bits in query  
 * Input: name : reference sequence string for output  
 *
 */
void helperReportMem(uint64_t &currRPos, uint64_t &currQPos, uint64_t totalRBits, uint64_t totalQBits, queryList* &currQueryMEMs, std::unordered_multimap <uint64_t, uint64_t> &currMEMs, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, mapObject &RefNpos, mapObject &QueryNpos, uint32_t &revComplement)
{
    /*
     * lRef and lQue are local variables for left extension of
     * reference and query sequence respectively. rRef and rQue
     * are their right counterparts.
     */
    uint64_t lRef=currRPos, lQue=currQPos; // Keeping lRef on currRPos-this makes offset computation simpler
    uint64_t offsetR=0,offsetQ=0;
    uint64_t rRef=currRPos+commonData::kmerSize, rQue=currQPos+commonData::kmerSize; // one character ahead of current match
    uint64_t currR=0, currQ=0;
    int i=0,j=0,mismatch=0;
    uint64_t matchSize=0;

    if (!(((QueryNpos.left==0x1)?true:QueryNpos.left<=lQue) && rQue<=QueryNpos.right)) 
        QueryFile.getKmerLeftnRightBoundForNs(lQue, QueryNpos); 

    if (!(((RefNpos.left==0x1)?true:RefNpos.left<=lRef) && rRef<=RefNpos.right)) 
        RefFile.getKmerLeftnRightBoundForNs(lRef, RefNpos); 

    if (RefNpos.right-((RefNpos.left==0x1)?0:RefNpos.left)+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    if (QueryNpos.right-((QueryNpos.left==0x1)?0:QueryNpos.left)+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    //match towards left
    while (lRef && lQue && ((QueryNpos.left==0x1)?true:QueryNpos.left<=lQue) && ((RefNpos.left==0x1)?true:RefNpos.left<=lRef))
    {
        if (!mismatch)
        {
            offsetR=(lRef)%DATATYPE_WIDTH;
            i=(lRef)/DATATYPE_WIDTH;
            offsetQ=(lQue)%DATATYPE_WIDTH;
            j=(lQue)/DATATYPE_WIDTH;
       
            if (offsetR > offsetQ)
                matchSize = offsetQ;
            else 
                matchSize = offsetR;

            if (!matchSize)
                matchSize=2;

            if ((QueryNpos.left!=0x1) && (matchSize > lQue-QueryNpos.left))
                matchSize = lQue-QueryNpos.left;
            
            if ((RefNpos.left!=0x1) && (matchSize > lRef-RefNpos.left))
                matchSize = lRef-RefNpos.left;
            
            if (!matchSize)
                break;

            /* 
             * There will never be case with offset=0 and i=0 because
             * i=0 happens only when lRef=0 and in that case we do not 
             * enter this loop.
             */
            currR = RefFile.binReads[offsetR?i:i-1];
            currR >>= DATATYPE_WIDTH-offsetR;
            currQ = QueryFile.binReads[offsetQ?j:j-1];
            currQ >>= DATATYPE_WIDTH-offsetQ;
        } 

        if((currR & global_mask_right[matchSize/2 - 1]) != (currQ &  global_mask_right[matchSize/2 - 1])) {
            if (matchSize==2)
                break;

            mismatch=1;
            matchSize/=2;
            if (matchSize%2)
                matchSize+=1;
        }else {
            lRef-=matchSize;
            lQue-=matchSize;
            if (mismatch) {
                if (matchSize==2) 
                    break;
                currR >>= matchSize;
                currQ >>= matchSize;
            }
        }
    }
    
    if (totalRBits-lRef+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;
    
    if (totalQBits-lQue+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    //match towards right
    mismatch=0;
    while ((rRef <= totalRBits) && (rQue <= totalQBits) && (rRef <= RefNpos.right) && (rQue <= QueryNpos.right)) 
    {
        if (!mismatch)
        {
            offsetR=rRef%DATATYPE_WIDTH;
            i=rRef/DATATYPE_WIDTH;
            offsetQ=rQue%DATATYPE_WIDTH;
            j=rQue/DATATYPE_WIDTH;

            if (offsetR > offsetQ)
                matchSize = DATATYPE_WIDTH-offsetR;
            else
                matchSize = DATATYPE_WIDTH-offsetQ;
    
            if (rRef+matchSize > totalRBits)
                matchSize = totalRBits-rRef;
        
            if (rQue+matchSize > totalQBits)
                matchSize = totalQBits-rQue;

            if (rQue+matchSize > QueryNpos.right)
                matchSize = QueryNpos.right-rQue;
            
            if (rRef+matchSize > RefNpos.right)
                matchSize = RefNpos.right-rRef;

            if(!matchSize)
                matchSize=2;

            
        
            currR = RefFile.binReads[i];
            currR <<= offsetR;
            currQ = QueryFile.binReads[j];
            currQ <<= offsetQ;
        }
         
        if((currR & global_mask_left[matchSize/2 - 1]) != (currQ &  global_mask_left[matchSize/2 - 1])) {
            if (matchSize==2){
                rRef-=2;
                rQue-=2;
                break;
            }
            mismatch=1;
            matchSize/=2;
            if (matchSize%2)
                matchSize+=1;
        }else {
            if (mismatch) {
                if (matchSize==2)
                    break;
            }
            if ((rRef == totalRBits) || (rQue == totalQBits))
                break;
            
            currR <<= matchSize;
            currQ <<= matchSize;
            rRef+=matchSize;
            rQue+=matchSize;
        }
    }

    /* Adjust rRef and rQue locations */

    if (rRef > RefNpos.right){
        rQue-=(rRef-RefNpos.right);
        rRef=RefNpos.right;
    }
    if (rQue > QueryNpos.right){
        rRef-=(rQue-QueryNpos.right);
        rQue=QueryNpos.right;
    }

    if (rRef > totalRBits){
        rQue-=(rRef-totalRBits);
        rRef=totalRBits;
    }
    if (rQue > totalQBits){
        rRef-=(rQue-totalQBits);
        rQue=totalQBits;
    }

    if (arrayTmpFile.writeMemInTmpFiles(lRef, rRef, lQue, rQue, QueryFile, RefFile, revComplement)) {
        uint64_t key = ((lRef << 32) | rRef);
        uint64_t value = ((lQue << 32) | rQue);
        currMEMs.insert(std::make_pair(key, value));
        currQueryMEMs->ListAdd(&currQueryMEMs, lQue, rQue, key);
    }
}

void reportMEM(Knode * &refHash, uint64_t totalBases, uint64_t totalQBases, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, uint32_t &revComplement, TSwitch sw)
{
  uint64_t totalQBits = CHARS2BITS(totalQBases);
  uint32_t copyBits=0;

  #pragma omp parallel num_threads( sw . T )
  {
      queryList *currQueryMEMs = NULL;
      unordered_multimap <uint64_t, uint64_t> currMEMs;
      uint64_t currKmer=0, j=0;
      int32_t offset=0;
      uint32_t first=1;
      int kmerWithNs=0;
      mapObject QueryNpos, RefNpos;
      vector<mapObject>::iterator it;
      it = upper_bound(QueryFile.blockOfNs.begin(), QueryFile.blockOfNs.end(), 0, mapObject()); 

      #pragma omp single
      {
      /*
       * Number of copy bits during query kmer processing depends on kmer size.
       */
      if (DATATYPE_WIDTH-commonData::kmerSize > 32 )
          copyBits=32; //16 characters
      else if (DATATYPE_WIDTH-commonData::kmerSize > 16)
          copyBits=16; //8 characters
      else
          copyBits=8; //4 characters

      /* If copyBits more than 8, the for loop parallelisation will give 
       * incorrect results - miss some Mems
       */
      if( sw . T > 1)   
          copyBits=8; //4 characters
      }   
      

      #pragma omp for 
      for (uint64_t currKmerPos=0; currKmerPos<=totalQBits; currKmerPos+=2)
      {
          if ((currKmerPos + commonData::kmerSize - 2) > totalQBits)
              continue;
        
          if(QueryFile.checkKmerForNs(currKmerPos, it)){
              kmerWithNs=1;
          }

          j=currKmerPos/DATATYPE_WIDTH;// current location in binReads 
          offset = currKmerPos%DATATYPE_WIDTH;
          if(first || !offset){
              currKmer = QueryFile.binReads[j];
              currKmer <<= offset;
              if(offset > DATATYPE_WIDTH-commonData::kmerSize)
                  currKmer |= ((QueryFile.binReads[j+1] & global_mask_left[offset/2-1])>>(DATATYPE_WIDTH-offset));
              first=0;
          }else
              currKmer <<= 2;

          if(offset  && !(offset % copyBits))
              currKmer |= ((QueryFile.binReads[j+1] & global_mask_left[offset/2-1])>>(DATATYPE_WIDTH-offset));

          if (kmerWithNs){
             /* Do not process this Kmer, Ns in it */
             kmerWithNs=0;
             continue;
          }
          /* Find the K-mer in the refHash */
          uint64_t *dataPtr=NULL;
          if (refHash->findKmer(currKmer & global_mask_left[commonData::kmerSize/2 - 1], dataPtr)) 
          {
              // We have a match
              for (uint64_t n=1; n<=dataPtr[0]; n++) {   
                  // Check if MEM has already been discovered, if not proces it
                  if (!(currQueryMEMs->checkRedundantMEM(&currQueryMEMs, dataPtr[n], currKmerPos, CHARS2BITS(totalBases), currMEMs)))
                      helperReportMem(dataPtr[n], currKmerPos, CHARS2BITS(totalBases), CHARS2BITS(totalQBases), currQueryMEMs, currMEMs, RefFile, QueryFile, arrayTmpFile, RefNpos, QueryNpos, revComplement);
              }
          }
      }
      currMEMs.clear();
      currQueryMEMs->ListFree(&currQueryMEMs);
  }  
}

void processQuery(Knode * &refHash, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, uint32_t &revComplement, TSwitch sw)
{
    QueryFile.clearFileFlag();
    QueryFile.resetCurrPos();
    for (int32_t i=0; i<commonData::d; i++) {
        if(QueryFile.readChunks()){
            reportMEM(refHash, RefFile.totalBases-1, QueryFile.totalBases-1, RefFile, QueryFile, arrayTmpFile, revComplement, sw);
            QueryFile.setCurrPos();
            QueryFile.clearMapForNs();
        }
        else
            break;
    }
    QueryFile.clearTmpString();
}

void processReference(seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, uint32_t &revComplement , TSwitch sw)
{
    uint64_t numberOfKmers=0,n=0;
    int hashTableSizeIndex=0;
    Knode *refHash;
    
    numberOfKmers = ceil((RefFile.totalBases-commonData::kmerSize/2+1)/((commonData::minMemLen/2-commonData::kmerSize/2 + 1)) + 1);

    /* Set the size of the hash table to the numberofKmers. */
    for (n=0; n<450; ++n)
    {
        if (hashTableSize[n] > 1.75*numberOfKmers)
        {
           hashTableSizeIndex = n;
           break;
        }
    }

    Knode::currHashTabSize = hashTableSize[hashTableSizeIndex];  //Store the size of the hash table.
    if (hashTableSizeIndex)
        Knode::prevHashTabSize = hashTableSize[hashTableSizeIndex-1];
    else
        Knode::prevHashTabSize = 3;

    /* Create the refHash for K-mers. */
    refHash = new Knode[Knode::currHashTabSize];

    buildRefHash(refHash, CHARS2BITS(RefFile.totalBases-1), RefFile);

    processQuery(refHash, RefFile, QueryFile, arrayTmpFile, revComplement, sw );

    delete [] refHash; 
}

bool is_numeric(const string &str)
{
    return all_of(str.begin(), str.end(), ::isdigit); 
}


int find_maximal_exact_matches( unsigned int l, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, TSwitch sw )
{

    fprintf ( stderr, " -Identifying maximal exact matches of minimum length %i\n", l );
    
    int32_t i=0, n=1;
    uint32_t options, revComplement=0;
    seqFileReadInfo RefFile, QueryFile;

    
    RefFile.openFile( string(sw.output_filename)+"_new_ref.fa" );
    QueryFile.openFile( string(sw.output_filename)+"_new_query.fa" );

    commonData::minMemLen = 2* l;
    if( l % 2 == 0 )
    	commonData::kmerSize = l; 
    else commonData::kmerSize = l + 1;

    sprintf(commonData::nucmer_path, "%s/%d_tmp", getenv("NUCMER_E_MEM_OUTPUT_DIRPATH")?getenv("NUCMER_E_MEM_OUTPUT_DIRPATH"):".",getpid());

    tmpFilesInfo arrayTmpFile(IS_MATCH_BOTH_DEF(options)?(2*NUM_TMP_FILES+2):NUM_TMP_FILES+2);
    arrayTmpFile.openFiles(ios::out|ios::binary, IS_MATCH_BOTH_DEF(options)?(2*NUM_TMP_FILES+2):NUM_TMP_FILES+2);

    RefFile.generateRevComplement(0); // This routine also computers size and num sequences
    QueryFile.generateRevComplement(0); // Reverse complement only for query

    arrayTmpFile.setNumMemsInFile(QueryFile.allocBinArray(), QueryFile.getNumSequences());
    RefFile.allocBinArray();
    RefFile.clearFileFlag();

    while (true)
    {
        for (i=0; i<commonData::d; i++) {
            if(RefFile.readChunks()){
                processReference(RefFile, QueryFile, arrayTmpFile, revComplement, sw);
                RefFile.setCurrPos();
                RefFile.clearMapForNs();
            }
            else
                break;
        }

        /*
         * Process MemExt list 
         */ 

        arrayTmpFile.mergeMemExtVector(revComplement);

        if (revComplement)
            break;
        if (IS_MATCH_BOTH_DEF(options)){
            SET_MATCH_BOTH(revComplement);
            //revComplement=1;
            RefFile.clearFileFlag();
            RefFile.resetCurrPos();
            RefFile.totalBases=0;
            QueryFile.setReverseFile();
            QueryFile.totalBases=0;
        }
        else
            break;
    }

    /*
     * Free up the allocated arrays
     */
    arrayTmpFile.closeFiles(IS_MATCH_BOTH_DEF(options)?(2*NUM_TMP_FILES):NUM_TMP_FILES);
    RefFile.destroy();
    QueryFile.destroy();

    /* 
     * Populate sequence information in vectors. Use this to get MEM
     * positions relative to the original sequences.
     */
    vector<seqData> refSeqInfo;
    vector<seqData> querySeqInfo;
    refSeqInfo.reserve(RefFile.getNumSequences());
    querySeqInfo.reserve(QueryFile.getNumSequences());
    RefFile.generateSeqPos(refSeqInfo);
    QueryFile.generateSeqPos(querySeqInfo);
    RefFile.closeFile();
    QueryFile.closeFile();

    arrayTmpFile.removeDuplicates(refSeqInfo, querySeqInfo, revComplement, q_grams, l);


    return 0;
}

