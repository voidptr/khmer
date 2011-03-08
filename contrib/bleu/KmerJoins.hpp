/*
 *  KmerJoins.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 3/7/11.
 *  Copyright 2011 VoidPtr. All rights reserved.
 *
 */
 
#include <sstream>
#include <iomanip>

//#include <iostream>

namespace bleu {
  
  using namespace khmer;
  using namespace std;
    
  class KmerJoins
  {
  private:
    vector<string> mReads;
    vector<string> mReadNames;
    int mK;
    string mKmer;
    int mNumberOfReads;
    int mMaxAlignmentLength;
    int mJoinCount;
    bool mIsAligned;
    bool mPairwiseDifferencesCalculated;
    bool mSegregatingSitesCalculated;
    bool mA1Calculated;
    
    int mTotalPairwiseDifferences;
    int mTotalPairwiseComparisons;
    vector< vector<int> > mPairwiseDifferences;
    vector< vector<int> > mPairwiseComparisons;

    int mNumberOfSegregatingSites;
    int mNumberOfComparableSegregatingSites;
    
    double mA1;
    
  public:
    KmerJoins( const string & aKmer, 
               const vector<string> & aReads, 
               const vector<string> & aReadNames, 
               const int aJoinCount ) :
               mReadNames( aReadNames ), 
               mReads( aReads ),
               mKmer( aKmer ),
               mK( aKmer.length() ),
               mNumberOfReads( aReads.size() ),
               mJoinCount( aJoinCount ),
               mIsAligned( false ),
               mPairwiseDifferencesCalculated( false ),
               mSegregatingSitesCalculated( false ),
               mA1Calculated( false )
    {
    }
    
  public:
  
    double Pi()
    {
      CalculatePairwiseDifferences();
      
      if ( mTotalPairwiseComparisons == 0 )
        return 0.0;
      
      return mTotalPairwiseDifferences / (double) mTotalPairwiseComparisons;
    }
    
    double S()
    {
      CalculateSegregatingSites();
      
      if ( mNumberOfComparableSegregatingSites == 0 )
        return 0.0;
      
      return mNumberOfSegregatingSites / (double) mNumberOfComparableSegregatingSites;
    }
    
    double A1()
    {
      CalculateA1();
      return mA1;
    }
    
    double E_Theta()
    {
      return ( S() / A1() );
    }
    
    string OutputKmerInfo()
    {
      stringstream lOut;
      
      lOut << OutputAlignment(true);

      lOut << "S: " << S() << endl;            
      lOut << "Pi: " << Pi() << endl;
      lOut << "E(Theta): " << E_Theta() << endl;
      
      return lOut.str();
    }
    
  private:
    void AlignReads()
    {
      if ( mIsAligned )
        return;
        
      // first, divine the alignment
      int lLongestAlignment = 0;


      for ( int i = 0; i < mNumberOfReads; ++i )
      {
        int lAlignmentStart = mReads[i].find( mKmer ); 
        
        if ( lAlignmentStart > lLongestAlignment )
          lLongestAlignment = lAlignmentStart;
      }

      mMaxAlignmentLength = 0;
      for ( int j = 0; j < mNumberOfReads; ++j )
      {  
        int leftpad = ( lLongestAlignment - mReads[j].find( mKmer ) );
        mReads[j].insert( 0, leftpad, ' ' );

        if ( mReads[j].length() > mMaxAlignmentLength )
          mMaxAlignmentLength = mReads[j].length();
      }
      
      mIsAligned = true;
    }
    
    void CalculatePairwiseDifferences()
    {
      if ( mPairwiseDifferencesCalculated )
        return;
      
      AlignReads();
            
      mPairwiseDifferences.resize( mNumberOfReads );
      mPairwiseComparisons.resize( mNumberOfReads );
      for ( int i = 0; i < mNumberOfReads; ++i )
      {
        mPairwiseDifferences[i].resize( mNumberOfReads );
        mPairwiseComparisons[i].resize( mNumberOfReads );
      }

      mTotalPairwiseDifferences = 0;
      mTotalPairwiseComparisons = 0;
      
      for (int k = 0; k < mNumberOfReads - 1; ++k)
      {                
        for (int j = k + 1; j < mNumberOfReads; ++j)
        {
          pair<int, int> lPairwiseDiff = Calculate_PairwiseDifferences( mReads[k], mReads[j] );
          
          mPairwiseComparisons[k][j] = lPairwiseDiff.first;
          mPairwiseDifferences[k][j] = lPairwiseDiff.second;
          
          mTotalPairwiseComparisons += lPairwiseDiff.first;
          mTotalPairwiseDifferences += lPairwiseDiff.second;
          
        }
      }
        
      mPairwiseDifferencesCalculated = true;
    }
    
    void CalculateSegregatingSites()
    {
      if ( mSegregatingSitesCalculated )
        return;
        
      AlignReads();
        
      mNumberOfSegregatingSites = 0;
      
      for ( int i = 0; i < mMaxAlignmentLength; ++i ) // go through each nucleotide, left to right
      {
        set<char> lNucleotides;
        for ( int j = 0; j < mNumberOfReads; ++j ) // go through the aligned reads, top to bottom.
        {
          if ( i < mReads[j].length() )
          {
            char lNucleotide = mReads[j][i];
            if ( lNucleotide != ' ' )
              lNucleotides.insert( lNucleotide );
          }
        }
        
        if ( lNucleotides.size() > 1 )
        {
          mNumberOfSegregatingSites++;
        }
      }
      
      mNumberOfComparableSegregatingSites = mMaxAlignmentLength - mK;
      
//      if ( aComparableSitesCount == 0 )
//        return 0.0;
//        
//      return (lSegregatingSitesCount / (double) aComparableSitesCount);
         
      mSegregatingSitesCalculated = true;
    }

    void CalculateA1()
    {
      if ( mA1Calculated )
        return;
    
      mA1 = 0;
      for ( int i = 1; i < mNumberOfReads; ++i )
      {
        mA1 += ( 1 / (double) i );
      }
      
      mA1Calculated = true;
    }
    
    string OutputAlignment(bool aWithNames=false)
    {
      stringstream lOut;
      
      int lNamePad = 0;
      if ( aWithNames )
      {
        for (int i = 0; i < mNumberOfReads; ++i)
        {
          if ( mReadNames[i].length() > lNamePad )
            lNamePad = mReadNames[i].length();
        }
        lNamePad += 2;
      }      
      
      if ( aWithNames )
        lOut << setw(lNamePad) << " ";
      lOut << GenerateAlignedKmer() << endl;
      for (int i = 0; i < mNumberOfReads; ++i )
      {
        if ( aWithNames )
          lOut << setw(lNamePad) << mReadNames[i] << " ";
        lOut << mReads[i] << endl;
      }
      if ( aWithNames )
        lOut << setw( lNamePad ) << " ";
      lOut << GenerateMismatchMarks() << endl;
      
      return lOut.str();
    }
    
    string GenerateMismatchMarks()
    {
      AlignReads();
    
      string lMismatchMarks;
    
      for ( int i = 0; i < mMaxAlignmentLength; ++i ) // go through each nucleotide, left to right
      {
        bool lMismatch = false;
        char lConsensusNuc = ' ';
        for ( int j = 0; j < mNumberOfReads; ++j ) // go through the aligned reads, top to bottom.
        {
          if ( i < mReads[j].length() )
          {
            char lNucleotide = mReads[j][i];
            if ( lNucleotide != ' ' )
            {
              if ( lConsensusNuc != ' ' && lNucleotide != lConsensusNuc )
              {
                lMismatch = true;
                break;
              }
              lConsensusNuc = lNucleotide;
            }
          }
        }
        
        if ( lMismatch )
          lMismatchMarks.append("#");
        else
          lMismatchMarks.append(" ");
      }
      
      return lMismatchMarks;
    }
    
    string GenerateAlignedKmer()
    {
      AlignReads();
    
      string lAlignedKmer = mKmer;
    
      // just grab the first one. the others should already be lined up.
      int lBareAlignedLocation = mReads[0].find( mKmer ); 
      
      lAlignedKmer.insert(0, lBareAlignedLocation, ' ');
      
      return lAlignedKmer;
    }

    
  private: /// really really private
    pair<int,int> Calculate_PairwiseDifferences( const string & aRead1, const string & aRead2 )
    {
      int lComparableSites = 0;
      int lDifferences = 0;
      
      int lCompareLength = aRead1.length();
      if ( aRead1.length() > aRead2.length() ) 
        lCompareLength = aRead2.length();
        
      for ( int i = 0; i < lCompareLength; ++i )
      {
        if ( aRead1[i] != ' ' && aRead2[i] != ' ')
        {
          lComparableSites++;
          if ( aRead1[i] != aRead2[i] )
            lDifferences++;
        }
      }
      
      lComparableSites -= mK; // we're only interested in the non-matched sections.
      
      return pair<int,int>(lComparableSites, lDifferences);
    }
  };
}