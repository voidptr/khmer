///*
// *  BitArray.h
// *  bleu
// *
// */
//
//namespace bleu {
//  
//
//  class BitArray {
//  private:
//    unsigned long long * bit_fields;
//    unsigned long long num_bits;
//    unsigned long long num_cells;
//    bool SemiNibbleCells;
//    unsigned int multiplier;
//    unsigned int mask;
//    unsigned int num_fields;
//  public:
//  //  BitArray() : bit_fields(NULL) 
//  //  { 
//  //    SemiNibbleCells = false;
//  //    ;
//  //  }
//
//    BitArray(const unsigned long long aNumCells, bool aAsSemiNibble)
//    {
//      SemiNibbleCells = aAsSemiNibble;
//    
//      num_bits = aNumCells;
//      num_cells = aNumCells;
//      
//      if ( aAsSemiNibble ) // only 32 addressible positions.
//      {
//        multiplier = 5;
//        mask = 31;
//        num_bits *= 2;
//      }
//      else // 64 addressible positions.
//      {
//        multiplier = 6;
//        mask = 64;
//      }
//      
//      num_fields = 1 + ((num_bits - 1) >> multiplier);
//      
//      bit_fields = new unsigned long long[ num_fields ];
//      memset( bit_fields, 0, (sizeof(unsigned long long) * num_fields) );
//    }
//  private:
//    ~BitArray() {
//      if (bit_fields != NULL) {
//        delete [] bit_fields;
//      }
//    }
//    inline unsigned long long GetField(const unsigned long long index) const { return (index >> multiplier); }
//    inline int GetFieldPos(const unsigned long long index) const { return index & mask; }
//  //  inline int GetFieldPosBit2(const unsigned long long index) const { return (index & (mask << 32)); };
//  public:
//    unsigned char GetBit( const unsigned long long aIndex ) const
//    {
//      assert (!SemiNibbleCells);
//      assert( aIndex < num_cells );
//    
//      const unsigned long long field_id = GetField(aIndex);
//      const int pos_id = GetFieldPos(aIndex);
//      
//      return (bit_fields[field_id] & (1 << pos_id)) != 0;
//    }
//    
//    void IncrementSemiNibble( const unsigned long long aIndex )
//    {
//      const unsigned long long field_id = GetField(aIndex);
//      const unsigned int field_pos1 = GetFieldPos(aIndex);
//
//      assert (SemiNibbleCells);
//        
//      bit_fields[field_id] |= ((bit_fields[field_id] & (1 << field_pos1)) << field_pos1) | (1 << field_pos1);          
//    }
//    
//    void CollapseArrayToSecondBit()
//    {
//      assert (SemiNibbleCells);
//        
//      unsigned long long * new_bit_fields = new unsigned long long[num_fields/2]; //only need half as much space now.
//        
//      for( unsigned long long new_field_index = 0, old_field_index = 0; old_field_index < num_fields; ++new_field_index, old_field_index+=2 )
//      {      
//        new_bit_fields[new_field_index] = (bit_fields[old_field_index] >>32) | (bit_fields[old_field_index+1] & ~4294967295ULL);
//      }
//      
//      delete bit_fields;
//      bit_fields = new_bit_fields;
//      
//      num_bits = num_cells;
//      multiplier = 6;
//      num_fields = 1 + ((num_bits - 1) >> multiplier);
//      SemiNibbleCells = false;
//      
//      //return true;
//    }
//    
//    // This technique counts the number of bits; it loops through once for each
//    // bit equal to 1.  This is reasonably fast for sparse arrays.
//    unsigned long long CountBits() const
//    {
//  //    const unsigned long long num_fields = GetNumFields(num_bits);
//      unsigned long long bit_count = 0;
//      
//      for (unsigned long long i = 0; i < num_fields; i++) 
//      {
//        unsigned int temp = bit_fields[i];
//        while (temp != 0) {
//          temp = temp & (temp - 1);
//          bit_count++;
//        }
//      }
//      return bit_count;
//    }
//
//    // This technique counts the number of bits; it loops through once for each
//    // bit equal to 1.  This is reasonably fast for sparse arrays.
//    unsigned long long CountBits(const unsigned long long start_bit, const unsigned long long stop_bit ) const
//    {
//      //const int num_fields = GetNumFields(stop_bit - start_bit);
//      unsigned long long bit_count = 0;
//      
//      unsigned long long start_field = GetField( start_bit );
//      unsigned long long stop_field = GetField( stop_bit );
//      unsigned int start_pos = GetFieldPos( start_bit );
//      unsigned int stop_pos = GetFieldPos( stop_bit );
//      
//      for (unsigned long long i = start_field; i <= stop_field; ++i) {
//        unsigned int temp = bit_fields[i];
//        
//        //cout << i << " " << start_field << endl;
//        
//        if ( i == start_field ) {
//          //cout << "SHIFTING RIGHT" << temp << " BY " << start_pos;
//          temp = temp >> start_pos;
//          //cout << " TO " << temp << endl;
//        }
//        
//        if ( i == stop_field ) {      
//          //cout << "SHIFTING LEFT" << temp << " BY " << (sizeof(unsigned int) * 8) - 1 - stop_pos;
//          temp = temp << (sizeof(unsigned long long) * 8) - 1 - stop_pos;
//          //cout << " TO " << temp << endl;
//        }
//
//        
//        while (temp != 0) {
//          temp = temp & (temp - 1);
//          bit_count++;
//        }
//      }
//      return bit_count;
//    }
//    
//  }; 
//}  