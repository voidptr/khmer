#ifndef BIT_ARRAY_H
#define BIT_ARRAY_H

#include <assert.h>

#include <iostream>

#include "tArray.h"

using namespace std;

// Class: cBitArray and cBitMatrix
// Desc: These classes handle an arbitrarily large array or matrix of bits,
//       and optimizes the operations on those bits to be as fast as possible.
//

// Constructors:
//  cBitArray()                            -- Assume a size zero array.
//  cBitArray(int in_size)                 -- Create an uninitialized array.
//  cBitArray(const cBitArray & in_array)  -- Copy Constructor

// Assignment and equality test:
//  cBitArray & operator=(const cBitArray & in_array)
//  bool operator==(const cBitArray & in_array) const

// Sizing:
//  int GetSize() const
//  void Resize(const int new_size)
//  void ResizeClear(const int new_size)

// Accessors:
//  void Set(int index, bool value)
//  bool Get(int index) const
//  bool operator[](int index) const
//  cBitProxy operator[](int index)
//  void Clear()
//  void SetAll()

// Printing:
//  void Print(ostream & out=cout) const
//  void PrintOneIDs(ostream & out=cout) const

// Bit play:
//  int CountBits()   -- Count 1s -- fast for sparse arrays.
//  int CountBits2()  -- Count 1s -- fast for arbitary arrays.
//  int FindBit1(int start_bit)   -- Return pos of first 1 after start_bit 

// Boolean math functions:
//  cBitArray NOT() const
//  cBitArray AND(const cBitArray & array2) const
//  cBitArray OR(const cBitArray & array2) const
//  cBitArray NAND(const cBitArray & array2) const
//  cBitArray NOR(const cBitArray & array2) const
//  cBitArray XOR(const cBitArray & array2) const
//  cBitArray EQU(const cBitArray & array2) const
//  cBitArray SHIFT(const int shift_size) const   -- positive for left shift, negative for right shift

//  const cBitArray & NOTSELF()
//  const cBitArray & ANDSELF(const cBitArray & array2)
//  const cBitArray & ORSELF(const cBitArray & array2)
//  const cBitArray & NANDSELF(const cBitArray & array2)
//  const cBitArray & NORSELF(const cBitArray & array2)
//  const cBitArray & XORSELF(const cBitArray & array2)
//  const cBitArray & EQUSELF(const cBitArray & array2)
//  const cBitArray & SHIFTSELF(const int shift_size) const

// Arithmetic:
//  cBitArray INCREMENTSELF()

// Operator overloads:
//  cBitArray operator~() const
//  cBitArray operator&(const cBitArray & ar2) const
//  cBitArray operator|(const cBitArray & ar2) const
//  cBitArray operator^(const cBitArray & ar2) const
//  cBitArray operator>>(const int) const
//  cBitArray operator<<(const int) const
//  const cBitArray & operator&=(const cBitArray & ar2)
//  const cBitArray & operator|=(const cBitArray & ar2)
//  const cBitArray & operator^=(const cBitArray & ar2)
//  const cBitArray & operator>>=(const int)
//  const cBitArray & operator<<=(const int)
//  cBitArray & operator++()     // prefix ++
//  cBitArray & operator++(int)  // postfix ++

// The following is an internal class used by cBitArray (and will be used in
// cBitMatrix eventually....).  It does not keep track of size, so this value
// must be passed in.

class cRawBitArray {
private:
  unsigned int * bit_fields;
  
  static const unsigned int position_masks[32];
  
  // Disallow default copy constructor and operator=
  // (we need to know the number of bits we're working with!)
  cRawBitArray(const cRawBitArray & ) { assert(false); }
  const cRawBitArray & operator=(const cRawBitArray & in_array)
  { assert(false); return *this; }
  
  inline unsigned long long GetNumFields(const unsigned long long num_bits) const { return 1 + ((num_bits - 1) >> 5); }
  inline unsigned long long GetField(const unsigned long long index) const { return index >> 5; }
  inline int GetFieldPos(const unsigned long long index) const { return index & 31; }
public:
  cRawBitArray() : bit_fields(NULL) 
  { 
    ;
  }
  ~cRawBitArray() {
    if (bit_fields != NULL) {
      delete [] bit_fields;
    }
  }
  
  void Zero(const unsigned long long num_bits) {
    const int num_fields = GetNumFields(num_bits);
    for (unsigned long long i = 0; i < num_fields; i++) {
      bit_fields[i] = 0;
    }    
  }
  
  void Ones(const unsigned long long num_bits) {
    const unsigned long long num_fields = GetNumFields(num_bits);
    for (unsigned long long i = 0; i < num_fields; i++) {
      bit_fields[i] = ~0;
    }    
    const int last_bit = GetFieldPos(num_bits);
    if (last_bit > 0) {
      bit_fields[num_fields - 1] &= (1 << last_bit) - 1;
    }
  }
  
  cRawBitArray(const unsigned long long num_bits) {
    const unsigned long long num_fields = GetNumFields(num_bits);
    bit_fields = new unsigned int[ num_fields ];
    Zero(num_bits);
  }
  
  // The Copy() method and the Copy Constructor must both be told how many
  // bits they are working with.
  void Copy(const cRawBitArray & in_array, const unsigned long long num_bits);
  cRawBitArray(const cRawBitArray & in_array, const unsigned long long num_bits)
  : bit_fields(NULL)
  {
    Copy(in_array, num_bits);
  }
  
  // For fast bit operations, we're not going to setup operator[]; instead
  // we're going to have a GetBit and a SetBit commamd.  For this raw version
  // we're also going to assume that the index is within range w/o any special
  // checks.
  bool GetBit(const unsigned long long index) const{
    const unsigned long long field_id = GetField(index);
    const int pos_id = GetFieldPos(index);
    
 //   unsigned int lVal = bit_fields[field_id];
  
 //       return (lVal & (1 << pos_id)) != 0;
    return (bit_fields[field_id] & (1 << pos_id)) != 0;
    //return (bit_fields[field_id] & (position_masks[pos_id])) != 0;
  }
  
  void SetBit(const unsigned long long index, const bool value) {
    const unsigned long long field_id = GetField(index);
    const int pos_id = GetFieldPos(index);
    const unsigned int pos_mask = position_masks[pos_id];
    
    if ( value )
      bit_fields[field_id] |= pos_mask;
    else
      bit_fields[field_id] &= ~pos_mask;
    
  }
  
  void IncrementSemiNibble(const unsigned long long index)
  {
    const unsigned long long lSemiNibbleIndex = index * 2;
//    assert( (index & 1) == 0 ); // is even.
    
    const unsigned long long field_id = GetField(lSemiNibbleIndex);
    const int pos_id = GetFieldPos(index);
    const unsigned int pos_mask = position_masks[pos_id];

//    
//    unsigned int lVal = bit_fields[field_id];
    
    bit_fields[field_id] |= ( (bit_fields[field_id] << 1) & ( pos_mask << 1 ) ) | ( bit_fields[field_id] | pos_mask );
    
//    if ( value )
//      bit_fields[field_id] |= pos_mask;
//    else
//      bit_fields[field_id] &= ~pos_mask;
  }
  
  bool IsEqual(const cRawBitArray & in_array, unsigned long long num_bits) const;
  
  void Resize(const unsigned long long old_bits, const unsigned long long new_bits);
  void ResizeSloppy(const unsigned long long new_bits);
  void ResizeClear(const unsigned long long new_bits);
  
  // Two different technique of bit counting...
  unsigned long long CountBits(const unsigned long long num_bits) const; // Better for sparse arrays
  unsigned long long CountBits(const unsigned long long start_bit, const unsigned long long stop_bit) const; // count the bits between a start and stop point
  unsigned long long CountBits2(const unsigned long long num_bits) const; // Better for dense arrays
  
  // Other bit-play
  unsigned long long FindBit1(const unsigned long long num_bits, const unsigned long long start_pos) const;
  tArray<int> GetOnes(const unsigned long long num_bits) const;
  void ShiftLeft(const unsigned long long num_bits, const unsigned long long shift_size); // Helper: call SHIFT with positive number instead
  void ShiftRight(const unsigned long long num_bits, const unsigned long long shift_size); // Helper: call SHIFT with negative number instead
  
  void Print(const unsigned long long num_bits, ostream & out=cout) const {
    for (unsigned long long i = 0; i < num_bits; i++) {
      out << GetBit(i);
    }
  }
  
  // prints in the accepted human readable low-to-hight = right-to-left format, taking bit 0 as low bit
  void PrintRightToLeft(const unsigned long long num_bits, ostream & out=cout) const {
    for (unsigned long long i = num_bits - 1; i >= 0; i--) {
      out << GetBit(i);
    }
  }
  
  void PrintOneIDs(const unsigned long long num_bits, ostream & out=cout) const {
    for (unsigned long long i = 0; i < num_bits; i++) {
      if (GetBit(i) == true) {
        out << i << " ";
      }
    }
  }
  
  // Fast bool operators where we uses this bit array as one of the 
  // inputs and the place to store the results.
  void NOT(const unsigned long long num_bits);
  void AND(const cRawBitArray & array2, const unsigned long long num_bits);
  void OR(const cRawBitArray & array2, const unsigned long long num_bits);
  void NAND(const cRawBitArray & array2, const unsigned long long num_bits);
  void NOR(const cRawBitArray & array2, const unsigned long long num_bits);
  void XOR(const cRawBitArray & array2, const unsigned long long num_bits);
  void EQU(const cRawBitArray & array2, const unsigned long long num_bits);
  void SHIFT(const unsigned long long num_bits, const unsigned long long shift_size);  // positive numbers for left and negative for right (0 does nothing)
  void INCREMENT(const unsigned long long num_bits);
  
  // Fast bool operators where we load all of the inputs and store the
  // results here.
  void NOT(const cRawBitArray & array1, const unsigned long long num_bits);
  void AND(const cRawBitArray & array1, const cRawBitArray & array2,
           const unsigned long long num_bits);
  void OR(const cRawBitArray & array1, const cRawBitArray & array2,
          const unsigned long long num_bits);
  void NAND(const cRawBitArray & array1, const cRawBitArray & array2,
            const unsigned long long num_bits);
  void NOR(const cRawBitArray & array1, const cRawBitArray & array2,
           const unsigned long long num_bits);
  void XOR(const cRawBitArray & array1, const cRawBitArray & array2,
           const unsigned long long num_bits);
  void EQU(const cRawBitArray & array1, const cRawBitArray & array2,
           const unsigned long long num_bits);
  void SHIFT(const cRawBitArray & array1, const unsigned long long num_bits, const unsigned long long shift_size);
  void INCREMENT(const cRawBitArray & array1, const unsigned long long num_bits);  // implemented for completeness, but unused by cBitArray
};


class cBitArray {
private:
  cRawBitArray bit_array;
  unsigned long long array_size;

  // Setup a bit proxy so that we can use operator[] on bit arrays as a lvalue.
  class cBitProxy {
  private:
    cBitArray & array;
    unsigned long long index;
  public:
    cBitProxy(cBitArray & _array, unsigned long long _idx) : array(_array), index(_idx) {;}

    inline cBitProxy & operator=(bool b);    // lvalue handling...
    inline operator bool() const;            // rvalue handling...
  };
  friend class cBitProxy;
public:
  cBitArray() : array_size(0) { ; }
  cBitArray(unsigned long long in_size) : bit_array(in_size), array_size(in_size) { ; }
  cBitArray(const cBitArray & in_array) : bit_array(in_array.bit_array, in_array.array_size), array_size(in_array.array_size) { ; }
  cBitArray(const cRawBitArray & in_array, int in_size) : bit_array(in_array, in_size), array_size(in_size) { ; }
  ~cBitArray()
  {
    ;
  }
  
  cBitArray & operator=(const cBitArray & in_array) {
    bit_array.Copy(in_array.bit_array, in_array.array_size);
    array_size = in_array.array_size;
    return *this;
  }

  bool operator==(const cBitArray & in_array) const {
    if (array_size != in_array.array_size) return false;
    return bit_array.IsEqual(in_array.bit_array, array_size);
  }

  unsigned long long GetSize() const { return array_size; }

  void Set(unsigned long long index, bool value) {
    assert(index < array_size);
    bit_array.SetBit(index, value);
  }
  
  void IncrementSemiNibble(unsigned long long index)
  {
    assert(index*2 < array_size);
    bit_array.IncrementSemiNibble(index);
  }

  bool Get(unsigned long long index) const {
    assert(index < array_size);
    return bit_array.GetBit(index);
  }

  bool operator[](unsigned long long index) const { return Get(index); }
  cBitProxy operator[](unsigned long long index) { return cBitProxy(*this, index); }

  void Clear() { bit_array.Zero(array_size); }
  void SetAll() { bit_array.Ones(array_size); }
  

  void Print(ostream & out=cout) const { bit_array.Print(array_size, out); }
  void PrintRightToLeft(ostream & out=cout) const { bit_array.PrintRightToLeft(array_size, out); }
  void PrintOneIDs(ostream & out=cout) const
    { bit_array.PrintOneIDs(array_size, out); }
  void Resize(const unsigned long long new_size) {
    bit_array.Resize(array_size, new_size);
    array_size = new_size;
  }
  void ResizeClear(const unsigned long long new_size) {
    bit_array.ResizeClear(new_size);
    array_size = new_size;
  }
  unsigned long long CountBits() const { return bit_array.CountBits(array_size); }
  unsigned long long CountBits( unsigned long long start_bit, unsigned long long stop_bit ) const { return bit_array.CountBits( start_bit, stop_bit ); }
  unsigned long long CountBits2() const { return bit_array.CountBits2(array_size); }

  unsigned long long FindBit1(unsigned long long start_bit=0) const
    { return bit_array.FindBit1(array_size, start_bit); }
  tArray<int> GetOnes() const { return bit_array.GetOnes(array_size); }

  // Boolean math functions...
  cBitArray NOT() const {
    cBitArray out_array;
    out_array.bit_array.NOT(bit_array, array_size);
    out_array.array_size = array_size;
    return out_array;
  }

  cBitArray AND(const cBitArray & array2) const {
    assert(array_size == array2.array_size);
    cBitArray out_array;
    out_array.bit_array.AND(bit_array, array2.bit_array, array_size);
    out_array.array_size = array_size;
    return out_array;
  }

  cBitArray OR(const cBitArray & array2) const {
    assert(array_size == array2.array_size);
    cBitArray out_array;
    out_array.bit_array.OR(bit_array, array2.bit_array, array_size);
    out_array.array_size = array_size;
    return out_array;
  }

  cBitArray NAND(const cBitArray & array2) const {
    assert(array_size == array2.array_size);
    cBitArray out_array;
    out_array.bit_array.NAND(bit_array, array2.bit_array, array_size);
    out_array.array_size = array_size;
    return out_array;
  }

  cBitArray NOR(const cBitArray & array2) const {
    assert(array_size == array2.array_size);
    cBitArray out_array;
    out_array.bit_array.NOR(bit_array, array2.bit_array, array_size);
    out_array.array_size = array_size;
    return out_array;
  }

  cBitArray XOR(const cBitArray & array2) const {
    assert(array_size == array2.array_size);
    cBitArray out_array;
    out_array.bit_array.XOR(bit_array, array2.bit_array, array_size);
    out_array.array_size = array_size;
    return out_array;
  }

  cBitArray EQU(const cBitArray & array2) const {
    assert(array_size == array2.array_size);
    cBitArray out_array;
    out_array.bit_array.EQU(bit_array, array2.bit_array, array_size);
    out_array.array_size = array_size;
    return out_array;
  }
  
  cBitArray SHIFT(const int shift_size) const {
    cBitArray out_array;
    out_array.bit_array.SHIFT(bit_array, array_size, shift_size);
    out_array.array_size = array_size;
    return out_array;
  }

  const cBitArray & NOTSELF() {
    bit_array.NOT(array_size);
    return *this;
  }

  const cBitArray & ANDSELF(const cBitArray & array2) {
    assert(array_size == array2.array_size);
    bit_array.AND(array2.bit_array, array_size);
    return *this;
  }

  const cBitArray & ORSELF(const cBitArray & array2) {
    assert(array_size == array2.array_size);
    bit_array.OR(array2.bit_array, array_size);
    return *this;
  }

  const cBitArray & NANDSELF(const cBitArray & array2) {
    assert(array_size == array2.array_size);
    bit_array.NAND(array2.bit_array, array_size);
    return *this;
  }

  const cBitArray & NORSELF(const cBitArray & array2) {
    assert(array_size == array2.array_size);
    bit_array.NOR(array2.bit_array, array_size);
    return *this;
  }

  const cBitArray & XORSELF(const cBitArray & array2) {
    assert(array_size == array2.array_size);
    bit_array.XOR(array2.bit_array, array_size);
    return *this;
  }

  const cBitArray & EQUSELF(const cBitArray & array2) {
    assert(array_size == array2.array_size);
    bit_array.EQU(array2.bit_array, array_size);
    return *this;
  }
  
  const cBitArray & SHIFTSELF(const int shift_size) {
    bit_array.SHIFT(array_size, shift_size);
    return *this;
  }
  
  cBitArray & INCREMENTSELF() {
    bit_array.INCREMENT(array_size);
    return *this;
  }
  

  // Operator overloads...
  cBitArray operator~() const { return NOT(); }
  cBitArray operator&(const cBitArray & ar2) const { return AND(ar2); }
  cBitArray operator|(const cBitArray & ar2) const { return OR(ar2); }
  cBitArray operator^(const cBitArray & ar2) const { return XOR(ar2); }
  cBitArray operator<<(const int shift_size) const { return SHIFT(shift_size); }
  cBitArray operator>>(const int shift_size) const { return SHIFT(-shift_size); }
  const cBitArray & operator&=(const cBitArray & ar2) { return ANDSELF(ar2); }
  const cBitArray & operator|=(const cBitArray & ar2) { return ORSELF(ar2); }
  const cBitArray & operator^=(const cBitArray & ar2) { return XORSELF(ar2); }
  const cBitArray & operator<<=(const int shift_size) { return SHIFTSELF(shift_size); }
  const cBitArray & operator>>=(const int shift_size) { return SHIFTSELF(-shift_size); }
  cBitArray & operator++() { return INCREMENTSELF(); }  // prefix ++
  cBitArray operator++(int) { cBitArray ans = *this; operator++(); return ans;}  // postfix ++

};

std::ostream & operator << (std::ostream & out, const cBitArray & bit_array);

cBitArray::cBitProxy & cBitArray::cBitProxy::operator=(bool b)
{
  array.Set(index, b);
  return *this;
}


cBitArray::cBitProxy::operator bool() const
{
  return array.Get(index);
}

#endif
