/*
 *  ReadUtilities.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 3/10/11.
 *
 */
 
#define is_valid_dna(ch) ( ch == 'A' || ch == 'a' || \
                           ch == 'C' || ch == 'c' || \
                           ch == 'G' || ch == 'g' || \
                           ch == 'T' || ch == 't' )

namespace bleu {
  class ReadUtilities {
  
  public:
    static bool check_read(const std::string & read, unsigned long long _ksize) 
    {
      if (read.length() < _ksize) {
        return false;
      }

      for (unsigned int i = 0; i < read.length(); i++)  {
        if (!is_valid_dna(read[i])) {
          return false;
        }
      }

      return true;
    }
  
  };
}

