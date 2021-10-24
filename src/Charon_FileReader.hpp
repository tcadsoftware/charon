
#ifndef _CHARON_FILEREADER_H_
#define _CHARON_FILEREADER_H_

namespace charon {

/**
 * @brief Abstract base class defining the interface for file readers.
 *
 * This class defines the interface for all generic file readers within
 * Charon. It's currently not utilized much, but may be useful if a
 * factory is later utilized to instantiate file reader objects.
 */
class FileReader {

public:

  virtual ~FileReader() {}

  /**
   * @brief The read member function for file readers.
   */
  virtual void read(void) = 0;
};

}

#endif // _CHARON_FILEREADER_H_
