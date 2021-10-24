
#ifndef _CHARON_ASCII_FILEREADER_H_
#define _CHARON_ASCII_FILEREADER_H_

#include <vector>

#include "Charon_FileReader.hpp"

namespace charon {

/**
 * @brief A class for reading files containing CSV-like ASCII data files.
 *
 * This class will read in a CSV-style ASCII data file, with an
 * arbitrary number of columns and rows, containing real values. Note
 * that the class does NOT have an extensive error-detection capability
 * so it is important that the file be consistently formatted, e.g.,
 * same number of columns in every row. The data columns can be
 * delimited by a comma, space, tab or semicolon.
 */
class ASCII_FileReader : public FileReader
{

public:

  /**
   * @brief Instantiate and perform the read of the specified ASCII
   * file.
   */
  ASCII_FileReader(std::string const& fileName);

  virtual ~ASCII_FileReader();

  /**
   * @brief Returns the number of columns in the associated ASCII file.
   */
  unsigned int columnCount() const {return num_cols;}

  /**
   * @brief Returns the number of rows of data in the associated ASCII
   * file.
   */
  unsigned int rowCount() const;

  /**
   * @brief Returns the data values as read from the file.
   *
   * The data values are stored by column and "index" is the column of
   * data to return.
   */
  std::vector<double>& dataValues(unsigned int const& index);

protected:

  /**
   * @brief Read the data.
   *
   * Read the data from the file.
   */
  virtual void read(void);

  /**
   * @brief Attempt to determine the number of columns of data within
   * the file.
   */
  void find_column_count();

private:

  /**
   * @brief The empty constructor is private since it's use is
   * discouraged.
   */
  ASCII_FileReader();

  std::string file_name;

  int num_cols;

  /**
   * @brief The data is currently stored in a vector of a vector of
   * doubles.
   */
  std::vector<std::vector<double> > value_holder;

};

}

#endif // _CHARON_ASCII_FILEREADER_H_
