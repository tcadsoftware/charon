
#include <fstream>
#include <string>
#include <vector>

#include <iostream>

#include "Teuchos_TestForException.hpp"

#include "Charon_ASCII_FileReader.hpp"

#define DELIMS " \t,;"

//========================================================================
//========================================================================
charon::ASCII_FileReader::ASCII_FileReader(std::string const& in_fileName) :
  file_name(in_fileName),
  num_cols(0)
{
  read();
}

//========================================================================
//========================================================================
unsigned int charon::ASCII_FileReader::rowCount() const
{
  std::vector<double> const& values = value_holder[0];

  return values.size();
}

//========================================================================
//========================================================================
std::vector<double>&
charon::ASCII_FileReader::dataValues(unsigned int const& index)
{
  std::vector<double>& values = value_holder[index];
  return values;
}

//========================================================================
//========================================================================
void charon::ASCII_FileReader::read()
{

  // Count the number of columns
  find_column_count();

  // Open the file
  std::ifstream file_stream(file_name.c_str());
  if (!file_stream) {
    std::string err_string = "Could not open file named: " + file_name;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, err_string);
  }

  // Initialize the storage for the values. Each column is stored in
  // its own vector
  for (int col=0; col < num_cols; ++col) {
    std::vector<double> values;
    value_holder.push_back(values);
  }

  // Read the complete file
  while (!file_stream.eof()) {

    for (int col=0; col < num_cols; ++col) {

      std::vector<double>& values = value_holder[col];

      double tmp_value;
      file_stream >> tmp_value;

      // This check avoids problems with blank lines and any other
      // garbage that can't be pushed into "tmp_value".
      if (file_stream)
        values.push_back(tmp_value);
    }
  }

}

//========================================================================
//========================================================================
void charon::ASCII_FileReader::find_column_count()
{

  // Open the file
  std::ifstream file_stream(file_name.c_str());
  if (!file_stream) {
    std::string err_string = "Could not open file named: " + file_name;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, err_string);
  }

  // Read an entire line
  std::string line;
  std::getline(file_stream, line);

  // Count entries/columns in the first row of the file
  std::string const delims(DELIMS);

  std::string::size_type beg_idx = line.find_first_not_of(delims);

  while (beg_idx != std::string::npos) {

    std::string::size_type end_idx = line.find_first_of(delims, beg_idx);

    if (end_idx == std::string::npos) {
      end_idx = line.length();
    }

    beg_idx = line.find_first_not_of(delims, end_idx);
    num_cols++;
  }

}

//========================================================================
//========================================================================
charon::ASCII_FileReader::ASCII_FileReader()
{
}

//========================================================================
//========================================================================
charon::ASCII_FileReader::~ASCII_FileReader()
{
}
