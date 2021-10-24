
///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>

// Teuchos
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"
#include "Teuchos_YamlParser_decl.hpp"

/**
 *  \brief Convert an XML input file to a YAML one.
 *
 *  This utility will take an XML input file (that is, a Teuchos::ParameterList
 *  that has been written out to XML) and convert it to an equivalent YAML
 *  input file.
 *
 *  You call this utility with a single argument, the name of the XML input
 *  file you wish to convert.  For instance, `./xmlToYaml myInputFile.xml`.
 *
 *  \throws `std::invalid_argument` If the wrong number of command line
 *                                  arguments are provided.
 *
 *  \returns 0 if successful.
 */
int
main(
  int   argc,
  char* argv[])
{
  using std::invalid_argument;
  using std::string;
  using Teuchos::getParametersFromXmlFile;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::writeParameterListToYamlFile;
  TEUCHOS_TEST_FOR_EXCEPTION(argc != 2, invalid_argument, "You must call "    \
    "this utility with a single argument, for instance, `./xmlToYaml "        \
    "myInputFile.xml`.")
  string inputFile(argv[1]), prefix(inputFile),
    lastFour(inputFile.substr(inputFile.length() - 4));
  if ((lastFour == ".xml") or (lastFour == ".XML"))
    prefix = inputFile.substr(0, inputFile.find_last_of("."));
  string outputFile(prefix + ".yaml");
  RCP<ParameterList> input = getParametersFromXmlFile(inputFile);
  writeParameterListToYamlFile(*input, outputFile);
  return 0;
} // end of main()
