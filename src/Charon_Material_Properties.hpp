
#ifndef CHARON_MATERIAL_PROPERTIES_HPP
#define CHARON_MATERIAL_PROPERTIES_HPP

#include <string>
#include <vector>
#include "Teuchos_ParameterList.hpp"

namespace charon
{

/**
 * @brief A class for storing the names of provided variables for
 * material properties.
 *
 * This class stores the default values of the semiconductor material
 * properties.  The default values may be changed by the user in the input .xml.
 *
 */

  class Material_Properties
  {
  public:

    // return the singleton instance of this class
    static Material_Properties& getInstance();

    // validate the material name, if not found, throw exception
    void validateMaterialName(const std::string& materialName);

    // set the property value for a given material and property name
    void setPropertyValue(const std::string& materialName,
                          const std::string& propertyName,
                          double propertyValue);

    // set the property value for a material-independent property name
    void setPropertyValue(const std::string& propertyName,
                          double propertyValue);

    // get the property value for a given material and property name
    double getPropertyValue(const std::string& materialName,
                            const std::string& propertyName);

    // get the property value for a material-independent property name
    double getPropertyValue(const std::string& propertyName);

    // get the type of a given material
    std::string getMaterialType(const std::string& materialName);

/*    // check if a material-dependent property value is changed by user
    bool isPropertyValueChanged(const std::string& materialName,
                                const std::string& propertyName);

    // check if a material-independent property is changed
    bool isPropertyValueChanged(const std::string& propertyName);
*/

  private:
    // private ctor invoked internally via getInstance()
    Material_Properties();

    // Set the properties for the various materials.
    void setSiliconParameters(Teuchos::ParameterList& p);
    void setSiO2Parameters(Teuchos::ParameterList& p);
    void setOxyNitrideParameters(Teuchos::ParameterList& p);
    void setGaAsParameters(Teuchos::ParameterList& p);
    void setInPParameters(Teuchos::ParameterList& p);
    void setAlGaAsParameters(Teuchos::ParameterList& p);
    void setInGaAsParameters(Teuchos::ParameterList& p);
    void setAlInAsParameters(Teuchos::ParameterList& p);
    void setGaAsPParameters(Teuchos::ParameterList& p);
    void setInGaPParameters(Teuchos::ParameterList& p);
    void setAlGaNParameters(Teuchos::ParameterList& p);
    void setAlNParameters(Teuchos::ParameterList& p);
    void setGaNParameters(Teuchos::ParameterList& p);
    void setTiO2Parameters(Teuchos::ParameterList& p);
    void setTantalumParameters(Teuchos::ParameterList& p);
    void setPlatinumParameters(Teuchos::ParameterList& p);
    void setPlatinumSemiParameters(Teuchos::ParameterList& p);
    void setTa2O5Parameters(Teuchos::ParameterList& p);
    void setTaOParameters(Teuchos::ParameterList& p);

    Teuchos::ParameterList pMaterials;

    // Teuchos::ParameterList pIsPropertyChanged;

    static bool mFirst;

  }; // End "class Material_Properties"

} // END "namespace charon"


#endif
