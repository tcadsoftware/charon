
#ifndef CHARON_MATERIAL_PROPERTIES_HPP
#define CHARON_MATERIAL_PROPERTIES_HPP

#include <string>
#include <vector>
#include "Teuchos_ParameterList.hpp"


namespace charon
{
  
  class Material_Properties;

  // Compound Material
  class CompoundMaterial {
  public:
    CompoundMaterial(const std::string& mat_name, 
                     Material_Properties& matProperty,
		     const std::string& mat_arity, const std::string& mat_cat, 
		     const std::string& mat1, const std::string& mat2, 
		     const std::string& mat3="");

    template<typename EvalT> typename EvalT::ScalarT
    compute_Eg(double xFrac, double yFrac, typename EvalT::ScalarT T);
   
    template<typename EvalT> typename EvalT::ScalarT
    compute_Chi(double xFrac, double yFrac, typename EvalT::ScalarT T);

    template<typename EvalT> typename EvalT::ScalarT
    compute_Eps(double xFrac, double yFrac, typename EvalT::ScalarT T);

    template<typename EvalT> typename EvalT::ScalarT
    compute_eDOS(double xFrac, double yFrac, typename EvalT::ScalarT T);

    template<typename EvalT> typename EvalT::ScalarT
    compute_hDOS(double xFrac, double yFrac, typename EvalT::ScalarT T);

    void setEg_Eg300(double b, double c, double Eg300_0, double Eg300_1) { 
      Eg_Eg300_b = b; Eg_Eg300_c = c; Eg_Eg300_0 = Eg300_0; Eg_Eg300_1 = Eg300_1;
    } 
    void setEg_Chi300(double b, double c, double Chi300_0, double Chi300_1) { 
      Eg_Chi300_b = b; Eg_Chi300_c = c; Eg_Chi300_0 = Chi300_0; Eg_Chi300_1 = Chi300_1;
    } 
    void setEg_alpha(double b, double c, double alpha_0, double alpha_1) { 
      Eg_alpha_b = b; Eg_alpha_c = c; Eg_alpha_0 = alpha_0; Eg_alpha_1 = alpha_1;
    } 
    void setEg_beta(double b, double c, double beta_0, double beta_1) { 
      Eg_beta_b = b; Eg_beta_c = c; Eg_beta_0 = beta_0; Eg_beta_1 = beta_1;
    } 
    void setRelPerm_Value(double b, double c, double Value_0, double Value_1) { 
      RelPerm_Value_b = b; RelPerm_Value_c = c; 
      RelPerm_Value_0 = Value_0; RelPerm_Value_1 = Value_1;
    } 
    void setDOS_Nc300(double Nc300_b, double Nc300_c, double Nc300_0, double Nc300_1) {
      DOS_Nc300_b = Nc300_b; DOS_Nc300_c = Nc300_c;
      DOS_Nc300_0 = Nc300_0; DOS_Nc300_1 = Nc300_1;
    }
    void setDOS_Nv300(double Nv300_b, double Nv300_c, double Nv300_0, double Nv300_1) {
      DOS_Nv300_b = Nv300_b; DOS_Nv300_c = Nv300_c;
      DOS_Nv300_0 = Nv300_0; DOS_Nv300_1 = Nv300_1;
    }
    void setDOS_Nc_F(double Nc_F_b, double Nc_F_c, double Nc_F_0, double Nc_F_1) {
      DOS_Nc_F_b = Nc_F_b; DOS_Nc_F_c = Nc_F_c;
      DOS_Nc_F_0 = Nc_F_0; DOS_Nc_F_1 = Nc_F_1;
    }
    void setDOS_Nv_F(double Nv_F_b, double Nv_F_c, double Nv_F_0, double Nv_F_1) {
      DOS_Nv_F_b = Nv_F_b; DOS_Nv_F_c = Nv_F_c;
      DOS_Nv_F_0 = Nv_F_0; DOS_Nv_F_1 = Nv_F_1;
    }


  private:
    std::string name_id;
    Material_Properties& matProperties;
    std::string arity;
    std::string category;
    std::string sideMat1;
    std::string sideMat2;
    std::string sideMat3;

    // bandgap side materials user overwrite
    double Eg_Eg300_0;
    double Eg_Eg300_1;
    double Eg_Chi300_0;
    double Eg_Chi300_1;
    double Eg_alpha_0;
    double Eg_alpha_1;
    double Eg_beta_0;
    double Eg_beta_1;
    // relative permittivity user overwrite
    double RelPerm_Value_0;
    double RelPerm_Value_1;
    // effective DOS overwrite
    double DOS_Nc300_0;
    double DOS_Nc300_1;
    double DOS_Nv300_0;
    double DOS_Nv300_1;
    double DOS_Nc_F_0;
    double DOS_Nc_F_1;
    double DOS_Nv_F_0;
    double DOS_Nv_F_1;

    // bandgap interpolation
    double Eg_Eg300_b;
    double Eg_Eg300_c;
    double Eg_Chi300_b;
    double Eg_Chi300_c;
    double Eg_alpha_b;
    double Eg_alpha_c;
    double Eg_beta_b;
    double Eg_beta_c;
    // relative permittivity interpolation
    double RelPerm_Value_b;
    double RelPerm_Value_c;
    // effective DOS interpolation
    double DOS_Nc300_b;
    double DOS_Nc300_c;
    double DOS_Nv300_b;
    double DOS_Nv300_c;
    double DOS_Nc_F_b;
    double DOS_Nc_F_c;
    double DOS_Nv_F_b;
    double DOS_Nv_F_c;
  }; // Compound Material


  // Binary Compound Material
  class BinaryCompoundMaterial : public CompoundMaterial {
  public:
    BinaryCompoundMaterial(const std::string& name, 
			   Material_Properties& matProperty,
			   const std::string& cat, const std::string& mat1, 
			   const std::string& mat2);
  private:
  }; 

  // Ternary Compound Material
  class TernaryCompoundMaterial : public CompoundMaterial {
  public:  
    TernaryCompoundMaterial(const std::string& name, 
                            Material_Properties& matProperty,
			    const std::string& cat, const std::string& mat1, 
			    const std::string& mat2);
  private:
  }; 

  // Quaternary Compound Material
  class QuaternaryCompoundMaterial : public CompoundMaterial {
  public:  
    QuaternaryCompoundMaterial(const std::string& name, 
                               Material_Properties& matProperty,
			       const std::string& cat, const std::string& mat1, 
			       const std::string& mat2, const std::string& mat3);
  private:
  }; 



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

    // check if a material is mole-fraction dependent
    bool hasMoleFracDependence(const std::string& materialName);

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

    // get the arity of a given material
    std::string getArityType(const std::string& materialName);

/*    // check if a material-dependent property value is changed by user
    bool isPropertyValueChanged(const std::string& materialName,
                                const std::string& propertyName);

    // check if a material-independent property is changed
    bool isPropertyValueChanged(const std::string& propertyName);
*/

    // add on-demand mole fraction material properties
    void registerMoleFracMaterial(const std::string& materialName);

    // get a mole fraction material
    const Teuchos::RCP<CompoundMaterial> getMoleFracMaterial(const std::string& materialName);

    void setupMoleFracBandgapParams(const std::string& materialName,
				    double Eg300_b, double Eg300_c,
				    double Chi300_b, double Chi300_c,
				    double alpha_b, double alpha_c,
				    double beta_b, double beta_c,
				    double Eg300_0, double Eg300_1,
                                    double Chi300_0, double Chi300_1,
				    double alpha_0, double alpha_1,
				    double beta_0, double beta_1);

    void setupMoleFracRelPermittivityParams(const std::string& materialName,
					    double Value_b, double Value_c,
					    double Value_0, double Value_1);

    void setupMoleFracDOSParams(const std::string&materialName,
				double Nc300_b, double Nc300_c,
				double Nv300_b, double Nv300_c,
				double Nc_F_b, double Nc_F_c,
				double Nv_F_b, double Nv_F_c,
				double Nc300_0, double Nc300_1,
				double Nv300_0, double Nv300_1,
				double Nc_F_0, double Nc_F_1,
				double Nv_F_0, double Nv_F_1);

    
  private:
    // private ctor invoked internally via getInstance()
    Material_Properties();

    // Set the properties for the various materials.
    void setSiliconParameters(Teuchos::ParameterList& p);
    void setGeParameters(Teuchos::ParameterList& p);
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
    void setAlxGa1minxNParameters(Teuchos::ParameterList& p);
    void setSi1minxGexParameters(Teuchos::ParameterList& p);
    void set4HSiCParameters(Teuchos::ParameterList& p); 

    Teuchos::ParameterList pMaterials;

    // Teuchos::ParameterList pIsPropertyChanged;

    static bool mFirst;

    std::map<std::string,Teuchos::RCP<CompoundMaterial>> moleFracMaterials;

  }; // End "class Material_Properties"

} // END "namespace charon"


#endif
