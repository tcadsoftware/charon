
#ifndef CHARON_EMPIRICALDAMAGE_DATA
#define CHARON_EMPIRICALDAMAGE_DATA

#include <map>

#include "Teuchos_RCP.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_MultiVariateParameter_decl.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

namespace charon {

/**
 * @brief Utility class for dealing with the empirical radiation damage model.
 *
 * This class is used for runs involving the empirical radiation damage
 * model. At present it's primary purpose is to track applied potentials
 * for the lookup of mu in the mu table of the empirical model.
 */
class EmpiricalDamage_Data
{
public:

  /**
   * @brief add a contact with a dirichlet condition to the list of
   * tracked contacts.
   *
   * This function adds a contact with a specified voltage to the list
   * of tracked contacts. Note that "contactID" must be
   * unique. Generally this is the sideset ID of the contact.
   */
  void addDirichletContact(std::string contactID,
                           double voltage)
  {
    // If it hasn't been added then do so here
    if (contactList.find(contactID) == contactList.end())
      contactList[contactID] = Teuchos::rcp(new DirichletContact(contactID, voltage));
  }

  /**
   * @brief add a contact with a constant-current condition to the list
   * of tracked contacts.
   *
   * This function adds a contact with a specified current to the list
   * of tracked contacts. Note that "contactID" must be
   * unique. Generally this is the sideset ID of the contact. The
   * voltage on such a contact can vary during the simulation.
   */
  void addConstantCurrentContact(std::string contactID,
                                 double initVoltage,
                                 Teuchos::RCP<panzer::GlobalData> globalData)
  {
    // If it hasn't been added then do so here
    if (contactList.find(contactID) == contactList.end())
      contactList[contactID] = Teuchos::rcp(new ConstantCurrentContact(contactID, initVoltage, globalData));
  }

  /**
   * @brief return the applied bias at the specified contact.
   */
  double getVoltage(std::string contactID)
  {
    return contactList[contactID]->getVoltage();
  }

private:

  /**
   * @brief interface for a contact
   */
  class GenericContact
  {
  public:
    virtual ~GenericContact(){;}
    virtual double getVoltage() = 0;

  };

  /**
   * @brief Class for contacts where a constant voltage is being set.
   */
  class DirichletContact : public GenericContact
  {
  public:

    /**
     * @brief just set the voltage to whatever was specified by the user
     */
    DirichletContact(std::string contact, double voltage)
    {
      contactID = contact;
      voltageValue = voltage;
    }

    ~DirichletContact(){;}

    double getVoltage()
    {
      return voltageValue;
    }

  private:
    std::string contactID;
    double voltageValue;
  };

  /**
   * @brief Class for contacts where a constant-current boundary condition is
   * applied.
   */
  class ConstantCurrentContact : public GenericContact
  {
  public:

    /**
     * @brief the GlobalData class is where the parameters are stored,
     * including applied potential during a constant current run.
     */
    ConstantCurrentContact(std::string contact,
                           double voltage,
                           Teuchos::RCP<panzer::GlobalData> glob)
    {
      contactID = contact;
      initialVoltage = voltage;
      globalData = glob;

      // parameter name in the list. This is hard coded as it is in
      // Charon_BCStrategy_Dirichlet_CurrentConstraint_impl.hpp. If it
      // changes it'll have to be changed here and in there.
      paramName = contactID + "ConstantCurrentVoltage";
    }

    ~ConstantCurrentContact(){;}

    double getVoltage()
    {
      // Don't need a non-residual type for this since it's simply used
      // to look up a value of mu in a table. No Jacobian entry is
      // associated with this calculation so just use the residual
      // type. This also allows this class to avoid being templated.
      Teuchos::RCP<panzer::ScalarParameterEntry<panzer::Traits::Residual> > voltparam =
        panzer::accessScalarParameter<panzer::Traits::Residual>(paramName, *globalData->pl);

      return voltparam->getValue();
    }

  private:
    std::string contactID;
    double initialVoltage;
    Teuchos::RCP<panzer::GlobalData> globalData;

    std::string paramName;
  };

  /**
   * @brief map to store the pointers to the contact calculation classes
   */
  std::map<std::string,Teuchos::RCP<class GenericContact> > contactList;
};

}

#endif // CHARON_EMPIRICALDAMAGE_DATA

