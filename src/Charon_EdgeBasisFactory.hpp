
#ifndef CHARON_EDGE_BASIS_FACTORY_HPP
#define CHARON_EDGE_BASIS_FACTORY_HPP

#include <sstream>
#include <string>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Intrepid2_Basis.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"

#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HCURL_TET_I1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

namespace charon {

  template <typename ScalarT, typename ArrayT>
  Teuchos::RCP<Intrepid2::Basis<ScalarT,ArrayT> >
  createEdgeBasis(const std::string& cell_name)
  {
    const std::string& name = cell_name;

    Teuchos::RCP<Intrepid2::Basis<ScalarT,ArrayT> > edge_basis;

    if (name == "Triangle_3")
      edge_basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TRI_I1_FEM<ScalarT,ArrayT> );
    else if (name == "Quadrilateral_4")
      edge_basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_QUAD_I1_FEM<ScalarT,ArrayT> );
    else if (name == "Tetrahedron_4")
      edge_basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TET_I1_FEM<ScalarT,ArrayT> );
    else if (name == "Hexahedron_8")
      edge_basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_HEX_I1_FEM<ScalarT,ArrayT> );

    if (Teuchos::is_null(edge_basis))
    {
      std::ostringstream s;
      s << "Failed to create edge basis for cell topology of " << name << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(edge_basis), std::runtime_error, s.str());
    }

    return edge_basis;
  }

}


#endif
