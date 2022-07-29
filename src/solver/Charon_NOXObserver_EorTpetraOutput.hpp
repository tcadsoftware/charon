
#ifndef CHARON_NOX_OBSERVER_EORTPETRA_OUTPUT_HPP
#define CHARON_NOX_OBSERVER_EORTPETRA_OUTPUT_HPP

#include "NOX_Abstract_PrePostOperator.H"

#include "Teuchos_RCP.hpp"
#include "Teuchos_dyn_cast.hpp"


#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_GlobalIndexer.hpp"

#include "NOX_Epetra_Vector.H"
#include "Epetra_Vector.h"

#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVector.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapOut.h"

#include "MatrixMarket_Tpetra.hpp"

#include <Kokkos_DefaultNode.hpp>

#include <vector>
#include <string>
#include <cstdio>

namespace charon {

  /** This class outputs the linear system solve at each newton step.
    * To check the convergence of the linear system, in matlab use
    *
      A = mmread(sprintf('J_%d.mm', index));
      f = mmread(sprintf('f_%d.mm', index));
      x = mmread(sprintf('x_%d.mm', index));
      px = mmread(sprintf('prev_x_%d.mm', index));
      % px = mmread(sprintf('x_%d.mm', index-1)); % an alternative

      norm(A*(x-px)+f)/norm(f) % should equal linear convergence
    *
    */
  class NOXObserver_EorTpetraOutput : public NOX::Abstract::PrePostOperator {

  public:

    NOXObserver_EorTpetraOutput()
       : prefix_(""), counter_(0)
    { }

    void runPreIterate(const NOX::Solver::Generic& /* solver */)
    {
    }

    void runPostIterate(const NOX::Solver::Generic& solver)
    {
      using Teuchos::dyn_cast;
      using Teuchos::RCP;

      std::stringstream ss;

      // extract Thyra objects
      const NOX::Thyra::Vector& x = Teuchos::dyn_cast<const NOX::Thyra::Vector>(solver.getSolutionGroup().getX());
      const NOX::Thyra::Vector& px = Teuchos::dyn_cast<const NOX::Thyra::Vector>(solver.getPreviousSolutionGroup().getX());
      const NOX::Thyra::Vector& f = dyn_cast<const NOX::Thyra::Vector>(solver.getPreviousSolutionGroup().getF());
      const RCP<const Thyra::LinearOpBase<double > > th_J = Teuchos::dyn_cast<const NOX::Thyra::Group>(solver.getPreviousSolutionGroup()).getJacobianOperator();
      const RCP<const Thyra::VectorBase<double> > th_x = x.getThyraRCPVector();
      const RCP<const Thyra::VectorBase<double> > th_px = px.getThyraRCPVector();
      const RCP<const Thyra::VectorBase<double> > th_f = f.getThyraRCPVector();

      const RCP<const Thyra::BlockedLinearOpBase<double> > blo_J = castOrCreateBlockedLinearOpBase(th_J);
      const RCP<const Thyra::ProductVectorBase<double> > blo_x = Thyra::castOrCreateProductVectorBase(th_x);
      const RCP<const Thyra::ProductVectorBase<double> > blo_px = Thyra::castOrCreateProductVectorBase(th_px);
      const RCP<const Thyra::ProductVectorBase<double> > blo_f = Thyra::castOrCreateProductVectorBase(th_f);

      int blockRows = blo_J->productRange()->numBlocks();
      int blockCols = blo_J->productDomain()->numBlocks();

      // sanity check
      TEUCHOS_ASSERT(blo_x->productSpace()->numBlocks()==blockCols);
      TEUCHOS_ASSERT(blo_f->productSpace()->numBlocks()==blockRows);

      std::vector<bool> x_written(blockCols,false), f_written(blockRows,false);

      bool requiresBlkStr = !(blockCols==1 && blockRows==1);
      for(int r=0;r<blockRows;r++) {
         // build block identifier string
         ss.str(""); ss << r << "_";
         std::string blkRowStr = requiresBlkStr ? ss.str() : "";

         // grab block row vector
         const RCP<const Thyra::VectorBase<double> > blk_f = blo_f->getVectorBlock(r);
         for(int c=0;c<blockCols;c++) {
            // build block identifier string
            ss.str(""); ss << c << "_";
            std::string blkColStr = requiresBlkStr ? ss.str() : "";

            // build block identifier string
            ss.str(""); ss << r << c << "_";
            std::string blkStr = requiresBlkStr ? ss.str() : "";

            // grab matrix block
            const RCP<const Thyra::LinearOpBase<double> > blk_J = blo_J->getBlock(r,c);

            // nothing to do
            if(blk_J==Teuchos::null)
               continue;

            // grab block column vector
            const RCP<const Thyra::VectorBase<double> > blk_x = blo_x->getVectorBlock(c);
            const RCP<const Thyra::VectorBase<double> > blk_px = blo_px->getVectorBlock(c);


            if(blk_J!=Teuchos::null) {
              ss.str(""); ss << prefix_ << "J_" << blkStr << counter_ << ".mm";
              writeMatrixToFile(ss.str(),blk_J);
            }

            // write solution vector
            if(!x_written[c]) {
               ss.str(""); ss << prefix_ << "x_" << blkColStr << counter_ << ".mm";
               writeVectorToFile(ss.str(),blk_x,blk_J,true);

               ss.str(""); ss << prefix_ << "prev_x_" << blkColStr << counter_ << ".mm";
               writeVectorToFile(ss.str(),blk_px,blk_J,true);

               x_written[c] = true;
            }

            // write residual vector
            if(!f_written[r]) {
               ss.str(""); ss << prefix_ << "f_" << blkRowStr << counter_ << ".mm";
               writeVectorToFile(ss.str(),blk_f,blk_J,false);

               f_written[r] = true;
            }
         }
      }

      // update counter
      counter_++;
    }

    void runPreSolve(const NOX::Solver::Generic& /* solver */) { }

    void runPostSolve(const NOX::Solver::Generic& /* solver */) { }

  private:

    /** Cast a Thyra operator to an epetra one, if case fails return null.
      */
    Teuchos::RCP<const Epetra_Operator> getEpetraOperator(const Thyra::LinearOpBase<double> & op) const
    {
      try {
        const Thyra::EpetraLinearOp & thyra_epetra_op = Teuchos::dyn_cast<const Thyra::EpetraLinearOp>(op);
        return thyra_epetra_op.epetra_op();
      }
      catch(std::bad_cast & bc) {
        return Teuchos::null;
      }
    }

    /** Write a vector to a file of a given name. The operator is required to define the Epetra_Map that is
      * used if this system is an epetra vector. Otherwise Tpetra does not need the operator.
      */
    void writeVectorToFile(const std::string & fileName,
                           const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
                           const Teuchos::RCP<const Thyra::LinearOpBase<double> > & J,bool useDomainMap,bool transpose=false) const
    {
      // now use tpetra to get domain map
      typedef Thyra::TpetraVector<double,int,panzer::GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> Thyra_TpVector;
      Teuchos::RCP<const Thyra_TpVector> th_tp_x = Teuchos::rcp_dynamic_cast<const Thyra_TpVector>(x);
      if(th_tp_x!=Teuchos::null) {
        typedef Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> sparse_matrix;

        Tpetra::MatrixMarket::Writer<sparse_matrix>::writeDenseFile(fileName,th_tp_x->getConstTpetraVector());

        return;
      }

      // first try to compute using Epetra
      const Teuchos::RCP<const Epetra_Operator> ep_J = getEpetraOperator(*J);
      if(ep_J!=Teuchos::null) {
        const Epetra_Map & ep_map = useDomainMap ? ep_J->OperatorDomainMap() : ep_J->OperatorRangeMap();
        const Teuchos::RCP<const Epetra_Vector> ep_x = Thyra::get_Epetra_Vector(ep_map, x);

        EpetraExt::VectorToMatrixMarketFile(fileName.c_str(),*ep_x);

        return;
      }

      // if the vector is a locally replicated, we can also print
      Teuchos::RCP<const Thyra::SpmdVectorBase<double> > spmdV
          = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(x);
      if(spmdV!=Teuchos::null && spmdV->spmdSpace()->isLocallyReplicated()) {
        Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<double> > spmdVS = spmdV->spmdSpace();
        Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > comm = spmdVS->getComm();

        if(comm->getRank()!=0) return; // only root processor needs to write

        Teuchos::ArrayRCP<const double> values;
        spmdV->getLocalData(Teuchos::ptrFromRef(values));

        TEUCHOS_ASSERT_EQUALITY(values.size(),spmdVS->dim());

        std::ofstream fout(fileName.c_str());
        fout << "%%MatrixMarket matrix array real general" << std::endl;
        if(transpose) {
          fout << 1 << " " << spmdVS->dim() << std::endl;
          for(int i=0;i<values.size();i++)
            fout << values[i] << " ";
          fout << std::endl;
        }
        else {
          fout << spmdVS->dim() << " " << 1 << std::endl;
          for(int i=0;i<values.size();i++)
            fout << values[i] << std::endl;
        }

        return;
      }

      // if the vector is not a locally replicated, we can still print
      if(spmdV!=Teuchos::null) {
        Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<double> > spmdVS = spmdV->spmdSpace();
        Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > comm = spmdVS->getComm();

        Teuchos::ArrayRCP<const double> values;
        spmdV->getLocalData(Teuchos::ptrFromRef(values));

        int offset = spmdVS->localOffset() + 1;

        FILE* f = 0;

        if(comm->getRank()==0) {
          f = fopen(fileName.c_str(),"w");
          fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
          if(transpose) {
            fprintf(f, "1 %ld %ld\n", spmdVS->dim(), spmdVS->dim());
          } else {
            fprintf(f, "%ld 1 %ld\n", spmdVS->dim(), spmdVS->dim());
          }
          fclose(f);
        }
        comm->barrier();
        f = fopen(fileName.c_str(),"a");

        if(transpose) {
          for(int i=0;i<values.size();i++)
            fprintf(f, "%d %d %22.16e\n", 1, i+offset, values[i]);
        }
        else {
          for(int i=0;i<values.size();i++)
            fprintf(f, "%d %d %22.16e\n", i+offset, 1, values[i]);
        }
        fclose(f);

        return;
      }

      std::cout << "X = " << Teuchos::describe(*x,Teuchos::VERB_EXTREME) << std::endl;
      TEUCHOS_ASSERT(false);
    }

    /** Write a operator to a file of a given name. This automatically handles the Epetra_RowMatrix versus
      * the Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,...> case.
      */
    void writeMatrixToFile(const std::string & fileName,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & J) const
    {
      if(Teuchos::rcp_dynamic_cast<const Thyra::DefaultZeroLinearOp<double> >(J)!=Teuchos::null)
        return;

      // test and conditionally write Tpetra
      typedef Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> th_tp_matrix;
      Teuchos::RCP<const th_tp_matrix> th_tp_J = Teuchos::rcp_dynamic_cast<const th_tp_matrix>(J);

      if(th_tp_J!=Teuchos::null) {
        typedef Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> sparse_matrix;

        Teuchos::RCP<const sparse_matrix> tp_J = Teuchos::rcp_dynamic_cast<const sparse_matrix>(th_tp_J->getConstTpetraOperator());
        Tpetra::MatrixMarket::Writer<sparse_matrix>::writeSparseFile(fileName,tp_J);

        return;
      }

      // test and conditionally write Epetra
      const Teuchos::RCP<const Epetra_Operator> ep_J = getEpetraOperator(*J);
      if(ep_J!=Teuchos::null) {
        Teuchos::RCP<const Epetra_RowMatrix> ep_rowJ = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(ep_J,true);
	std::string fileNameMap="epetraMap.map";
	EpetraExt::BlockMapToMatrixMarketFile (fileNameMap.c_str(),ep_rowJ->Map());
        EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(),*ep_rowJ);

        return;
      }

      // this next bit handles the case when the matrix is actually a locally replicated vector (not multi!)
      // first thing is to determine if the transpose has been taken, if so that will be passed on to the
      // writeVectorToFile method

      double J_scalar = 0.0;
      Thyra::EOpTransp J_transp = Thyra::NOTRANS;
      Teuchos::RCP<const Thyra::LinearOpBase<double> > op_J;
      Thyra::unwrap( J, &J_scalar, &J_transp, &op_J );
      TEUCHOS_ASSERT(J_scalar==1.0);

      // if the vector is a locally replicated, we can also print
      Teuchos::RCP<const Thyra::MultiVectorBase<double> > vec = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(op_J);
      if(vec!=Teuchos::null) {
        writeVectorToFile(fileName,vec->col(0),J,true,J_transp!=Thyra::NOTRANS);
        return;
      }

      std::cout << "J = " << Teuchos::describe(*op_J,Teuchos::VERB_EXTREME) << std::endl;
      TEUCHOS_ASSERT(false); // this means no matrix is recognized
    }

    Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> >
    castOrCreateBlockedLinearOpBase(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
    {
       TEUCHOS_ASSERT(A!=Teuchos::null);

       // try to cast it to blocked linear op
       Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > blo
          = Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(A);

       if(blo!=Teuchos::null)
          return blo;

       // build a 1x1 blocked linear operator and cast it to the base type
       return Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(Thyra::block1x1(A),true);
    }

    std::string prefix_;
    std::size_t counter_;

  };
}

#endif
