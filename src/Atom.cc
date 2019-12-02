#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <omp.h>
#include "IMSRG.hh"
#include "PhysicalConstants.hh"
#include "Parameters.hh"
#include "Atom.hh"

namespace Atom
{
  int MainAtom(int argc, char** argv)
  {
    Parameters parameters(argc,argv);
    if (parameters.help_mode) return 0;

    std::string inputtbme = parameters.s("2bme");
    std::string reference = parameters.s("reference");
    std::string valence_space = parameters.s("valence_space");
    std::string custom_valence_space = parameters.s("custom_valence_space");
    std::string basis = parameters.s("basis");
    std::string basis_type = parameters.s("basis_type");
    std::string method = parameters.s("method");
    std::string flowfile = parameters.s("flowfile");
    std::string intfile = parameters.s("intfile");
    std::string core_generator = parameters.s("core_generator");
    std::string valence_generator = parameters.s("valence_generator");
    std::string fmt2 = parameters.s("fmt2");
    std::string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
    std::string scratch = parameters.s("scratch");
    std::string use_brueckner_bch = parameters.s("use_brueckner_bch");
    std::string valence_file_format = parameters.s("valence_file_format");
    std::string occ_file = parameters.s("occ_file");
    std::string goose_tank = parameters.s("goose_tank");
    std::string write_omega = parameters.s("write_omega");
    std::string IMSRG3 = parameters.s("IMSRG3");
    std::string freeze_occupations = parameters.s("freeze_occupations");
    std::string hunter_gatherer = parameters.s("hunter_gatherer");
    std::string relativistic_correction = parameters.s("relativistic_correction");
    bool use_NAT_occupations = (parameters.s("use_NAT_occupations")=="true") ? true : false;
    int eMax = parameters.i("emax");
    int lmax = parameters.i("lmax"); // so far I only use this with atomic systems.
    int file2e1max = parameters.i("file2e1max");
    int file2e2max = parameters.i("file2e2max");
    int file2lmax = parameters.i("file2lmax");
    int targetMass = parameters.i("A");
    int nsteps = parameters.i("nsteps");

    double hw = parameters.d("hw");
    double smax = parameters.d("smax");
    double ode_tolerance = parameters.d("ode_tolerance");
    double dsmax = parameters.d("dsmax");
    double ds_0 = parameters.d("ds_0");
    double domega = parameters.d("domega");
    double omega_norm_max = parameters.d("omega_norm_max");
    double denominator_delta = parameters.d("denominator_delta");
    double eta_criterion = parameters.d("eta_criterion");

    std::vector<std::string> opnames = parameters.v("Operators");
    std::vector<std::string> opsfromfile = parameters.v("OperatorsFromFile");
    std::vector<std::string> opnamesPT1 = parameters.v("OperatorsPT1");
    std::vector<std::string> opnamesRPA = parameters.v("OperatorsRPA");
    std::vector<std::string> opnamesTDA = parameters.v("OperatorsTDA");

    std::vector<Operator> ops;
    std::vector<std::string> spwf = parameters.v("SPWF");

    std::ifstream test;
    // test 2bme file
    if (inputtbme != "none" and fmt2.find("oakridge")==std::string::npos and fmt2 != "schematic" )
    {
      test.open(inputtbme);
      ////    if( not test.good() and fmt2!="oakridge_binary")
      //    if( not test.good() and  fmt2.find("oakridge")== std::string::npos)
      if( not test.good() )
      {
        std::cout << "trouble reading " << inputtbme << "  fmt2 = " << fmt2 << "   exiting. " << std::endl;
        return 1;
      }
      test.close();
    }

    ReadWrite rw;
    rw.SetScratchDir(scratch);
    // Test whether the scratch directory exists and we can write to it.
    // This is necessary because otherwise you get garbage for transformed operators and it's
    // not obvious what went wrong.
    if ( method=="magnus" and  scratch != "" and scratch!= "/dev/null" and scratch != "/dev/null/")
    {
      std::string testfilename = scratch + "/_this_is_a_test_delete_me";
      std::ofstream testout(testfilename);
      testout << "PASSED" << std::endl;
      testout.close();
      std::remove( testfilename.c_str() );
      if ( not testout.good() )
      {
        std::cout << "WARNING in " << __FILE__ <<  " failed test write to scratch directory " << scratch;
        if (opnames.size()>0 )
        {
          std::cout << "   dying now. " << std::endl;
          exit(EXIT_FAILURE);
        }
        else std::cout << std::endl;
      }
    }
    if ( (method=="magnus") and (scratch=="/dev/null" or scratch=="/dev/null/") )
    {
      if ( opnames.size() > 0 )
      {
        std::cout << "WARNING!!! using Magnus with scratch = " << scratch << " but you're also trying to transform some operators: ";
        for (auto opn : opnames ) std::cout << opn << " ";
        std::cout << "   dying now." << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    //  ModelSpace modelspace;

    if (custom_valence_space!="") // if a custom space is defined, the input valence_space is just used as a name
    {
      if (valence_space=="") // if no name is given, then just name it "custom"
      {
        parameters.string_par["valence_space"] = "custom";
        flowfile = parameters.DefaultFlowFile();
        intfile = parameters.DefaultIntFile();
      }
      valence_space = custom_valence_space;
    }

    ModelSpace modelspace;
    modelspace.SetLmax(lmax);
    modelspace.InitAtomicSpace(eMax, basis_type, reference, valence_space);
    if (nsteps < 0) nsteps = modelspace.valence.size()>0 ? 2 : 1;
    Operator Hbare = Operator(modelspace,0,0,0,2);
    Hbare.SetHermitian();
    if ( goose_tank == "true" or goose_tank == "True") Commutator::SetUseGooseTank(true);
    //
    std::cout << "Reading interactions..." << std::endl;
    if (inputtbme != "none")
    {
      if (fmt2 == "tokyo") rw.ReadTokyoAtomic(inputtbme,Hbare);
    }

    if (inputtbme == "none")
    {
      using PhysConst::M_ELECTRON;
      using PhysConst::M_NUCLEON;
      int Z = modelspace.GetTargetZ() ;
      Hbare -= Z*imsrg_util::VCentralCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;
      Hbare += imsrg_util::VCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;  // convert oscillator length from fm with nucleon mass to nm with electon mass (in eV).
      Hbare += imsrg_util::KineticEnergy_Op(modelspace); // Don't need to rescale this, because it's related to the oscillator frequency, which we input.
      Hbare /= PhysConst::HARTREE; // Convert to Hartree
    }

    std::cout << "Creating HF" << std::endl;
    HFMBPT hf(Hbare); // HFMBPT inherits from HartreeFock, so no harm done here.
    if (freeze_occupations == "false" )  hf.UnFreezeOccupations();
    std::cout << "Solving" << std::endl;
    hf.Solve();

    //  Operator HNO;
    Operator& HNO = Hbare;
    if (basis == "HF" and method !="HF")
    {
      HNO = hf.GetNormalOrderedH();
    }
    else if (basis == "NAT") // we want to use the natural orbital basis
    {
      hf.UseNATOccupations( use_NAT_occupations );
      hf.GetNaturalOrbitals();
      HNO = hf.GetNormalOrderedHNAT();

      if (use_NAT_occupations)
      {
        hf.FillLowestOrbits();
        std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
        HNO = HNO.UndoNormalOrdering();
        hf.UpdateReference();
        modelspace.SetReference(modelspace.core); // change the reference
        std::cout << "Doing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
        HNO = HNO.DoNormalOrdering();
      }
    }
    else if (basis == "oscillator")
    {
      HNO = Hbare.DoNormalOrdering();
    }

    if (method != "HF")
    {
      std::cout << "Perturbative estimates of gs energy:" << std::endl;
      double EMP2 = HNO.GetMP2_Energy();
      std::cout << "EMP2 = " << EMP2 << std::endl;
      std::array<double,3> Emp_3 = HNO.GetMP3_Energy();
      double EMP3 = Emp_3[0]+Emp_3[1]+Emp_3[2];
      std::cout << "E3_pp = " << Emp_3[0] << "  E3_hh = " << Emp_3[1] << " E3_ph = " << Emp_3[2] << "   EMP3 = " << EMP3 << std::endl;
      std::cout << "To 3rd order, E = " << HNO.ZeroBody+EMP2+EMP3 << std::endl;
    }

    if (basis=="HF" or basis=="NAT")
    {
      std::cout << basis << " Single particle energies and wave functions:" << std::endl;
      hf.PrintSPEandWF();
      std::cout << std::endl;
    }

    if ( method == "HF" or method == "MP3")
    {
      HNO.PrintTimes();
      return 0;
    }

    IMSRGSolver imsrgsolver(HNO);
    imsrgsolver.SetReadWrite(rw);
    imsrgsolver.SetEtaCriterion(eta_criterion);
    imsrgsolver.max_omega_written = 500;
    bool brueckner_restart = false;
    if (hunter_gatherer=="true") imsrgsolver.SetHunterGatherer( true);

    if (method == "NSmagnus") // "No split" magnus
    {
      omega_norm_max=50000;
      method = "magnus";
    }
    if (method.find("brueckner") != std::string::npos)
    {
      if (method=="brueckner2") brueckner_restart=true;
      if (method=="brueckner1step")
      {
        nsteps = 1;
        core_generator = valence_generator;
      }
      use_brueckner_bch = "true";
      omega_norm_max=500;
      method = "magnus";
    }

    if (use_brueckner_bch == "true" or use_brueckner_bch == "True")
    {
      //    Hbare.SetUseBruecknerBCH(true);
      //    HNO.SetUseBruecknerBCH(true);
      Commutator::SetUseBruecknerBCH(true);
      std::cout << "Using Brueckner flavor of BCH" << std::endl;
    }
    if (IMSRG3 == "true")
    {
      Commutator::SetUseIMSRG3(true);
      std::cout << "Using IMSRG(3) commutators. This will probably be slow..." << std::endl;
    }
    imsrgsolver.SetMethod(method);
    //  imsrgsolver.SetHin(Hbare);
    imsrgsolver.SetHin(HNO);
    imsrgsolver.SetSmax(smax);
    imsrgsolver.SetFlowFile(flowfile);
    imsrgsolver.SetDs(ds_0);
    imsrgsolver.SetDsmax(dsmax);
    imsrgsolver.SetDenominatorDelta(denominator_delta);
    imsrgsolver.SetdOmega(domega);
    imsrgsolver.SetOmegaNormMax(omega_norm_max);
    imsrgsolver.SetODETolerance(ode_tolerance);
    if (denominator_delta_orbit != "none")
      imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

    imsrgsolver.SetGenerator(core_generator);
    if (core_generator.find("imaginary")!=std::string::npos or core_generator.find("wegner")!=std::string::npos )
    {
      if (ds_0>1e-2)
      {
        ds_0 = 1e-4;
        dsmax = 1e-2;
        imsrgsolver.SetDs(ds_0);
        imsrgsolver.SetDsmax(dsmax);
      }
    }

    imsrgsolver.Solve();

    if (IMSRG3 == "true")
    {
      std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
    }

    if (brueckner_restart)
    {
      arma::mat newC = hf.C * arma::expmat( -imsrgsolver.GetOmega(0).OneBody  );
      //     if (input3bme != "none") Hbare.SetParticleRank(3);
      HNO = hf.GetNormalOrderedH(newC);
      imsrgsolver.SetHin(HNO);
      imsrgsolver.s = 0;
      imsrgsolver.Solve();
    }

    if (nsteps > 1 and valence_space != reference) // two-step decoupling, do core first
    {
      if (method == "magnus") smax *= 2;

      imsrgsolver.SetGenerator(valence_generator);
      std::cout << "Setting generator to " << valence_generator << std::endl;
      modelspace.ResetFirstPass();
      if (valence_generator.find("imaginary")!=std::string::npos or valence_generator.find("wegner")!=std::string::npos)
      {
        if (ds_0>1e-2)
        {
          ds_0 = 1e-4;
          dsmax = 1e-2;
          imsrgsolver.SetDs(ds_0);
          imsrgsolver.SetDsmax(dsmax);
        }
      }
      imsrgsolver.SetSmax(smax);
      imsrgsolver.Solve();
    }

    // Transform all the operators
    if (method == "magnus")
    {
      if (ops.size()>0) std::cout << "transforming operators" << std::endl;
      for (size_t i=0;i<ops.size();++i)
      {
        std::cout << opnames[i] << " " << std::endl;
        ops[i] = imsrgsolver.Transform(ops[i]);
        std::cout << " (" << ops[i].ZeroBody << " ) " << std::endl;
        //      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");
      }
      std::cout << std::endl;
      // increase smax in case we need to do additional steps
      smax *= 1.5;
      imsrgsolver.SetSmax(smax);
    }
    if (method == "flow" or method == "flow_RK4" )
    {
      for (size_t i=0;i<ops.size();++i)
      {
        ops[i] = imsrgsolver.GetOperator(i+1);  // the zero-th operator is the Hamiltonian
      }
    }

    // If we're doing targeted/ensemble normal ordering
    // we now re-normal order wrt to the core
    // and do any remaining flow.
    ModelSpace ms2(modelspace);
    bool renormal_order = false;
    if (modelspace.valence.size() > 0 )
      //  if (modelspace.valence.size() > 0 or basis=="NAT")
    {
      renormal_order = modelspace.holes.size() != modelspace.core.size();
      if (not renormal_order)
      {
        for (auto c : modelspace.core)
        {
          if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (std::abs(1-modelspace.GetOrbit(c).occ)>1e-6))
          {
            renormal_order = true;
            break;
          }
        }
      }
    }
    if ( renormal_order )
    {

      HNO = imsrgsolver.GetH_s();

      int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
      std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
      HNO = HNO.UndoNormalOrdering();

      ms2.SetReference(ms2.core); // change the reference
      HNO.SetModelSpace(ms2);

      std::cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << std::endl;
      HNO = HNO.DoNormalOrdering();

      // More flowing is unnecessary, since things should stay decoupled.
      imsrgsolver.SetHin(HNO);
      //    imsrgsolver.SetEtaCriterion(1e-4);
      //    imsrgsolver.Solve();
      // Change operators to the new basis, then apply the rest of the transformation
      std::cout << "Final transformation on the operators..." << std::endl;
      int iop = 0;
      for (auto& op : ops)
      {
        std::cout << opnames[iop++] << std::endl;
        op = op.UndoNormalOrdering();
        op.SetModelSpace(ms2);
        op = op.DoNormalOrdering();
        // transform using the remaining omegas
        op = imsrgsolver.Transform_Partial(op,nOmega);
      }
    }

    // Write the output

    // If we're doing a shell model interaction, write the
    // interaction files to disk.
    if (modelspace.valence.size() > 0)
    {
      if (valence_file_format == "antoine") // this is still being tested...
      {
        rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
        rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
      }
      std::cout << "Writing files: " << intfile << std::endl;
      //rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
      //rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
      rw.WriteTokyo(imsrgsolver.GetH_s(),intfile+".snt", "");

      if (method == "magnus" or method=="flow_RK4")
      {
        for (index_t i=0;i<ops.size();++i)
        {
          if ( ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1) and (ops[i].GetNumberLegs()%2==0) )
          {
            //rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
            rw.WriteTokyo(ops[i],intfile+opnames[i]+".snt", "op");
          }
          else if ( ops[i].GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
          {
            //            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
            rw.WriteDaggerOperator( ops[i], intfile+opnames[i]+".dag",opnames[i]);
          }
          else
          {
            //rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
            //rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
            rw.WriteTensorTokyo(intfile+opnames[i]+"_2b.snt",ops[i]);
          }
        }
      }
    }
    else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
    {
      std::cout << "Core Energy = " << std::setprecision(6) << imsrgsolver.GetH_s().ZeroBody << std::endl;
      for (index_t i=0;i<ops.size();++i)
      {
        Operator& op = ops[i];
        std::cout << opnames[i] << " = " << ops[i].ZeroBody << std::endl;
        if ( opnames[i] == "Rp2" )
        {
          int Z = modelspace.GetTargetZ();
          int A = modelspace.GetTargetMass();
          std::cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
          std::cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << std::endl;
        }
        if ((op.GetJRank()>0) or (op.GetTRank()>0)) // if it's a tensor, you probably want the full operator
        {
          std::cout << "Writing operator to " << intfile+opnames[i]+".op" << std::endl;
          rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
        }
      }
    }


    //  std::cout << "Made it here and write_omega is " << write_omega << std::endl;
    if (write_omega == "true" or write_omega == "True")
    {
      std::cout << "writing Omega to " << intfile << "_omega.op" << std::endl;
      rw.WriteOperatorHuman(imsrgsolver.Omega.back(),intfile+"_omega.op");
    }


    if (IMSRG3 == "true")
    {
      std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
    }
    Hbare.PrintTimes();

    return 0;

  }
}