/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                  ____                                         ///
///        _________________           _____________/   /\               _________________        ///
///       /____/_____/_____/|         /____/_____/ /___/  \             /____/_____/_____/|       ///
///      /____/_____/__G_ /||        /____/_____/|/   /\  /\           /____/_____/____ /||       ///
///     /____/_____/__+__/|||       /____/_____/|/ G /  \/  \         /____/_____/_____/|||       ///
///    |     |     |     ||||      |     |     |/___/   /\  /\       |     |     |     ||||       ///
///    |  I  |  M  |     ||/|      |  I  |  M  /   /\  /  \/  \      |  I  |  M  |     ||/|       ///
///    |_____|_____|_____|/||      |_____|____/ + /  \/   /\  /      |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |   /___/   /\  /  \/       |     |     |     ||||       ///
///    |  S  |  R  |     ||/|      |  S  |   \   \  /  \/   /        |  S  |  R  |  G  ||/|       ///
///    |_____|_____|_____|/||      |_____|____\ __\/   /\  /         |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |     \   \  /  \/          |     |     |     ||||       ///
///    |     |  +  |     ||/       |     |  +  |\ __\/   /           |     |  +  |  +  ||/        ///
///    |_____|_____|_____|/        |_____|_____|/\   \  /            |_____|_____|_____|/         ///
///                                               \___\/                                          ///
///                                                                                               ///
///           imsrg++ : Interface for performing standard IMSRG calculations.                     ///
///                     Usage is imsrg++  option1=value1 option2=value2 ...                       ///
///                     To get a list of options, type imsrg++ help                               ///
///                                                                                               ///
///                                                      - Ragnar Stroberg 2016                   ///
///                                                                                               ///
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
//    imsrg++.cc, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"
#include "PhysicalConstants.hh"
#include "version.hh"

struct OpFromFile {
   std::string file2name,file3name,opname;
   int j,p,t,r; // J rank, parity, dTz, particle rank
};

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
  std::cout << "######  imsrg++ build version: " << version::BuildVersion() << std::endl;

  Parameters parameters(argc,argv);
  if (parameters.help_mode) return 0;

  std::string inputtbme = parameters.s("2bme");
  std::string input3bme = parameters.s("3bme");
  std::string input3bme_type = parameters.s("3bme_type");
  std::string no2b_precision = parameters.s("no2b_precision");
  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");
  std::string custom_valence_space = parameters.s("custom_valence_space");
  std::string basis = parameters.s("basis");
  std::string method = parameters.s("method");
  std::string flowfile = parameters.s("flowfile");
  std::string intfile = parameters.s("intfile");
  std::string core_generator = parameters.s("core_generator");
  std::string valence_generator = parameters.s("valence_generator");
  std::string fmt2 = parameters.s("fmt2");
  std::string fmt3 = parameters.s("fmt3");
  std::string input_op_fmt = parameters.s("input_op_fmt");
  std::string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
  std::string LECs = parameters.s("LECs");
  std::string scratch = parameters.s("scratch");
  std::string valence_file_format = parameters.s("valence_file_format");
  std::string occ_file = parameters.s("occ_file");
  std::string density_file = parameters.s("density_file");
  std::string me1j_file = parameters.s("me1j_file");
  std::string me2jp_file = parameters.s("me2jp_file");
  std::string physical_system = parameters.s("physical_system");
  std::string denominator_partitioning = parameters.s("denominator_partitioning");
  std::string NAT_order = parameters.s("NAT_order");

  bool approx_3f2 = parameters.s("approx_3f2") == "true";
  bool use_brueckner_bch = parameters.s("use_brueckner_bch") == "true";
  bool nucleon_mass_correction = parameters.s("nucleon_mass_correction") == "true";
  bool relativistic_correction = parameters.s("relativistic_correction") == "true";
  bool IMSRG3 = parameters.s("IMSRG3") == "true";
  bool imsrg3_n7 = parameters.s("imsrg3_n7") == "true";
  bool imsrg3_mp4 = parameters.s("imsrg3_mp4") == "true";
  bool imsrg3_at_end = parameters.s("imsrg3_at_end") == "true";
  bool imsrg3_no_qqq = parameters.s("imsrg3_no_qqq") == "true";
  bool write_omega = parameters.s("write_omega") == "true";
  bool freeze_occupations = parameters.s("freeze_occupations")=="true";
  bool discard_no2b_from_3n = parameters.s("discard_no2b_from_3n")=="true";
  bool hunter_gatherer = parameters.s("hunter_gatherer") == "true";
  bool goose_tank = parameters.s("goose_tank") == "true";
  bool discard_residual_input3N = parameters.s("discard_residual_input3N")=="true";
  bool use_NAT_occupations = (parameters.s("use_NAT_occupations")=="true") ? true : false;
  bool order_NAT_by_energy = (parameters.s("order_NAT_by_energy")=="true") ? true : false;
  bool store_3bme_pn = (parameters.s("store_3bme_pn")=="true");
  bool only_2b_eta = (parameters.s("only_2b_eta")=="true");
  bool only_2b_omega = (parameters.s("only_2b_omega")=="true");
  bool perturbative_triples = (parameters.s("perturbative_triples")=="true");
  bool brueckner_restart = false;
  bool write_HO_ops = parameters.s("write_HO_ops") == "true";  // added by Antoine Belley
  bool write_HF_ops = parameters.s("write_HF_ops") == "true";  // added by Antoine Belley

  int eMax = parameters.i("emax");
  int lmax = parameters.i("lmax"); // so far I only use this with atomic systems.
  int E3max = parameters.i("e3max");
  int lmax3 = parameters.i("lmax3");
  int targetMass = parameters.i("A");
  int nsteps = parameters.i("nsteps");
  int file2e1max = parameters.i("file2e1max");
  int file2e2max = parameters.i("file2e2max");
  int file2lmax = parameters.i("file2lmax");
  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");
  int atomicZ = parameters.i("atomicZ");
  int emax_unocc = parameters.i("emax_unocc");
  int eMax_imsrg = parameters.i("emax_imsrg");
  int e2Max_imsrg = parameters.i("e2max_imsrg");
  int e3Max_imsrg = parameters.i("e3max_imsrg");
  int eMax_3body_imsrg = parameters.i("emax_3body_imsrg");
//  if ( not ( eMax_imsrg==-1 and e2Max_imsrg==-1 and e3Max_imsrg==-1 ) )
//  {
//    if ( eMax_imsrg==-1 ) eMax_imsrg = eMax;
//    if ( e2Max_imsrg==-1 ) e2Max_imsrg = 2*eMax_imsrg;
//    if ( e3Max_imsrg==-1 ) e3Max_imsrg = std::min( E3max, 3*eMax_imsrg);
//  }
////  if (e2Max_imsrg==-1 and eMax_imsrg != -1) e2Max_imsrg = 2*eMax_imsrg;
////  if (e3Max_imsrg==-1 and eMax_imsrg != -1) e3Max_imsrg = std::min(E3max, 3*eMax_imsrg);

  double hw = parameters.d("hw");
  double smax = parameters.d("smax");
  double ode_tolerance = parameters.d("ode_tolerance");
  double dsmax = parameters.d("dsmax");
  double ds_0 = parameters.d("ds_0");
  double domega = parameters.d("domega");
  double omega_norm_max = parameters.d("omega_norm_max");
  double denominator_delta = parameters.d("denominator_delta");
  double BetaCM = parameters.d("BetaCM");
  double hwBetaCM = parameters.d("hwBetaCM");
  double eta_criterion = parameters.d("eta_criterion");
  double hw_trap = parameters.d("hw_trap");
  double dE3max = parameters.d("dE3max");
  double OccNat3Cut = parameters.d("OccNat3Cut");
  double threebody_threshold = parameters.d("threebody_threshold");

  std::vector<std::string> opnames = parameters.v("Operators");
  std::vector<std::string> opsfromfile = parameters.v("OperatorsFromFile");
  std::vector<std::string> opnamesPT1 = parameters.v("OperatorsPT1");
  std::vector<std::string> opnamesRPA = parameters.v("OperatorsRPA");
  std::vector<std::string> opnamesTDA = parameters.v("OperatorsTDA");

  std::vector<Operator> ops;
  std::vector<std::string> spwf = parameters.v("SPWF");

  using PhysConst::PROTON_RCH2;
  using PhysConst::NEUTRON_RCH2;
  using PhysConst::DARWIN_FOLDY;

  // test 2bme file
  if (inputtbme != "none" and fmt2.find("oakridge")==std::string::npos and fmt2 != "schematic" )
  {
    if( not std::ifstream(inputtbme).good() )
    {
      std::cout << "trouble reading " << inputtbme << "  fmt2 = " << fmt2 << "   exiting. " << std::endl;
      return 1;
    }
  }
  // test 3bme file
  if (input3bme != "none")
  {
    if( not std::ifstream(input3bme).good() )
    {
      std::cout << "trouble reading " << input3bme << " exiting. " << std::endl;
      return 1;
    }
  }


  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  rw.Set3NFormat( fmt3 );

  // Test whether the scratch directory exists and we can write to it.
  // This is necessary because otherwise you get garbage for transformed operators and it's
  // not obvious what went wrong.
  if ( (method == "magnus") and  ( (opnames.size() + opsfromfile.size()) > 0 )  )
  {
    if ( scratch=="/dev/null" or scratch=="/dev/null/")
    {
      std::cout << "ERROR!!! using Magnus with scratch = " << scratch << " but you're also trying to transform some operators. Dying now. " << std::endl;
      exit(EXIT_FAILURE);
    }
    else if ( scratch != "" )
    {
      std::string testfilename = scratch + "/_this_is_a_test_delete_me";
      std::ofstream testout(testfilename);
      testout << "PASSED" << std::endl;
      testout.close();

      // now read it back.
      std::ifstream testin(testfilename);
      std::string checkpassed;
      testin >> checkpassed;
      if ( (checkpassed != "PASSED") or ( not testout.good() ) or ( not testin.good() ) )
      {
        std::cout << "ERROR in " << __FILE__ <<  " failed test write to scratch directory " << scratch << " that's bad. Dying now." << std::endl;
        exit(EXIT_FAILURE);
      }

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


  ModelSpace modelspace = ( reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space) );

  modelspace.SetE3max(E3max);
  modelspace.SetLmax(lmax);
  modelspace.SetdE3max(dE3max);
  modelspace.SetOccNat3Cut(OccNat3Cut);



  if (nsteps < 0) // default to 1 step for single ref, 2 steps for valence decoupling
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);

  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();
  std::cout << "Reading interactions..." << std::endl;


  if (inputtbme != "none")
  {
    rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
    std::cout << "done reading 2N" << std::endl;
  }

  // Read in the 3-body file
  if (Hbare.particle_rank >=3)
  {
    if(input3bme_type == "full")
    {
      rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
    }
    if(input3bme_type == "no2b")
    {
      Hbare.ThreeBody.SetMode("no2b");
      if (no2b_precision == "half")  Hbare.ThreeBody.SetMode("no2bhalf");

      Hbare.ThreeBody.ReadFile( {input3bme}, {file3e1max, file3e2max, file3e3max, file3e1max} );
      rw.File3N = input3bme;
    }
    std::cout << "done reading 3N" << std::endl;
  }

  if (store_3bme_pn)
  {
    Hbare.ThreeBody.TransformToPN();
  }

  Hbare += imsrg_util::Trel_Op(modelspace);

  // Add a Lawson center of mass term. If hwBetaCM is specified, use that frequency, otherwise use the basis frequency
  if (std::abs(BetaCM)>1e-6)
  {
    if (hwBetaCM < 0) hwBetaCM = modelspace.GetHbarOmega();
    std::ostringstream hcm_opname;
    hcm_opname << "HCM_" << hwBetaCM;
    Hbare += BetaCM * imsrg_util::OperatorFromString( modelspace, hcm_opname.str());
  }

  std::cout << "Creating HF" << std::endl;
  HFMBPT hf(Hbare); // HFMBPT inherits from HartreeFock, so this works for HF and NAT bases.
  Operator rho = Operator(modelspace,0,0,0,1);
  Operator& HNO = Hbare; // The reference & means we overwrite Hbare and save some memory
  int hno_particle_rank = 2;

  rw.Read_me1j(density_file, rho, eMax, eMax);
  if(me1j_file != "none" and me2jp_file != "none"){
    HNO.Erase();
    rw.Read_me1j(me1j_file, HNO, eMax, eMax);
    rw.Read_me2jp(me2jp_file, HNO, eMax, 2*eMax, eMax);
    //std::map<std::array<int, 4>,double> hole_map;
    //for ( auto i : modelspace.all_orbits) {
    //  Orbit& oi = modelspace.GetOrbit(i);
    //  hole_map[{oi.n, oi.l, oi.j2, oi.tz2}] = rho.OneBody(i);
    //}
    //modelspace.Init(eMax, hole_map, valence_space);
  }
  else {
    arma::mat C;
    arma::vec Occ;
    bool success = false;
    success = arma::eig_sym(Occ, C, rho.OneBody); // eigenvalues of rho (i.e. occupations) are in ascending order
    Occ = arma::reverse(Occ);
    C = arma::reverse(C, 1);
    std::map<std::array<int, 4>,double> hole_map;
    for ( auto i : modelspace.all_orbits) {
      Orbit& oi = modelspace.GetOrbit(i);
      hole_map[{oi.n, oi.l, oi.j2, oi.tz2}] = Occ(i);
    }
    modelspace.Init(eMax, hole_map, valence_space);
    hf.Occ = Occ;

    hf.rho = rho.OneBody;
    hf.BuildMonopoleV();
    if(Hbare.GetParticleRank()>2) hf.BuildMonopoleV3();

    hf.DiagonalizeRho();
    if(NAT_order=="close_to_1"){
      std::cout << "Ordering NAT orbits so that the transformation close to 1..." << std::endl;
      for (auto& it : Hbare.OneBodyChannels){
        arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end()));
        arma::mat CNAT_chan = hf.C_HF2NAT.submat(orbvec, orbvec);
        arma::mat tmp = hf.C_HF2NAT.submat(orbvec, orbvec);
        arma::vec occ = hf.Occ.elem(orbvec);
        arma::vec tmp_occ = hf.Occ.elem(orbvec);
        for(int i=0; i<CNAT_chan.n_rows; i++){
          arma::rowvec v = CNAT_chan.row(i);
          arma::uword idx = arma::index_max(arma::abs(v));
          tmp.col(i) = CNAT_chan.col(idx);
          tmp_occ(i) = occ(idx);
        }
        hf.C_HF2NAT.submat(orbvec, orbvec) = tmp;
        hf.Occ.elem(orbvec) = tmp_occ;
      }
    }
    hf.C_HO2NAT = hf.C_HF2NAT;
    HNO = hf.GetNormalOrderedH(hf.C_HO2NAT, hno_particle_rank);
  }
  for (auto& it : HNO.OneBodyChannels){
    arma::uvec orbvec(std::vector<index_t>(it.second.begin(),it.second.end()));
    std::cout << "F: " << std::endl;
    arma::mat tmp = HNO.OneBody.submat(orbvec, orbvec);
    std::cout << tmp << std::endl;
    std::cout << "rho: " << std::endl;
    tmp = rho.OneBody.submat(orbvec, orbvec);
    std::cout << tmp << std::endl;
  }
  std::cout << basis << " Single particle energies and wave functions:" << std::endl;
  std::cout << std::fixed << std::setw(3) << "i" << ": " << std::setw(3) << "n" << " " << std::setw(3) << "l" << " "
       << std::setw(3) << "2j" << " " << std::setw(3) << "2tz" << "   " << std::setw(12) << "SPE" << " " << std::setw(12) << "occ."
       << " " << std::setw(12) << "occNAT" << "   |   " << " overlaps" << std::endl;
  for ( auto i : modelspace.all_orbits )
  {
    Orbit& oi = modelspace.GetOrbit(i);
    std::cout << std::fixed << std::setw(3) << i << ": " << std::setw(3) << oi.n << " " << std::setw(3) << oi.l << " "
         << std::setw(3) << oi.j2 << " " << std::setw(3) << oi.tz2 << "   " << std::setw(12) << std::setprecision(6) << HNO.OneBody(i,i) << " " << std::setw(12) << oi.occ << " " << std::setw(12) << oi.occ_nat << "   | ";
    for (int j : Hbare.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) // j runs over HO states
    {
      std::cout << std::setw(9) << hf.C_HO2NAT(j,i) << "  ";  // C is <HO|NAT>
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  HNO -= BetaCM * 1.5*hwBetaCM; // This is just the zero-body piece. The other stuff was added earlier.
  std::cout << "Hbare 0b = " << std::setprecision(8) << HNO.ZeroBody << std::endl;
  ModelSpace modelspace_imsrg = modelspace;

  std::cout << "done with perterbative stuff, method = " << method << std::endl;
  // Calculate all the desired operators. If we're using magnus, we'll do this after the flow is over
  if ( method != "magnus" )
  {

    for (auto& opname : opnames)
    {
        ops.emplace_back( imsrg_util::OperatorFromString(modelspace, opname) );
    }

    // Calculate first order perturbative correction to some operators, if that's what we asked for.
    // Strictly speaking, it doesn't make much sense to do this and then proceed with the IMSRG calculation,
    // but I'm not here to tell people what to do...
    for (auto& opnamept1 : opnamesPT1 )
    {
      ops.emplace_back( imsrg_util::FirstOrderCorr_1b( imsrg_util::OperatorFromString(modelspace,opnamept1)   , HNO ) );
      opnames.push_back( opnamept1+"PT1" );
    }
    for (auto& opnametda : opnamesTDA )
    {  // passing the argument "TDA" just sets the phhp and hpph blocks to zero in the RPA calculation
      ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnametda)   , HNO, "TDA" ) );
      opnames.push_back( opnametda+"TDA" );
    }
    for (auto& opnamerpa : opnamesRPA )
    {
      ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnamerpa)   , HNO, "RPA" ) );
      opnames.push_back( opnamerpa+"RPA" );
    }


   if (ops.size()>0)
   {
     std::cout << "operators to transform: " << std::endl;
     for ( auto& opn : opnames ) std::cout << opn << " ";
     std::cout << std::endl;
   }

//  for (auto& op : ops)
   for (size_t i=0;i<ops.size();++i)
   {
     if (ops[i].GetJRank()==0 and (ops[i].GetTRank()!=0 or ops[i].GetParity()!=0) )
     {
         std::cout << "Before doing HF transformation, making op " << i << " " << opnames[i] << " not reduced. " << std::endl;
         ops[i].MakeNotReduced();
     }
//     std::cout << "Before transforming  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
      // We don't transform a DaggerHF, because we want the a^dagger to already refer to the HF basis.
     if ((basis == "HF") and (opnames[i].find("DaggerHF") == std::string::npos)  )
     {
       ops[i] = hf.TransformToHFBasis(ops[i]);
     }
     else if ((basis == "NAT") and (opnames[i].find("DaggerHF") == std::string::npos)  )
     {
       ops[i] = hf.TransformHOToNATBasis(ops[i]);
     }
//     std::cout << "After transforming  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
     ops[i] = ops[i].DoNormalOrdering();
//     std::cout << "Before normal ordering  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
       std::cout << basis << " expectation value  " << opnames[i] << "  " << ops[i].ZeroBody << std::endl;
     if (method == "MP3")
     {
       double dop = ops[i].MP1_Eval( HNO );
       std::cout << "Operator 1st order correction  " << dop << "  ->  " << ops[i].ZeroBody + dop << std::endl;
     }
    if ( opnames[i] == "Rp2" )
    {
      double Rp2 = ops[i].ZeroBody;
      int Z = modelspace.GetTargetZ();
      int A = modelspace.GetTargetMass();
      std::cout << " HF point proton radius = " << sqrt( Rp2 ) << std::endl;
      std::cout << " HF charge radius = " << ( abs(Rp2)<1e-6 ? 0.0 : sqrt( Rp2 + PROTON_RCH2 + NEUTRON_RCH2*(A-Z)/Z + DARWIN_FOLDY) ) << std::endl;
    }
   }// for ops.size


  }// if method != "magnus"



  if (method == "FCI")
  {
   if ( valence_file_format == "tokyo" )
   {
      HNO = HNO.UndoNormalOrdering();
      for (size_t i=0; i<ops.size();i++)
      {
         ops[i] = ops[i].UndoNormalOrdering();
      }


      modelspace.SetReference("vacuum");
      rw.WriteTokyo(HNO,intfile+".snt", "");
      // Haven't yet implemented FCI operators for Tokyo format. I should do this...
      for (size_t i=0; i<ops.size();i++)
      {
         if (ops[i].GetJRank()==0 and ops[i].GetTRank()==0 )
         {
           rw.WriteTokyo(ops[i], intfile + "_" + opnames[i] + ".snt","");
         }
         else
         {
          rw.WriteTensorTokyo(intfile+opnames[i]+".snt",ops[i]);
         }
      }
   }
   else // Write in NuShellX Format
   {
     // we want the 1b piece to be diagonal in the vacuum NO representation
      HNO = HNO.UndoNormalOrdering();
      double previous_zero_body = HNO.ZeroBody;
      modelspace.SetReference("vacuum");
      HartreeFock hfvac(HNO);
      hfvac.Solve();

  //    Operator Hvac = hfvac.GetNormalOrderedH();
      HNO = hfvac.GetNormalOrderedH();
      std::cout << "HNO had zero body = " << HNO.ZeroBody << "  and I add " << previous_zero_body << " to it. " << std::endl;
      HNO.ZeroBody += previous_zero_body;

      rw.WriteNuShellX_int(HNO,intfile+".int");
      rw.WriteNuShellX_sps(HNO,intfile+".sp");

      std::cout << "NO wrt vacuum. One Body term is hopfully still diagonal?" << std::endl << HNO.OneBody << std::endl;

      for (index_t i=0;i<ops.size();++i)
      {
        ops[i] = ops[i].UndoNormalOrdering();
        if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
        {
          rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
        }
        else
        {
          rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
          rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
        }
      }
    }
    HNO.PrintTimes();
    return 0;
  }
  if (only_2b_omega)
  {
    std::cout << " Restricting the Magnus operator Omega to be 2b." << std::endl;
    BCH::SetOnly2bOmega(only_2b_omega);
  }


  // We may want to use a smaller model space for the IMSRG evolution than we used for the HF step.
  // This is most effective when using natural orbitals or when including 3-body operators.
//  ModelSpace modelspace_imsrg = ( reference=="default" ? ModelSpace(eMax_imsrg,e2Max_imsrg,e3Max_imsrg,valence_space) : ModelSpace(eMax_imsrg,e2Max_imsrg,e3Max_imsrg,reference,valence_space) );
//  ModelSpace modelspace_imsrg = modelspace;
  if ( (eMax_imsrg != -1) or (e2Max_imsrg != -1) or (e3Max_imsrg != -1) or (eMax_3body_imsrg != -1))
  {

//     if ( eMax_imsrg==-1 ) eMax_imsrg = eMax;
//     if ( e2Max_imsrg==-1 ) e2Max_imsrg = 2*eMax_imsrg;
//     if ( e3Max_imsrg==-1 ) e3Max_imsrg = std::min( E3max, 3*eMax_imsrg);
//     if ( eMax_3body_imsrg==-1) eMax_3body_imsrg = eMax_imsrg;
//
////     ModelSpace modelspace_imsrg = modelspace;
//     std::cout << "Truncating modelspace for IMSRG calculation: emax e2max e3max  ->  " << eMax_imsrg << " " << e2Max_imsrg << " " << e3Max_imsrg << std::endl;
//     modelspace_imsrg.SetEmax( eMax_imsrg);
//     modelspace_imsrg.SetE2max( e2Max_imsrg);
//     modelspace_imsrg.SetE3max( e3Max_imsrg);
//     modelspace_imsrg.SetEmax3Body( eMax_3body_imsrg );
//     modelspace_imsrg.Init( eMax_imsrg, reference, valence_space);
//   //  if (emax_unocc>0) modelspace_imsrg.SetEmaxUnocc(emax_unocc);
//     if (physical_system == "atomic") modelspace_imsrg.InitSingleSpecies(eMax_imsrg, reference, valence_space);
//     if (occ_file != "none" and occ_file != "" ) modelspace_imsrg.Init_occ_from_file(eMax_imsrg,valence_space,occ_file);
////     if (physical_system == "atomic") modelspace_imsrg.InitSingleSpecies(eMax_imsrg, eMax_imsrg, e3Max_imsrg, reference, valence_space);
////     if (occ_file != "none" and occ_file != "" ) modelspace_imsrg.Init_occ_from_file(eMax_imsrg,e2Max_imsrg,e3Max_imsrg,valence_space,occ_file);
//
//
//     // If the occupations in modelspace were different from the naive filling, we want to keep those.
//     std::map<index_t,double> hole_map;
//     for ( auto& i_new : modelspace_imsrg.all_orbits )
//     {
//        Orbit& oi_new = modelspace_imsrg.GetOrbit(i_new);
//        index_t i_old = modelspace.GetOrbitIndex( oi_new.n, oi_new.l, oi_new.j2, oi_new.tz2 );
//        Orbit& oi_old = modelspace.GetOrbit(i_old);
//        hole_map[i_new] = oi_old.occ;
//     }
//     modelspace_imsrg.SetReference( hole_map );

     /// If HNO has a 3N piece, we already did the truncation while transforming to the HF basis
     /// so we don't want to do that again. Kludgey solution, make a temporary 2N operator, truncate and copy.
     if (HNO.GetParticleRank() < 3)
     {
       HNO = HNO.Truncate(modelspace_imsrg);
       if (IMSRG3) // we'll want a 3N structure for IMSRG3
       {
           // Always do IMSRG(3) in pn mode. SetMode also calls Allocate.
           HNO.ThreeBody.SetMode("pn");
           HNO.SetParticleRank(3);
       }
     }
     else
     {
       Operator Htmp2b = Operator(modelspace, 0,0,0,2);
       Htmp2b.OneBody = HNO.OneBody;
       Htmp2b.TwoBody = HNO.TwoBody;
       Htmp2b = Htmp2b.Truncate(modelspace_imsrg);
       HNO.OneBody = Htmp2b.OneBody;
       HNO.TwoBody = Htmp2b.TwoBody;
     }

//     HNO = HNO.Truncate(modelspace_imsrg);
//     if (IMSRG3) {
//       HNO.ThreeBody.SwitchToPN_and_discard();
//     }

//     modelspace = modelspace_imsrg;  // this could cause some confusion later on...
//    hf.PrintSPEandWF();
  }
  else
  {
    std::cout << "Im here " << __LINE__ << " particle rank is " << HNO.GetParticleRank() << " IMSRG3 is " << IMSRG3 << std::endl;
    HNO.SetModelSpace(modelspace_imsrg);
    if (HNO.GetParticleRank()<3 and IMSRG3) {
        HNO.ThreeBody.SetMode("pn");
        HNO.SetParticleRank(3);
    std::cout << "Im here " << __LINE__ << " particle rank is " << HNO.GetParticleRank() << "  pn mode? " << HNO.ThreeBody.Is_PN_Mode() << std::endl;
//      HNO.ThreeBody.SwitchToPN_and_discard();
    }
  }

 // After truncating, get the perturbative energies again to see how much things changed.
  if (eMax_imsrg != eMax)
  {
    std::cout << "Perturbative estimates of gs energy:" << std::endl;
    double EMP2 = HNO.GetMP2_Energy();
    double EMP2_3B = HNO.GetMP2_3BEnergy();
    std::cout << "EMP2 = " << EMP2 << std::endl;
    std::cout << "EMP2_3B = " << EMP2_3B << std::endl;
    std::cout << "To 2nd order, E = " << HNO.ZeroBody + EMP2 + EMP2_3B << std::endl;
    std::array<double,3> Emp_3 = HNO.GetMP3_Energy();
    double EMP3 = Emp_3[0]+Emp_3[1]+Emp_3[2];
    std::cout << "E3_pp = " << Emp_3[0] << "  E3_hh = " << Emp_3[1] << " E3_ph = " << Emp_3[2] << "   EMP3 = " << EMP3 << std::endl;
    std::cout << "To 3rd order, E = " << HNO.ZeroBody + EMP2 + EMP3 + EMP2_3B << std::endl;
  }

  if ( method == "MP3" )
  {
    HNO.PrintTimes();
    return 0;
  }


//// Now we're ready do to the IMSRG calculation.

  IMSRGSolver imsrgsolver(HNO);
//  imsrgsolver.SetHin(HNO); // necessary?
  imsrgsolver.SetReadWrite(rw);
  imsrgsolver.SetMethod(method);
  imsrgsolver.SetDenominatorPartitioning(denominator_partitioning);
  imsrgsolver.SetEtaCriterion(eta_criterion);
  imsrgsolver.GetGenerator().SetOnly2bEta(only_2b_eta);
  imsrgsolver.max_omega_written = 500;
  imsrgsolver.SetHunterGatherer( hunter_gatherer );
  imsrgsolver.SetPerturbativeTriples(perturbative_triples);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(dsmax);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if(approx_3f2) {
    imsrgsolver.SetHunterGatherer(true);
    BCH::SetUseFactorizedCorrection(true);
    Commutator::FactorizedDoubleCommutator::SetUse_1b_Intermediates(true);
    Commutator::FactorizedDoubleCommutator::SetUse_2b_Intermediates(false);
  }
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  BCH::SetUseBruecknerBCH(use_brueckner_bch);
  Commutator::SetUseIMSRG3(IMSRG3);
  Commutator::SetUseIMSRG3N7(imsrg3_n7);
  Commutator::SetUseIMSRG3_MP4(imsrg3_mp4);
  Commutator::SetIMSRG3Noqqq(imsrg3_no_qqq);
  if (use_brueckner_bch)
  {
    std::cout << "Using Brueckner flavor of BCH" << std::endl;
  }
  if (IMSRG3)
  {
    std::cout << "Using IMSRG(3) commutators. This will probably be slow..." << std::endl;
  }
  if (imsrg3_n7)
  {
    std::cout << "  only including IMSRG3 commutator terms that scale up to n7" << std::endl;
  }
  if ( threebody_threshold > 1e-12 )
  {
    std::cout << "skipping IMSRG(3) commutator terms if norm of either operator is below " << threebody_threshold << std::endl;
  }


  // Here's where we need to have the operators
//  std::cout << "MADE IT TO LINE " << __LINE__ << std::endl;


  if (method == "flow" or method == "flow_RK4" )
  {
    for (auto& op : ops )  imsrgsolver.AddOperator( op );
    std::cout << " Added ops. FlowingOps.size = " << imsrgsolver.FlowingOps.size() << std::endl;
  }

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

  if (IMSRG3)
  {
    std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
  }
  if ( perturbative_triples and method=="magnus" )
  {
//    modelspace.SetdE3max(dE3max);
//    modelspace.SetOccNat3Cut(OccNat3Cut);
//    size_t nstates_kept = modelspace.CountThreeBodyStatesInsideCut();
//    std::array<size_t,2> nstates = modelspace.CountThreeBodyStatesInsideCut();
//    std::cout << "Truncations: dE3max = " << dE3max << "   OccNat3Cut = " << std::scientific << OccNat3Cut << "  ->  number of 3-body states kept:  " << nstates[0] << " out of " << nstates[1] << std::endl << std::fixed;
//    double dE_triples = imsrgsolver.GetPerturbativeTriples();
    double dE_triples = imsrgsolver.CalculatePerturbativeTriples();
    std::cout << "Perturbative triples:  " << std::setw(16) << std::setprecision(8) << dE_triples << " -> " << imsrgsolver.GetH_s().ZeroBody + dE_triples << std::endl;
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
//    modelspace.ResetFirstPass();
    modelspace_imsrg.ResetFirstPass();
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


  if ( imsrg3_at_end )
  {
    if ( method.find("magnus") != std::string::npos )
    {
      std::cout << "Performing final BCH transformation at the IMSRG(3) level" << std::endl;

//      modelspace.SetdE3max(dE3max);
//      modelspace.SetOccNat3Cut(OccNat3Cut);
//      int new_E3max = std::min(modelspace.GetE3max(), int( std::ceil(3*std::max( modelspace.GetEFermi()[-1], modelspace.GetEFermi()[+1])+dE3max)));
      modelspace_imsrg.SetdE3max(dE3max);
      modelspace_imsrg.SetOccNat3Cut(OccNat3Cut);
      int new_E3max = std::min(modelspace_imsrg.GetE3max(), int( std::ceil(3*std::max( modelspace_imsrg.GetEFermi()[-1], modelspace_imsrg.GetEFermi()[+1])+dE3max)));
      std::cout << "Setting new E3max = " << new_E3max << std::endl;
//      modelspace.SetE3max(  new_E3max);
      modelspace_imsrg.SetE3max(  new_E3max);

//      Operator H3(modelspace,0,0,0,3);
      Operator H3(modelspace_imsrg,0,0,0,3);
      std::cout << "Constructed H3" << std::endl;
      H3.ZeroBody = HNO.ZeroBody;
      H3.OneBody = HNO.OneBody;
      H3.TwoBody = HNO.TwoBody;
      HNO = H3;
      std::cout << "Replacing HNO" << std::endl;
      std::cout << "Hbare Three Body Norm is " << Hbare.ThreeBodyNorm() << std::endl;
      HNO.ThreeBody.SwitchToPN_and_discard();


      Commutator::SetUseIMSRG3(true);
      Commutator::SetUseIMSRG3N7(imsrg3_n7);
//      Operator H_with_3 = imsrgsolver.Transform(  *(imsrgsolver.H_0) );
      Operator H_with_3 = imsrgsolver.Transform( HNO );

      // Now throw away the residual 3-body so we don't need to keep it after re-normal ordering
      H_with_3.SetNumberLegs(4);
      H_with_3.SetParticleRank(2);
      imsrgsolver.FlowingOps[0] = H_with_3;
    }
    else
    {
      std::cout << "selected imsrg3_at_end, but method != magnus, so I don't know what to do. Ignoring." << std::endl;
    }

  }


/*
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
*/


  // If we're doing targeted/ensemble normal ordering
  // we now re-normal order wrt to the core
  // and do any remaining flow.
//  ModelSpace ms2(modelspace);
  ModelSpace ms2(modelspace_imsrg);
  ms2.SetReference(ms2.core); // change the reference
  bool renormal_order = false;
  Operator Hs = imsrgsolver.GetH_s();

  Hs = Hs.UndoNormalOrdering();
  //hf.C_HO2NAT = hf.C_HO2NAT.t(); // NAT --> HO
  //Hs = hf.TransformHOToNATBasis(Hs);
  rw.Write_me1j(me1j_file + ".out", Hs, Hs.modelspace->GetEmax(), Hs.modelspace->GetLmax());
  rw.Write_me2jp(me2jp_file + ".out", Hs, Hs.modelspace->GetEmax(), Hs.modelspace->GetE2max(), Hs.modelspace->GetLmax());
  //hf.C_HO2NAT = hf.C_HO2NAT.t(); // HO --> NAT


//  if (modelspace.valence.size() > 0 )
  if (modelspace_imsrg.valence.size() > 0 )
//  if (modelspace.valence.size() > 0 or basis=="NAT")
  {
//    renormal_order = modelspace.holes.size() != modelspace.core.size();
    renormal_order = modelspace.holes.size() != modelspace_imsrg.core.size();
    if (not renormal_order)
    {
//      for (auto c : modelspace.core)
      for (auto c : modelspace_imsrg.core)
      {
//         if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (std::abs(1-modelspace.GetOrbit(c).occ)>1e-6))
         if ( (find( modelspace_imsrg.holes.begin(), modelspace_imsrg.holes.end(), c) == modelspace_imsrg.holes.end()) or (std::abs(1-modelspace_imsrg.GetOrbit(c).occ)>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
  if ( renormal_order )
  {
    Operator Hs = HNO;
    if(approx_3f2){
      Commutator::FactorizedDoubleCommutator::SetUse_1b_Intermediates(true);
      Commutator::FactorizedDoubleCommutator::SetUse_2b_Intermediates(true);
      Hs = imsrgsolver.Transform(HNO);
      double dE_triple = imsrgsolver.CalculatePerturbativeTriples();
      std::cout << "Perturbative triples: " << std::setw(16) << std::setprecision(8) << dE_triple << std::endl;
      Hs.ZeroBody += dE_triple;
    }
    else{
      Hs = imsrgsolver.GetH_s();
    }



//    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
//    std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
    std::cout << "Undoing NO wrt A=" << modelspace_imsrg.GetAref() << " Z=" << modelspace_imsrg.GetZref() << std::endl;
    std::cout << "Before doing so, the spes are " << std::endl;
//    for ( auto i : modelspace.all_orbits ) std::cout << "  " << i << " : " << Hs.OneBody(i,i) << std::endl;
    for ( auto i : modelspace_imsrg.all_orbits ) std::cout << "  " << i << " : " << Hs.OneBody(i,i) << std::endl;
    if (IMSRG3)
    {
      std::cout << "Re-normal-ordering wrt the core. For now, we just throw away the 3N at this step." << std::endl;
      Hs.SetNumberLegs(4);
      Hs.SetParticleRank(2);
    }

    Hs = Hs.UndoNormalOrdering();
    Hs.SetModelSpace(ms2);
    std::cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << std::endl;
    Hs = Hs.DoNormalOrdering();


    imsrgsolver.FlowingOps[0] = Hs;


// More flowing is unnecessary, since things should stay decoupled.
//    imsrgsolver.SetHin(HNO);
//    imsrgsolver.SetEtaCriterion(1e-4);
//    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
//    std::cout << "Final transformation on the operators..." << std::endl;
//    int iop = 0;
//    for (auto& op : ops)
//    {
//      std::cout << opnames[iop++] << std::endl;
//      op = op.UndoNormalOrdering();
//      op.SetModelSpace(ms2);
//      op = op.DoNormalOrdering();
//      // transform using the remaining omegas
//      op = imsrgsolver.Transform_Partial(op,nOmega);
//    }
  }


  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
//  if (modelspace.valence.size() > 0)
  if (modelspace_imsrg.valence.size() > 0)
  {
    if (valence_file_format == "antoine") // this is still being tested...
    {
      rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
//      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace_imsrg.GetAref(),modelspace_imsrg.GetZref());
    }
    std::cout << "Writing files: " << intfile << std::endl;
    if (valence_file_format == "tokyo")
    {
     rw.WriteTokyo(imsrgsolver.GetH_s(),intfile+".snt", "");
    }
    else
    {
      rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
      rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
    }

//    if (method == "magnus" or method=="flow_RK4")
//    {
//       for (index_t i=0;i<ops.size();++i)
//       {
//          if ( ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1) and (ops[i].GetNumberLegs()%2==0) )
//          {
//            if (valence_file_format == "tokyo")
//            {
//              rw.WriteTokyo(ops[i],intfile+opnames[i]+".snt", "op");
//            }
//            else
//            {
//              rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
//            }
//          }
//          else if ( ops[i].GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
//          {
////            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
//            rw.WriteDaggerOperator( ops[i], intfile+opnames[i]+".dag",opnames[i]);
//          }
//          else
//          {
//            if (valence_file_format == "tokyo")
//            {
//              rw.WriteTensorTokyo(intfile+opnames[i]+"_2b.snt",ops[i]);
//            }
//            else
//            {
//              rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
//              rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
//            }
//          }
//       }
//    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    double dE_triple = 0.0;
    if(approx_3f2){
      Commutator::FactorizedDoubleCommutator::SetUse_1b_Intermediates(true);
      Commutator::FactorizedDoubleCommutator::SetUse_2b_Intermediates(true);
      double dE_triple = imsrgsolver.CalculatePerturbativeTriples();
    }
    std::cout << "Core Energy = " << std::setprecision(6) << imsrgsolver.GetH_s().ZeroBody + dE_triple << std::endl;
    if ( method != "magnus")
    {
      for (index_t i=0;i<ops.size();++i)
      {
//        Operator& op = ops[i];
        Operator& op = imsrgsolver.FlowingOps[i+1]; // the first operator is the Hamiltonian
        std::cout << opnames[i] << " = " << op.ZeroBody << std::endl;
        if ( opnames[i] == "Rp2" )
        {
//           int Z = modelspace.GetTargetZ();
//           int A = modelspace.GetTargetMass();
           int Z = modelspace_imsrg.GetTargetZ();
           int A = modelspace_imsrg.GetTargetMass();
           std::cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
           std::cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + PROTON_RCH2 + NEUTRON_RCH2*(A-Z)/Z + DARWIN_FOLDY) << std::endl;
        }
        if ((op.GetJRank()>0) or (op.GetTRank()>0)) // if it's a tensor, you probably want the full operator
        {
          std::cout << "Writing operator to " << intfile+opnames[i]+".op" << std::endl;
          rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
        }
      }
    }
  }




  Hbare.PrintTimes();

  return 0;
}

