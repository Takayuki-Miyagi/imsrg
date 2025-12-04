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

  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");
  std::string custom_valence_space = parameters.s("custom_valence_space");
  std::string basis = parameters.s("basis");
  std::string intfile = parameters.s("intfile");
  std::string input_op_fmt = parameters.s("input_op_fmt");
  std::string scratch = parameters.s("scratch");
  std::string valence_file_format = parameters.s("valence_file_format");
  std::string occ_file = parameters.s("occ_file");

  bool approx_3f2 = parameters.s("approx_3f2") == "true";
  bool write_omega = parameters.s("write_omega") == "true";

  int eMax = parameters.i("emax");
  int lmax = parameters.i("lmax"); // so far I only use this with atomic systems.
  int E3max = parameters.i("e3max");
  int lmax3 = parameters.i("lmax3");
  int targetMass = parameters.i("A");
  int nsteps = parameters.i("nsteps");
  double hw = parameters.d("hw");
  double smax = parameters.d("smax");

  std::vector<std::string> opnames = parameters.v("Operators");
  std::vector<std::string> opsfromfile = parameters.v("OperatorsFromFile");

  std::vector<Operator> ops;

  // unpack the awkward input format for reading an operator from file, and put it into a struct.
  // the format should look like OpName^j_t_p_r^/path/to/2bfile^/path/to/3bfile  if particle rank of Op is 2-body, then 3bfile is not needed.
  std::vector< OpFromFile> opsfromfile_unpacked;
  // If we're reading in other operators, make sure those are ok too
  for (auto& tag : opsfromfile)
  {
     std::istringstream ss(tag);
     std::string opname,qnumbers,f2name,f3name="";

     OpFromFile opff;

     getline(ss,opname,'^');
     getline(ss,qnumbers,'^');
     getline(ss,f2name,'^');
     if ( not ss.eof() )  getline(ss,f3name,'^');
     opff.opname = opname;
     opff.file2name = f2name;
     opff.file3name = f3name;

      ss.str(qnumbers);
      ss.clear();
      std::string tmp;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.j;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.t;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.p;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.r;

      std::cout << "Parsed tag. opname = " << opff.opname << "  " << opff.j << " " << opff.t << " " << opff.p << " " << opff.r << "   file2 = " << opff.file2name   << "    file3 = " << opff.file3name << std::endl;

      // now make sure the files exist before we add them to the list.

//     if( not std::ifstream(f2name).good() )
     if( not std::ifstream(opff.file2name).good() )
     {
//       std::cout << "trouble reading " << f2name << " exiting. " << std::endl;
       std::cout << "trouble reading " << opff.file2name << " exiting. " << std::endl;
       return 1;
     }

     if ( opff.file3name != "") // is there a 3-body file too?
     {
//       getline(ss,f3name,'^');
//       if( not std::ifstream(f3name).good() )
       if( not std::ifstream(opff.file3name).good() )
       {
         std::cout << "trouble reading " << opff.file3name << " exiting. " << std::endl;
//         std::cout << "trouble reading " << f3name << " exiting. " << std::endl;
         return 1;
       }
     }
     // if the files look good, then add it to the list
     opsfromfile_unpacked.push_back( opff );
  }



  ReadWrite rw;
  rw.SetScratchDir(scratch);

  // Test whether the scratch directory exists and we can write to it.
  // This is necessary because otherwise you get garbage for transformed operators and it's
  // not obvious what went wrong.
  if ( opnames.size() + opsfromfile.size() > 0)
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
      intfile = parameters.DefaultIntFile();
    }
    valence_space = custom_valence_space;
  }


  ModelSpace modelspace = ( reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space) );

//  std::cout << __LINE__ << "  constructed modelspace " << std::endl;
  modelspace.SetE3max(E3max);
  modelspace.SetLmax(lmax);

  if (occ_file == "none" or occ_file == "" )
  {
    std::cout << "You need to pass occ_file, exiting." << std::endl;
    return 1;
  }
  modelspace.Init_occ_from_file(eMax,valence_space,occ_file);

  modelspace.SetHbarOmega(hw);
  if (targetMass>0) modelspace.SetTargetMass(targetMass);
  if (lmax3>0) modelspace.SetLmax3(lmax3);

//  std::cout << "Making the Hamiltonian..." << std::endl;
  Operator Hbare = Operator(modelspace,0,0,0,2);
  Hbare.SetHermitian();

  std::cout << "Creating HF" << std::endl;
  HFMBPT hf(Hbare); // HFMBPT inherits from HartreeFock, so this works for HF and NAT bases.
  hf.C.load(intfile + "C.mat");

  IMSRGSolver imsrgsolver(Hbare);
  imsrgsolver.SetReadWrite(rw);

  for ( auto& opff : opsfromfile_unpacked)
  {
    opnames.push_back( opff.opname + "_FROMFILE");
  }

  int count_from_file =0;
  if (opnames.size()>0) std::cout << "transforming operators" << std::endl;

  for (size_t i=0;i<opnames.size();++i)
  {
    auto opname = opnames[i];
    std::cout << i << ": " << opname << " " << std::endl;

    Operator op;

    if ( opname.find("_FROMFILE") != std::string::npos)
    {
      OpFromFile& opff = opsfromfile_unpacked[count_from_file];
      std::cout << "reading " << opff.opname << " with " << opff.j << " " << opff.t << " " << opff.p << " " << opff.r << "  from file " << opff.file2name << std::endl;
      op = Operator(modelspace, opff.j, opff.t, opff.p, opff.r );
      if (opff.r>2) op.ThreeBody.Allocate();
      if ( input_op_fmt == "navratil" )
      {
        rw.Read2bCurrent_Navratil( opff.file2name, op );
      }
      else if ( input_op_fmt == "miyagi" )
      {
        if (opff.file2name != "")
        {
          Operator optmp = rw.ReadOperator2b_Miyagi( opff.file2name, modelspace );
          op.OneBody = optmp.OneBody;
          op.TwoBody = optmp.TwoBody;
        }
//        if ( opff.r>2 and opff.file3name != "")  rw.Read_Darmstadt_3body( opff.file3name, op,  file3e1max,file3e2max,file3e3max);
      }
      count_from_file++;
      opname = opff.opname; // Get rid of the _FROMFILE bit.
    }
    else
    {
      op = imsrg_util::OperatorFromString( modelspace, opname );
    }
    //      Operator op = imsrg_util::OperatorFromString( modelspace, opname );

    if ( op.GetJRank()==0 and ( op.GetTRank()!=0 or op.GetParity()!=0 ) )
    {
      std::cout << "Before doing HF, making " << opname << "  not reduced" << std::endl;
      op.MakeNotReduced();
    }

    op = hf.TransformToHFBasis(op).DoNormalOrdering();
    std::cout << "   HF: " << op.ZeroBody << std::endl;

    int istart =-1;
    int nomega = 0;
    bool checkpoint = true;
    int n_op_written = 0;
    while (true)
    {
      std::string omega_file = intfile + "_Omega_" + std::to_string(nomega);
      std::string file_dump = intfile + "_" + opname + "_" + std::to_string(nomega);
      if(not std::ifstream(omega_file).good()) break;
      if(std::ifstream(file_dump).good()) {
        std::ifstream f(file_dump, std::fstream::binary);
        if(checkpoint) {
          op.ReadBinary(f);
          checkpoint = false;
          istart = nomega;
          n_op_written += 1;
        }
      }
      nomega += 1;
    }
    if(istart>-1) std::cout << "check point: " << intfile << "_" << opname << "_" << std::to_string(istart) << " was found" << std::endl;

    int i_read = 0;
    for(int i_read=istart+1; i_read<nomega; i_read++)
    {
      std::string fn_omega = intfile + "_Omega_" + std::to_string(i_read);
      std::string fn_dump = intfile + "_" + opname + "_" + std::to_string(i_read);
      std::ifstream f(fn_omega, std::fstream::binary);
      std::ofstream fop(fn_dump, std::fstream::binary);
      imsrgsolver.Omega.back().ReadBinary(f);
      std::cout << "Transforming with " << intfile << "_Omega_" << std::to_string(i_read) << std::endl;
      op = imsrgsolver.Transform(op);
      std::cout << "Saving " << intfile + "_" + opname + "_" + std::to_string(i_read) << std::endl;
      op.WriteBinary(fop);
      n_op_written += 1;
      if(n_op_written > 2){
        std::string fn_old = intfile + "_" + opname + "_" + std::to_string(i_read-2);
        if(std::filesystem::exists(fn_old)) std::filesystem::remove(fn_old);
        n_op_written -= 1;
      }
    }

    bool renormal_order = false;
    renormal_order = modelspace.holes.size() != modelspace.core.size();
    if (renormal_order)
    {
      op = op.UndoNormalOrdering();
      op = op.DoNormalOrderingCore();
    }

    if ( ((op.GetJRank()+op.GetTRank()+op.GetParity())<1) and (op.GetNumberLegs()%2==0) )
    {
      std::cout << "writing scalar files " << std::endl;
      if (valence_file_format == "tokyo")
      {
        rw.WriteTokyo(op,intfile+"_"+opname+".snt", "op");
      }
      else
      {
        rw.WriteNuShellX_op(op,intfile+opname+".int");
      }
    }
    else
    {
      std::cout << "writing tensor files " << std::endl;
      if (valence_file_format == "tokyo")
      {
        if (op.GetJRank()==0 and (op.GetTRank()!=0 or op.GetParity()!=0) )
        {
          op.MakeReduced();
        }

        rw.WriteTensorTokyo(intfile+"_"+opname+".snt",op);
      }
      else
      {
        rw.WriteTensorOneBody(intfile+opname+"_1b.op",op,opname);
        rw.WriteTensorTwoBody(intfile+opname+"_2b.op",op,opname);
      }
    }

  }// for opnames
}


