
#ifndef ReadWrite_h
#define ReadWrite_h 1
#define BOOST_IOSTREAMS_NO_LIB 1

#include <map>
#include <string>
#include "Operator.hh"

using namespace std;


class ReadWrite
{
 public:
   ~ReadWrite();
   ReadWrite();
   void ReadSettingsFile(  string filename);
   void ReadTBME_Oslo( string filename, Operator& Hbare);
//   void ReadBareTBME( string filename, Operator& Hbare);
   void ReadBareTBME_Jason( string filename, Operator& Hbare);
   void ReadBareTBME_Navratil( string filename, Operator& Hbare);
   void ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int E1max, int E2max, int lmax);
   template<class T> void ReadBareTBME_Darmstadt_from_stream( T & infile, Operator& Hbare, int E1max, int E2max, int lmax);
//   void ReadBareTBME_Darmstadt_from_stream( istream & infile, Operator& Hbare, int E1max, int E2max, int lmax);
   void Read_Darmstadt_3body( string filename, Operator& Hbare, int E1max, int E2max, int E3max);
   template<class T>void Read_Darmstadt_3body_from_stream( T & infile, Operator& Hbare, int E1max, int E2max, int E3max);
//   void Read_Darmstadt_3body_from_stream( istream & infile, Operator& Hbare, int E1max, int E2max, int E3max);
   void GetHDF5Basis( ModelSpace* modelspace, string filename, vector<array<int,5>>& Basis );
   void Read3bodyHDF5( string filename, Operator& op);
   void WriteNuShellX_sps( Operator& op, string filename);
   void WriteNuShellX_int( Operator& op, string filename);
   void WriteAntoine_int( Operator& op, string filename); // <- not implemented yet...
   void WriteOperator(Operator& op, string filename);
   void ReadOperator(Operator& op, string filename); 

   std::map<string,string> InputParameters;

   bool InGoodState(){return goodstate;};
   bool doCoM_corr;
   void SetCoMCorr(bool b){doCoM_corr = b;cout <<"Setting com_corr to "<< b << endl;};
   bool goodstate;

};



/// Wrapper class so we can treat a vector of floats like a stream, using the extraction operator >>.
/// This is used for the binary version of ReadWrite::Read_Darmstadt_3body_from_stream().
class VectorStream 
{
 public:
  VectorStream(vector<float>& v) : vec(v), i(0) {};
  VectorStream& operator>>(float& x) { x = vec[i++]; return (VectorStream&)(*this);}
  bool good(){ return i<vec.size(); };
  void getline(char[], int) {}; // Don't do nuthin'.
 private:
  vector<float>& vec;
  long long unsigned int i;
};

#endif

