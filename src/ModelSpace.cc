#include "ModelSpace.hh"
#include "AngMom.hh"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

Orbit::Orbit()
{
}

Orbit::Orbit(int nn, int ll, int jj2, int ttz2, int hhvq, float e=0.0)
{
   n = nn;
   l = ll;
   j2 = jj2;
   tz2 = ttz2;
   hvq = hhvq;
   spe = e;
}

Orbit::Orbit(const Orbit& orb)
{
  n = orb.n;
  l = orb.l;
  j2 = orb.j2;
  tz2 = orb.tz2;
  hvq = orb.hvq;
  spe = orb.spe;
}

void Orbit::Set(int nn, int ll, int jj2, int ttz2, int hhvq, float e)
{
   n = nn;
   l = ll;
   j2 = jj2;
   this->tz2 = ttz2;
   hvq = hhvq;
   spe = e;
}


//************************************************************************

Ket::Ket(ModelSpace * modelspace, int pp, int qq)
{
   ms = modelspace;
   p = pp;
   q = qq;
   Orbit * op = ms->GetOrbit(p);
   Orbit * oq = ms->GetOrbit(q);
   parity = (op->l + oq->l)%2;
   Tz = (op->tz2 + oq->tz2)/2;
   Jmin = abs(op->j2 - oq->j2)/2;
   Jmax = (op->j2 + oq->j2)/2;
   Jstep = 1;
   if (p==q) // Pauli principle
   { 
      Jmax--;
      Jstep++;
   }
}

int Ket::Phase(int J)
{
   int exponent = (ms->GetOrbit(p)->j2 + ms->GetOrbit(q)->j2)/2 + J + 1;
   return 1- 2*(exponent%2);
}


//************************************************************************

TwoBodyChannel::TwoBodyChannel()
{}

TwoBodyChannel::TwoBodyChannel(int j, int p, int t, ModelSpace *ms)
{
  Initialize(JMAX*2*t + JMAX*p + j, ms);
}

TwoBodyChannel::TwoBodyChannel(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}

void TwoBodyChannel::Initialize(int N, ModelSpace *ms)
{
   J = N%JMAX;
   parity = (N/JMAX)%2;
   Tz = (N/(2*JMAX)-1);
   modelspace = ms;
   NumberKets = 0;
   int nk = modelspace->GetNumberKets();
   KetMap.resize(nk,-1); // set all values to -1
   for (int i=0;i<nk;i++)
   {
      Ket *ket = modelspace->GetKet(i);
      if ( CheckChannel_ket(ket) )
      {
         KetMap[i] = NumberKets;
         KetList.push_back(i);
         int hp = ms->GetOrbit(ket->p)->hvq;
         int hq = ms->GetOrbit(ket->q)->hvq;
         if (hp==0 and hq==0) KetIndex_hh.push_back(NumberKets);
         else if (hp>0 and hq>0) KetIndex_pp.push_back(NumberKets);
         else KetIndex_ph.push_back(NumberKets);
         NumberKets++;
      }
   }
   // Set up projectors which are used in the commutators
   Proj_pp = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   Proj_hh = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   //Proj_ph_cc = arma::mat(NumberKets, NumberKets, arma::fill::zeros);
   Proj_ph_cc = arma::mat(2*NumberKets, 2*NumberKets, arma::fill::zeros);
   for (int i=0;i<NumberKets;i++)
   {
      Ket *ket = GetKet(i);
      int hvqa = modelspace->GetOrbit(ket->p)->hvq;
      int hvqb = modelspace->GetOrbit(ket->q)->hvq;
      int j2a = modelspace->GetOrbit(ket->p)->j2;
      int j2b = modelspace->GetOrbit(ket->q)->j2;
      if ( hvqa ==0 and hvqb==0 )
      {
         Proj_hh(i,i) = 1;
      }
      if ( hvqa>0 and hvqb>0 )
      {
         Proj_pp(i,i) = 1;
      }
      if ( hvqa==0 and hvqb>0 )
      {
//         Proj_ph_cc(i,i) = 1-2*((int(j2a+j2b)/2)%2);
//         Proj_ph_cc(i+NumberKets,i+NumberKets) = -(1-2*((int(j2a+j2b)/2)%2));
         Proj_ph_cc(i,i) = 1;
         Proj_ph_cc(i+NumberKets,i+NumberKets) = -1;
      }
      if ( hvqa>0 and hvqb==0 )
      {
         //Proj_ph_cc(i,i) = 1;
         Proj_ph_cc(i,i) = -(1);
         //Proj_ph_cc(i,i) = -(1-2*((int(j2a+j2b)/2)%2));
      }
   }

}


void TwoBodyChannel::Copy( const TwoBodyChannel& rhs)
{
   J                 = rhs.J;
   parity            = rhs.parity;
   Tz                = rhs.Tz;
   modelspace        = rhs.modelspace;
   NumberKets        = rhs.NumberKets;
   Proj_hh           = rhs.Proj_hh;
   Proj_pp           = rhs.Proj_pp;
   Proj_ph_cc        = rhs.Proj_ph_cc;
   KetMap            = rhs.KetMap;
   KetList           = rhs.KetList;
   KetIndex_pp       = rhs.KetIndex_pp;
   KetIndex_ph       = rhs.KetIndex_ph;
   KetIndex_hh       = rhs.KetIndex_hh;
}


int TwoBodyChannel::GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];}; 
Ket * TwoBodyChannel::GetKet(int i) const { return modelspace->GetKet(KetList[i]);}; // get pointer to ket using local index


bool TwoBodyChannel::CheckChannel_ket(int p, int q) const
{
   if ((p==q) and (J%2 != 0)) return false; // Pauli principle
   Orbit * op = modelspace->GetOrbit(p);
   Orbit * oq = modelspace->GetOrbit(q);
   if ((op->l + oq->l)%2 != parity) return false;
   if ((op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}

//************************************************************************

TwoBodyChannel_CC::TwoBodyChannel_CC()
{}

TwoBodyChannel_CC::TwoBodyChannel_CC(int j, int p, int t, ModelSpace *ms)
{
  Initialize(JMAX*2*t + JMAX*p + j, ms);
}

TwoBodyChannel_CC::TwoBodyChannel_CC(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}



bool TwoBodyChannel_CC::CheckChannel_ket(int p, int q) const
{
//   if ((p==q) and (J%2 != 0)) return false; // Pauli principle
   Orbit * op = modelspace->GetOrbit(p);
   Orbit * oq = modelspace->GetOrbit(q);
   if ((op->l + oq->l)%2 != parity) return false;
   if (abs(op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}


//************************************************************************

ModelSpace::ModelSpace()
{
//   nCore = 0;
   norbits = 0;
   maxj = 0;
   hbar_omega=20;
   target_mass = 16;
}


ModelSpace::ModelSpace(const ModelSpace& ms)
{
//   nCore = 0;
   norbits = 0;
   maxj = 0;
   hbar_omega = ms.hbar_omega;
   target_mass = ms.target_mass;
   int norbits = ms.GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      AddOrbit( Orbit(*ms.GetOrbit(i)) );
   }
   SetupKets();
//   hole = hole;
//   particle=particle;
//   valence = valence;
//   qspace = qspace;
}


ModelSpace ModelSpace::operator=(const ModelSpace& ms)
{
   return ModelSpace(ms);
}

void ModelSpace::AddOrbit(Orbit orb)
{
   int ind = Index1(orb.n, orb.l, orb.j2, orb.tz2);
   if (Orbits.size() <= ind) Orbits.resize(ind+1,NULL);
   Orbits[ind] = (new Orbit(orb));
   norbits = Orbits.size();
   if (orb.j2 > maxj)
   {
      maxj = orb.j2;
      //nTwoBodyChannels = 2*3*(maxj+1);
      nTwoBodyChannels = 2*3*(JMAX);
   }
//   if (orb.hvq == 0) nCore+=orb.j2+1; // 0 means hole (ie core), 1 means valence, 2 means outside the model space
   if (orb.hvq == 0) holes.push_back(ind);
   if (orb.hvq > 0) particles.push_back(ind);
   if (orb.hvq == 1) valence.push_back(ind);
   if (orb.hvq == 2) qspace.push_back(ind);
}

int ModelSpace::GetTwoBodyChannelIndex(int j, int p, int t)
{
   return (t+1)*2*JMAX + p*JMAX + j;
}


void ModelSpace::SetupKets()
{
   int index = 0;

   for (int p=0;p<norbits;p++)
   {
     for (int q=p;q<norbits;q++)
     {
        //index = q*(q+1)/2 + p;
        //index = q*norbits + p;
        index = Index2(p,q);
        if (index >= Kets.size()) Kets.resize(index+1,NULL);
        Kets[index] = new Ket(this,p,q);
     }
   }
   for (int ch=0;ch<nTwoBodyChannels;++ch)
   {
      TwoBodyChannels.push_back(TwoBodyChannel(ch,this));
      TwoBodyChannels_CC.push_back(TwoBodyChannel_CC(ch,this));
      cout << "ch = " << ch << " " << TwoBodyChannels[ch].J << " " << TwoBodyChannels[ch].parity << " " << TwoBodyChannels[ch].Tz << " nKets = " << TwoBodyChannels[ch].GetNumberKets()
      << "  " << TwoBodyChannels_CC[ch].GetNumberKets() << endl;
   }
}




//double ModelSpace::GetSixJ(int j1, int j2, int j3, int J1, int J2, int J3)
double ModelSpace::GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3)
{
// Don't really need to store all of them, only need to store
// unique combinations. Implement this if it becomes a speed/storage issue
//   int jlist[6] = {j1,j2,j3,J1,J2,J3};
//   int klist[6];
//   kmin = =std::min_element(jlist,jlist+6);
   // use 2J in the key so we don't have to worry about half-integers
   int k1 = int(2*j1);
   int k2 = int(2*j2);
   int k3 = int(2*j3);
   int K1 = int(2*J1);
   int K2 = int(2*J2);
   int K3 = int(2*J3);
   long int key = 10000000000*k1 + 100000000*k2 + 1000000*k3 + 10000*K1 + 100*K2 + K3;
//   long int key = 10000000000*j1 + 100000000*j2 + 1000000*j3 + 10000*J1 + 100*J2 + J3;
   map<long int,double>::iterator it = SixJList.find(key);
   if ( it != SixJList.end() )  return it->second;

   double sixj = AngMom::SixJ(j1,j2,j3,J1,J2,J3);
   SixJList[key] = sixj;
   return sixj;

}


double ModelSpace::GetNineJ(double j1,double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
{
   double ninej = 0;
   //int ph = 1-2*(abs(int(J-j1))%2);
   int ph = 1-2*(int(2*j1+J)%2);
//   cout << "NineJ:   { " << j1  << " " << j2  << " " << J12 << " }" << endl
//        << "         { " << j3  << " " << j4  << " " << J34 << " }" << endl
//        << "         { " << J13 << " " << J24 << " " << J   << " }" << endl
//        << endl;
//   cout << "ph = " << ph << endl;
   for (float g = fabs(J-j1); g<=J+j1; g+=1)
   {
      ninej +=  ph * (2*g+1)
                * GetSixJ(j1,j2,J12,J34,J,g)
                * GetSixJ(j3,j4,J34,j2,g,J24)
                * GetSixJ(J13,J24,J,g,j1,j3);
//                cout << "g = " << g << "  Six-J's: "
//                     << GetSixJ(j1,j2,J12,J34,J,g) << " " 
//                     << GetSixJ(j3,j4,J34,j2,g,J24) << " " 
//                     << GetSixJ(J13,J24,J,g,j1,j3) << endl;
   }
   return ninej;
}
