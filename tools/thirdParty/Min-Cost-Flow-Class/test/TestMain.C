/*--------------------------------------------------------------------------*/
/*-------------------------- File TestMain.C -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing MCF solvers derived from MCFClass. According to the
 * command line parameters, a "trusted" solver and a solver still being
 * tested are constructed. An instance of a MCF in DIMACS format is read from
 * file. The problem is then repeatedly solved with several changes in
 * costs/capacities/deficits and arcs openings/closures. The same operations
 * are performed on the two solvers, and the results are printed out. If the
 * results don't match, then at least one of the two solvers is incorrect. If
 * the results match, chances are the two solvers are correct.
 *
 * \version 2.01
 *
 * \date 24 - 04 - 2017
 *
 * \author Alessandro Bertolini \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 1996 - 2017 by Antonio Frangioni, Claudio Gentile
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINES -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define HAVE_CSCL2 1
// > 0 if the CS2 class is available

#define HAVE_CPLEX 1
// > 0 if the MCFCplex class is available

#define HAVE_MFSMX 1
// > 0 if the MCFSimplex class is available

#define HAVE_MFZIB 1
// > 0 if the MCFZIP class is available

#define HAVE_RELAX 1
// > 0 if the RelaxIV class is available

#define HAVE_SPTRE 1
// > 0 if the SPTree class is available
// NOTE: SPTree cannot solve most MCF instances, so this may result in errors
//       even if SPTree and the other solvers are "correct"

#define NMS_IS_USED 0

// if NMS_IS_USED > 0, then the Chg****() routines are fed with a
// non-consecutive set of names; otherwise, all the involved arcs are
// consecutive

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <iomanip>

#if( HAVE_CSCL2 )
 #include "CS2.h"
#endif

#if( HAVE_CPLEX )
 #include "MCFCplex.h"
#endif

#if( HAVE_MFSMX )
 #include "MCFSimplex.h"
#endif

#if( HAVE_MFZIB )
 #include "MCFZIB.h"
#endif

#if( HAVE_RELAX )
 #include "RelaxIV.h"
#endif

#if( HAVE_SPTRE )
 #include "SPTree.h"
#endif


/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace MCFClass_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/

/*!!
template<class T>
static inline T min( const T x , const T y )
{
 return( x <= y ? x : y );
 }
!!*/

/*--------------------------------------------------------------------------*/

/*!!
template<class T>
static inline T max( const T x , const T y )
{
 return( x >= y ? x : y );
 }
!!*/

/*--------------------------------------------------------------------------*/

template<class T>
static inline void Str2Sthg( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/

static inline MCFClass *CreateProb( int wSlvr , int Optns )
{
 MCFClass *mcf = NULL;  // unknown solver, or the required solver is not
                        // available due to the macroes settings
 switch( wSlvr ) {
  #if( HAVE_RELAX )  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( 0 ): {
    RelaxIV *rlx = new RelaxIV();
    #if( AUCTION )
     if( Optns )
      rlx->SetPar( RelaxIV::kAuction , MCFClass::kYes );
    #endif
    mcf = rlx;
    cout << "RelaxIV";
    break;
    }
  #endif
  #if( HAVE_SPTRE )  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( 1 ): mcf = new SPTree();
              cout << "SPTree";
	      break;
  #endif
  #if( HAVE_CPLEX )  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( 2 ): {
    MCFCplex *cpx = new MCFCplex();
    if( Optns >= 0 )
     cpx->SetPar( CPX_PARAM_NETPPRIIND , Optns );
    mcf = cpx;
    cout << "MCFCplex";
    break;
    }
  #endif
  #if( HAVE_MFZIB )  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( 3 ): {
    MCFZIB *zib = new MCFZIB();
    bool PrmlSmplx = Optns & 1;
    char Prcng;
    switch( Optns / 2 ) {
     case( 0 ): Prcng = char( MCFZIB::kDantzig ); break;
     case( 1 ): Prcng = char( MCFZIB::kFrstElA ); break;
     default:   Prcng = char( MCFZIB::kMltPrPr );
     }
    if( ( ! PrmlSmplx ) && ( Prcng == MCFZIB::kDantzig ) )
     Prcng = char( MCFZIB::kMltPrPr );
    zib->SetAlg( PrmlSmplx , Prcng );
    mcf = zib;
    cout << "MCFZIB";
    break;
    }
  #endif
  #if( HAVE_CSCL2 )  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( 4 ): mcf = new CS2();
              cout << "CS2";
	      break;
  #endif
  #if( HAVE_MFSMX )  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( 5 ): {
    MCFSimplex *mcfs = new MCFSimplex();
    bool PrmlSmplx = Optns & 1;
    char Prcng;
    switch( Optns / 2 ) {
     case( 0 ): Prcng = char( MCFSimplex::kDantzig ); break;
     case( 1 ): Prcng = char( MCFSimplex::kFirstEligibleArc ); break;
     default:   Prcng = char( MCFSimplex::kCandidateListPivot );
     }
    if( ( ! PrmlSmplx ) && ( Prcng == MCFSimplex::kDantzig ) )
     Prcng = char( MCFSimplex::kCandidateListPivot );
    mcf = mcfs;
    cout << "MCFSimplex";
    break;
    }
  #endif
  default:    mcf = NULL;
              cout << "unknown";

  }  // end( switch( wSlvr ) )- - - - - - - - - - - - - - - - - - - - - - - -

 return( mcf );
 }

/*--------------------------------------------------------------------------*/

static inline void PrintResult( MCFClass *mcf )
{
 switch( mcf->MCFGetStatus() ) {
  case( MCFClass::kOK ):         cout << mcf->MCFGetFO();
                                 break;
  case( MCFClass::kUnfeasible ): cout << "        +INF";
                                 break;
  case( MCFClass::kUnbounded ):  cout << "        -INF";
                                 break;
  default:                       cout << "      Error!";
  }
 }

/*--------------------------------------------------------------------------*/

static inline void SolveMCF( MCFClass *mcf1 , MCFClass *mcf2 ) 
{
 try {
  mcf1->SolveMCF();
  cout << "MCF1 = ";
  PrintResult( mcf1 );
  mcf2->SolveMCF();
  cout << ", MCF2 = ";
  PrintResult( mcf2 );
  cout << endl;
  }
 catch( exception &e ) {
  cerr << "MCF1: " << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int seed = 1;
 MCFClass::Index n_change = 10;
 MCFClass::Index n_repeat = 10;
 int reoptmz = 1;
 int optns1 = 0;
 int optns2 = 0;
 int wmcf1;
 int wmcf2;

 ifstream iFile;

 switch( argc ) {
  case( 10 ): Str2Sthg( argv[ 9 ] , seed );
  case( 9 ):  Str2Sthg( argv[ 8 ] , n_change );
  case( 8 ):  Str2Sthg( argv[ 7 ] , n_repeat );
  case( 7 ):  Str2Sthg( argv[ 6 ] , reoptmz );
  case( 6 ):  Str2Sthg( argv[ 5 ] , optns2 );
  case( 5 ):  Str2Sthg( argv[ 4 ] , optns1 );
  case( 4 ):  iFile.open( argv[ 1 ] );
              if( ! iFile ) {
	       cerr << "ERROR: opening input file " << argv[ 1 ] << endl;
	       return( 1 );
	       }
	      Str2Sthg( argv[ 2 ] , wmcf1 );
	      Str2Sthg( argv[ 3 ] , wmcf2 );
              break;
  default: cerr <<
	   "Usage: MCFSolve <input file> <MCF1> <MCF2> [<optns1>] [<optns2>]"
		<< endl <<
           "       [<reoptmz>][<# repeats>][<# changes>][<seed>]"
		<< endl <<
	   " MCFx: 0 = Relax, 1 = SPTree, 2 = Cplex, 3 = MCFZIB, 4 = CS2"
                << ", 5 = MCFSimplex" << endl <<
	   " optnsx: Relax   : > 0 uses Auction"
		<< endl <<
	   "         Cplex   : network pricing parameter"
		<< endl <<
	   "         ZIB     : 1st bit == 1 ==> primal +"
		<< endl <<
	   "                   0 = Dantzig, 2 = First Eligible, 4 = MPP"
	        << endl <<
	   "         Simplex : 1st bit == 1 ==> primal +"
		<< endl <<
	   "                   0 = Dantzig, 2 = First Eligible, 4 = MPP"
	        << endl;
	   return( 1 );
  }

 // construction of the objects - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cout << argv[ 1 ];

 // construct mcf1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cout << ", MCF1 = ";

 MCFClass *mcf1 = CreateProb( wmcf1 , optns1 );
 if( ! mcf1 ) {
  cerr << "Solver type " << wmcf1 << " unavailable" << endl;
  return( 1 );
  }

 if( ! reoptmz )
  mcf1->SetPar( MCFClass::kReopt , MCFClass::kNo );

 mcf1->SetMCFTime();  // do timing

 // construct mcf2- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cout << ", MCF2 = ";

 MCFClass *mcf2 = CreateProb( wmcf2 , optns2 );
 if( ! mcf2 ) {
  cerr << "Solver type " << wmcf2 << " unavailable" << endl;
  return( 1 );
  }

 if( ! reoptmz )
  mcf2->SetPar( MCFClass::kReopt , MCFClass::kNo );

 mcf2->SetMCFTime();  // do timing

 // load the network- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 try {
  mcf1->LoadDMX( iFile );
  }
 catch( exception &e ) {
  cerr << "MCF1: " << e.what() << endl;
  return( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  return( 1 );
  }

 iFile.close();
 iFile.open( argv[ 1 ] );

 try {
  mcf2->LoadDMX( iFile );
  }
 catch( exception &e ) {
  cerr << "MCF2: " << e.what() << endl;
  return( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  return( 1 );
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFClass::cIndex n = mcf1->MCFn();
 MCFClass::cIndex m = mcf1->MCFm();

 cout << ", n = " << n << ", m = " << m << endl;
 if( n_change > m )
  n_change = m;

 // compute min/max cost & max deficit- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFClass::CNumber c_max = Inf<MCFClass::CNumber>();  // max cost
 MCFClass::CNumber c_min = - c_max;                   // min cost

 for( MCFClass::Index i = 0 ; i < m ; i++ )
  #if( DYNMC_MCF )
   if( ! mcf1->IsClosedArc( i ) )
  #endif
   {
    MCFClass::cCNumber ci = mcf1->MCFCost( i );
    if( ci < c_min )
     c_min = ci;

    if( ci > c_max )
     c_max = ci;
    }

 bool nzdfct = false;

 for( MCFClass::Index i = 0 ; i < n ; )
  if( mcf1->MCFDfct( i++ ) > 0 ) {
   nzdfct = true;
   break;
   }

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cout << "First call:\t\t ";
 cout.setf( ios::scientific, ios::floatfield );
 cout << setprecision( 6 );

 SolveMCF( mcf1 , mcf2 );
 
 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now, for n_repeat times:
 // - n_change costs are changed, then the two problems are re-solved;
 // - n_change capacities are changed, then the two problems are re-solved,
 //   then the original capacities are restored;
 // - if the problem is not a circulation problem, 2 deficits are modified
 //   (adding and subtracting the same number), then the two problems are
 //   re-solved, then the original deficits are restored;
 // - if the relative methods are implemented:
 //   = n_change arcs are closed, then the two problems are re-solved;
 //   = the same arcs arcs are re-opened, then the two problems are re-solved

 MCFClass::Index_Set nms = new MCFClass::Index[ n_change + 1 ];
 #if( ! NMS_IS_USED )
  MCFClass::Index strt, stp;
 #endif

 MCFClass::CRow newcsts = new MCFClass::CNumber[ n_change ];
 MCFClass::FRow newcaps = new MCFClass::FNumber[ max( n_change , n ) ];
 MCFClass::FRow oldcaps = new MCFClass::FNumber[ max( n_change , n ) ];

 while( n_repeat-- ) {
  // compute new costs- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  MCFClass::Index i = 0;
  for( ; i < n_change ; i++ ) {
   newcsts[ i ] = c_min + MCFClass::CNumber( drand48() * ( c_max - c_min ) );

   #if( NMS_IS_USED )
    if( i )
     nms[ i ] = nms[ i - 1 ] +
                max( MCFClass::Index( 1 ) ,
		     MCFClass::Index( drand48() * ( m - nms[ i - 1 ] ) )
		     );
    else
     nms[ i ] = MCFClass::Index( drand48() * ( m / 2 ) );

    if( nms[ i ] == m - 1 )
     break;
   #endif
   }

  #if( NMS_IS_USED )
   nms[ i ] = Inf<MCFClass::Index>();
  #else
   strt = MCFClass::Index( drand48() * ( m - n_change ) );
   stp = strt + n_change;
  #endif

  // change costs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( NMS_IS_USED )
   mcf1->ChgCosts( newcsts , nms );
   mcf2->ChgCosts( newcsts , nms );
  #else
   mcf1->ChgCosts( newcsts , NULL , strt , stp );
   mcf2->ChgCosts( newcsts , NULL , strt , stp );
  #endif

  // re-solve the problems - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "Changing costs:\t\t ";
  SolveMCF( mcf1 , mcf2 );

  // compute new capacities - - - - - - - - - - - - - - - - - - - - - - - - -

  for( i = 0 ; i < n_change ; i++ ) {
   newcaps[ i ] = 1;

   /*!! This would appear to be more fair, but it very seldom changes
        anything in the solution:

   newcaps[ i ] = cap_min + FNumber( drand48() * ( cap_max - cap_min ) );
   !!*/

   #if( NMS_IS_USED )
    if( i )
     nms[ i ] = nms[ i - 1 ] +
                max( MCFClass::Index( 1 ) ,
		     MCFClass::Index( drand48() * ( m - nms[ i - 1 ] ) )
		     );
    else
     nms[ i ] = MCFClass::Index( drand48() * ( m / 2 ) );

    if( nms[ i ] == m - 1 )
     break;
   #endif
   }

  #if( NMS_IS_USED )
   nms[ i ] = Inf<MCFClass::Index>();
  #else
   strt = MCFClass::Index( drand48() * ( m - n_change ) );
   stp = strt + n_change;
  #endif

  // save the old capacities- - - - - - - - - - - - - - - - - - - - - - - - -

  #if( NMS_IS_USED )
   mcf1->MCFUCaps( oldcaps , nms );
  #else
   mcf1->MCFUCaps( oldcaps , NULL , strt , stp );
  #endif

  // change capacities- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( NMS_IS_USED )
   mcf1->ChgUCaps( newcaps , nms );
   mcf2->ChgUCaps( newcaps , nms );
  #else
   mcf1->ChgUCaps( newcaps , NULL , strt , stp );
   mcf2->ChgUCaps( newcaps , NULL , strt , stp );
  #endif

  // re-solve the problems- - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "Changing capacities:\t ";
  SolveMCF( mcf1 , mcf2 );

  // revert back to the old capacities- - - - - - - - - - - - - - - - - - - -

  #if( NMS_IS_USED )
   mcf1->ChgUCaps( oldcaps , nms );
   mcf2->ChgUCaps( oldcaps , nms );
  #else
   mcf1->ChgUCaps( oldcaps , NULL , strt , stp );
   mcf2->ChgUCaps( oldcaps , NULL , strt , stp );
  #endif

  // compute new deficits - - - - - - - - - - - - - - - - - - - - - - - - - -
  // since n >= 2, oldcaps[] and newcaps[] have at least 2 elements

  if( nzdfct ) {  // ... if there are nonzero deficits at all
   mcf1->MCFDfcts( oldcaps );

   do
    i = MCFClass::Index( drand48() * n );  // select one node with positive
   while( oldcaps[ i ] <= 0 );             // deficit (one must exist)

   newcaps[ 0 ] = oldcaps[ i ];
   nms[ 0 ] = i;

   do
    i = MCFClass::Index( drand48() * n );  // select one node with negative
   while( oldcaps[ i ] >= 0 );             // deficit (one must exist)

   newcaps[ 1 ] = oldcaps[ i ];
   nms[ 1 ] = i;
   nms[ 2 ] = Inf<MCFClass::Index>();

   // remember the old values of the deficit in oldcaps[]
   oldcaps[ 0 ] = newcaps[ 0 ];
   oldcaps[ 1 ] = newcaps[ 1 ];

   MCFClass::FNumber Dlt = min( newcaps[ 0 ] , - newcaps[ 1 ] );
   Dlt = max( MCFClass::FNumber( 1 ) ,
	      MCFClass::FNumber( drand48() * Dlt ) );
   newcaps[ 0 ] -= Dlt;
   newcaps[ 1 ] += Dlt;

   // change deficits - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   mcf1->ChgDfcts( newcaps , nms );
   mcf2->ChgDfcts( newcaps , nms );

   // re-solve the problems - - - - - - - - - - - - - - - - - - - - - - - - -

   cout << "Changing deficits:\t ";
   SolveMCF( mcf1 , mcf2 );

   // revert back to the old deficits - - - - - - - - - - - - - - - - - - - -

   mcf1->ChgDfcts( oldcaps , nms );
   mcf2->ChgDfcts( oldcaps , nms );

   }  // end( if( dfct_max > 0 ) )

  #if( DYNMC_MCF >= 2 )
   // closing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( i = 0 ; i < n_change ; i++ ) {
    if( i )
     nms[ i ] = nms[ i - 1 ] +
      max( MCFClass::Index( 1 ) ,
	   MCFClass::Index( drand48() * ( m - nms[ i - 1 ] ) )
	   );
    else
     nms[ i ] = MCFClass::Index( drand48() * ( m / 2 ) );

    mcf1->CloseArc( nms[ i ] );
    mcf2->CloseArc( nms[ i ] );

    if( nms[ i ] == m - 1 ) {
     i++;
     break;
     }
    }

   nms[ i ] = Inf<MCFClass::Index>();

   // re-solve the problems - - - - - - - - - - - - - - - - - - - - - - - - -

   cout << "Closing arcs:\t\t ";
   SolveMCF( mcf1 , mcf2 );

   // re-opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( i = 0 ; nms[ i ] < Inf<MCFClass::Index>() ; i++ ) {
    mcf1->OpenArc( nms[ i ] );
    mcf2->OpenArc( nms[ i ] );
    }

   // re-solve the problems - - - - - - - - - - - - - - - - - - - - - - - - -

   cout << "Re-opening arcs:\t ";
   SolveMCF( mcf1 , mcf2 );
  #endif

  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double tu , ts;
 mcf1->TimeMCF( tu , ts );
 cout << "Time: MCF1 = " << tu + ts << ", MCF2 = ";
 mcf2->TimeMCF( tu , ts );
 cout << tu + ts << endl;
 
 // destroy objects and vectors - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] oldcaps;
 delete[] newcaps;
 delete[] newcsts;
 delete[] nms;

 delete mcf2;
 delete mcf1;

 // terminate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*----------------------- End File TestMain.C ------------------------------*/
/*--------------------------------------------------------------------------*/
