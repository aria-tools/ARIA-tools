/*--------------------------------------------------------------------------*/
/*------------------------- File RelaxIV.C ---------------------------------*/
/*--------------------------------------------------------------------------*/
/* Linear Min Cost Flow problems solver, based on the RELAXIV code by
 * D. Bertsekas and P. Tseng, as described in
 *
 *    Bertsekas, Dimitri P., and Paul Tseng.
 *    "RELAX-IV: A faster version of the RELAX code for solving minimum
 *     cost flow problems." (1994), Report LIDS-P-2276, MIT.
 *
 * Conforms to the standard (MCF) interface defined in MCFClass.h.
 *
 * \version 1.83
 *
 * \date 14 - 04 - 2019
 *
 * \author <b>(original FORTRAN code)</b> \n
 *         Dimitri P. Bertsekas \n
 *         Lab. for Information and Decision Systems \n
 *         Massachusetts Institute of Technology \n
 *
 * \author <b>(original FORTRAN code)</b> \n
 *         Paul Tseng \n
 *         Department of Mathematics \n
 *         University of Washington \m
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * Copyright &copy 1996 - 2019 by Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "RelaxIV.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace MCFClass_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static MCFClass::cIndex npasslim = 2;

/* Number of single-node iterations that are attempted on all nodes before
   the first multinode iteration is allowed. */

static MCFClass::cIndex tp     = 10;
static MCFClass::cIndex ts_den = 15;

/* An adaptive strategy is used to decide whether to continue the scanning
   process after a multinode price change. The thresold parameter tp tells
   what is a "small number of nonzero deficit nodes". The thresold parameter
   ts, defined as ts = n / ts_den, is used to decide whether or not to
   continue scanning even after a price change. */

static int it_aug = 8; 

/* If the number of iteration is greater of it_aug * the number of AugFlow()
   calls (and/or if other conditions are verified), then another multinode
   iteration is performed. */

static MCFClass::cIndex maxdns = 10;

/* In the standard initialization, the number of passes is controlled by the
   density of the graph: if m / n > maxdns then 2 passes are done, otherwise
   3 passes are done. */

#if( AUCTION )
 static const int factor    = 3;
 static const int npassauct = 1;
 static const int maxdf     = 8;

 /* Auction parameters:
    - factor determines by how much eps is reduced at each minimization;
    - npassauct determines how many auction scaling iterations are
      performed, that is how many times eps is divided by factor;
    - maxdf is used to set the initial value of eps in auction function:
      eps = max( 1 , ( maxcost - mincost ) / maxdf ), where maxcost and
      mincost are, respectively, the maximum and the mininum reduced cost
      at the beginning of auction function. */

 static cCNumber C_LARGE = Inf<CNumber>() / 4;
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- "PRIVATE" MACROS ------------------------------*/
/*--------------------------------------------------------------------------*/

#define P_ALLOC ( AUCTION || ( DYNMC_MCF_RIV > 1 ) )

#if( RELAXIV_STATISTICS )
 #define pp( x ) x++
#else
 #define pp( x )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF RIVState -------------------------*/
/*--------------------------------------------------------------------------*/

RelaxIV::RIVState::RIVState( MCFClass::cIndex m )
{
 Flow    = new RelaxIV::FNumber[ m ];
 RedCost = new CNumber[ m ];
 }

RelaxIV::RIVState::~RIVState()
{
 delete[] Flow;
 delete[] RedCost;
 }

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF RelaxIV --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

RelaxIV::RelaxIV( cIndex nmx , cIndex mmx )
         :
         MCFClass( nmx , mmx )
{
 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( nmax && mmax )
  MemAlloc();
 else
  nmax = mmax = 0;

 // other initializations - - - - - - - - - - - - - - - - - - - - - - - - - -

 InstCntr++;

 }  // end( RelaxIV )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void RelaxIV::LoadNet( cIndex nmx , cIndex mmx , cIndex pn , cIndex pm ,
                       cFRow pU , cCRow pC , cFRow pDfct , cIndex_Set pSn ,
                       cIndex_Set pEn )
{
 // allocating and deallocating memory- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( nmx != nmax ) || ( mmx != mmax ) ) {
  if( nmax && mmax ) {
   MemDeAlloc();
   nmax = mmax = 0;
   }

  if( mmx && nmx ) {
   nmax = nmx;
   mmax = mmx;
   MemAlloc();
   }
  }

 if( ( ! nmax ) || ( ! mmax ) ) {  // just sit down in the corner and wait
  nmax = mmax = 0;
  return;
  }

 // now setting up data - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 n = pn;
 m = pm;

 if( pDfct ) {  // copy the deficits- - - - - - - - - - - - - - - - - - - - -
  FRow tB = B + n;
  for( pDfct += n ; tB > B ; )
   *(tB--) = *(--pDfct);
  }
 else         // deficits are all-0
  for( FRow tB = B + n ; tB > B ; )
   *(tB--) = 0;

 bool InfCap = false;
 if( pU ) {  // copy the capacities - - - - - - - - - - - - - - - - - - - - -
  FRow tCap = Cap + m;
  for( pU += m ; tCap > Cap ; )
   if( ( (*tCap--) = *(--pU) ) == Inf<FNumber>() ) {
    InfCap = true;
    break;
    }

  for( ; tCap > Cap ; )
   (*tCap--) = *(--pU);
  }
 else {    // capacities are all-INF
  InfCap = true;
  for( FRow tCap = Cap + m ; tCap > Cap ; )
   (*tCap--) = Inf<FNumber>();
  }

 if( pC ) {  // copy the costs- - - - - - - - - - - - - - - - - - - - - - - -
  CRow tC = C + m;
  CRow tRC = RC + m;
  #if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
   FRow tCap = Cap + m;
   for( pC += m ; tC > C ; tCap-- ) {
    CNumber ttC = *(--pC);
    if( ttC == Inf<CNumber>() ) {
     *tCap = 0;
     ttC = 0;
     }

    *(tC--) = *(tRC--) = ttC;
    }
  #else
   for( pC += m ; tC > C ; ) {
    cCNumber ttC = *(--pC);
    if( ( *(tRC--) = ttC ) < Inf<CNumber>() )
     *(tC--) = ttC;
    else
     *(tC--) = 0;
    }
  #endif
  }
 else {  // costs are all-0 - - - - - - - - - - - - - - - - - - - - - - - - -
  CRow tRC = RC + m;
  for( CRow tC = C + m ; tC > C ; )
   *(tC--) = *(tRC--) = 0;
  }

 if( InfCap ) {  // make all capacities finite- - - - - - - - - - - - - - - -
  FNumber maxcap = 0;
  for( FRow tB = B + n ; tB > B ; tB-- )
   if( *tB > 0 )
    maxcap += *tB;

  FRow tCap = Cap + m;
  for( CRow tC = C + m ; tC > C ; tCap-- )
   if( *(tC--) < 0 ) {
    if( *tCap == Inf<FNumber>() )
     throw(
      MCFException( "RelaxIV::LoadNet(): C[ i , j ] < 0 and U[ i , j ] = INF"
                    ) );

    maxcap += *tCap;
    }

  for( tCap = Cap + m ; tCap > Cap ; tCap-- )
   if( *tCap == Inf<FNumber>() )
    *tCap = maxcap;
  }

 #if( SAME_GRPH_RIV )
  if( ! Startn[ 1 ] )  // Startn[] and Endn[] have not been initialized yet -
 #endif
  {
   Index_Set tEn = Endn + m;
   Index_Set tSn = Startn + m;
   for( pSn += m , pEn += m ; tSn > Startn ; ) {
    *(tSn--) = *(--pSn) + USENAME0;
    *(tEn--) = *(--pEn) + USENAME0;
    }
   }

 #if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
  if( FIn[ 1 ] == Inf<Index>() )  // FS and BS have not been initialized yet-
 #endif
  {
   // clean up the FS and BS information- - - - - - - - - - - - - - - - - - -

   Index_Set tOu = FOu + n;
   for( Index_Set tIn = FIn + n ; tIn > FIn ; )
    *(tIn--) = *(tOu--) = 0;

   // now construct the FS and BS structures- - - - - - - - - - - - - - - - -

   for( Index j = 0 ; j++ < m ; )
    #if( ( ! SAME_GRPH_RIV ) || DYNMC_MCF_RIV )
     if( RC[ j ] == Inf<CNumber>() )
      X[ j ] = 0;
     else
    #endif
     {
      Index i = Startn[ j ];

      NxtOu[ j ] = FOu[ i ];
      FOu[ i ] = j;

      NxtIn[ j ] = FIn[ i = Endn[ j ] ];
      FIn[ i ] = j;
      }
   }

 #if( P_ALLOC )
  PiOwnr = NULL;
 #endif
 #if( AUCTION )
  crash = false;
 #endif

 #if( RELAXIV_STATISTICS )
  iter = nmultinode = num_augm = num_ascnt = 0;
  #if( AUCTION )
   nsp = 0;
  #endif
 #endif

 #if( DYNMC_MCF_RIV > 2 )
  ffp = 0;
 #endif

 status = MCFClass::kUnSolved;

 }  // end( LoadNet )

/*--------------------------------------------------------------------------*/

void RelaxIV::PreProcess( void )
{
 Index_Set tFOu = FOu + n;
 Index_Set tFIn = FIn + n;
 for( FRow tB = B + n ; tB > B ; tB-- ) {
  FNumber dfctn = *tB;
  FNumber tcap = 0;
  Index arc = *tFOu;

  while( arc ) {
   tcap += Cap[ arc ];
   arc = NxtOu[ arc ];
   }

  FNumber cap = tcap + dfctn;

  if( LTZ( cap , EpsDfct ) ) {  // problem is unfeasible
   status = MCFClass::kUnfeasible;
   error_node = tB - B;
   error_info = 1;
   return;
   }

  tcap = 0;
  for( arc = *(tFIn--) ; arc ; ) {
   if( cap < Cap[ arc ] )
    Cap[ arc ] = cap;

   tcap += Cap[ arc ];
   arc = NxtIn[ arc ];
   }

  cap = tcap - dfctn;

  if( LTZ( cap , EpsDfct ) ) {  // problem is unfeasible
   status = MCFClass::kUnfeasible;
   error_node = tB - B;
   error_info = 2;
   return;
   }

  for( arc = *(tFOu--) ; arc ; ) {
   if( cap < Cap[ arc ] )
    Cap[ arc ] = cap;

   arc = NxtOu[ arc ];
   }
  }

 status = MCFClass::kUnSolved;

 }  // end( PreProcess )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

void RelaxIV::SolveMCF( void )
{
 if( MCFt )
  MCFt->Start();

 #if( P_ALLOC )
  PiOwnr = NULL;
 #endif

 FO = Inf<FONumber>();
 iter = num_augm = 0;
 #if( RELAXIV_STATISTICS )
  nmultinode = num_ascnt = 0;
 #endif

 if( status )  // initializations are skipped if status == 0- - - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  status = MCFClass::kOK;

  // prices and flows are initialized by either calling the auction()
  // routine or by performing only single-node iterations

  #if( AUCTION )
   if( crash )
    auction();
   else
  #endif 
    init_standard();

  if( status )
   return;

  init_tree();
  }

 // initialize other variables- - - - - - - - - - - - - - - - - - - - - - - -

 Bool_Vec tmark = mark + n;
 for( Bool_Vec tscan = scan + n ; tscan > scan ; )
  *(tmark--) = *(tscan--) = false;

 // an adaptive strategy is used to decide whether to continue the scanning
 // process after a multinode price change: the thresold parameters that
 // control this strategy are tp and ts, that is set in the next line

 cIndex ts = n / ts_den;

 // initialize the queue of nodes with nonzero deficit- - - - - - - - - - - -

 Index node = 2;
 for( Index_Set tnxtq = queue ; node <= n ; )
  *(++tnxtq) = node++;

 queue[ lastq = n ] = 1;
 
 FNumber deficit;
 Index ndfct = node = n;
 Index nnonz = 0;
 Index nlabel = 0;
 Index npass = 0;

 bool Switch = false;

 for(;;)  // main loop, repeated until there are unbalanced nodes - - - - - -
 {        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bool posit = false;

  for(;;) {  // advancing the queue of nonzero deficit nodes- - - - - - - - -
   prvnde = node;
   node = queue[ node ];
   deficit = Dfct[ node ];

   if( node == lastq ) {
    ndfct = nnonz;
    nnonz = 0;
    lastq = prvnde;
    npass++;
    }

   // deleting a node from the queue- - - - - - - - - - - - - - - - - - - - -

   if( ETZ( deficit , EpsDfct ) ) {
    Index nxtnode = queue[ node ];

    if( node == nxtnode ) {
     posit = true;  // condition for termination of SolveMCF
     break;
     }
    else {
     queue[ node ] = 0;
     queue[ prvnde ] = node = nxtnode;
     }
    }
   else   // selected a node for the current relaxation iteration
    break;
   }

  if( posit )
   break;  // terminate main loop

  iter++;
  nnonz++;
  bool quit;

  if( ( posit = GTZ( deficit , EpsDfct ) ) ) {
   // attempt a single node iteration from node with positive deficit- - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   bool pchange = false;
   FNumber indef = deficit;
   FNumber delx = svblncdarcs( node , tfstou , tnxtou , tfstin , tnxtin );

   // end of initial node scan - - - - - - - - - - - - - - - - - - - - - - -

   bool cond1;

   if( GEZ( delx - deficit , EpsDfct ) ) {  // if no price change is
    quit = ( deficit < indef );             // possible, skip do loop
    cond1 = true;
    }
   else {
    // search along the ascent direction for the best price by checking the
    // slope of the dual cost at successive breakpoints; first, compute the
    // the distance to the next breakpoint

    CNumber delprc = nxtbrkpt( FIn + node , NxtIn , FOu + node , NxtOu );

    for(;;) {
     if( GTZ( deficit - delx , EpsDfct ) && ( delprc == Inf<CNumber>() ) ) {
      error_node = node;
      error_info = 5;
      status = MCFClass::kUnfeasible;
      return;
      }

     if( ! ETZ( delx , EpsDfct ) ) {
      // skip flow adjustment if there is no flow to modify
      // adjust the flow on the balanced arcs incident to node to maintain
      // complementary slackness after the price change

      Index j = nb_pos;
      Index_Set t_save = save;

      for( ; j-- ; ) {
       Index arc = *(t_save++);
       Index t2 = Endn[ arc ];
       FNumber f = X[ arc ];
       Dfct[ t2 ] += f;

       if( ! queue[ t2 ] ) {
        queue[ prvnde ] = t2;
        queue[ t2 ] = node;
        prvnde = t2;
        }

       U[ arc ] += f;
       X[ arc ] = 0;
       }

      for( j = m - nb_neg , t_save = save + m ; j-- ; ) {
       Index arc = *(--t_save);
       Index t2 = Startn[ arc ];
       FNumber f = U[ arc ];
       Dfct[ t2 ] += f;

       if( ! queue[ t2 ] ) {
        queue[ prvnde ] = t2;
        queue[ t2 ] = node;
        prvnde = t2;
        }

       X[ arc ] += f;
       U[ arc ] = 0;
       }

      deficit -= delx;

      }  // end if( delx == 0 )

     if( delprc == Inf<CNumber>() ) {
      quit = true;
      cond1 = false;
      break;
      }

     // node corresponds to a dual ascent direction: decrease the price of
     // node by delprc and compute the stepsize to the next breakpoint in
     // the dual cost

     pchange = true;
     delx = dascnt( node , delprc , FOu , NxtOu , FIn , NxtIn );

     if( GEZ( delx - deficit , EpsDfct ) ) {  // if no price change is
      quit = ( deficit < indef );              // possible, exit do loop
      cond1 = true;
      break;
      }
     }  // end for( ever )
    }  // end else( delx != deficit )

   // perform flow augmentation at node- - - - - - - - - - - - - - - - - - -

   if( cond1 ) {
    Index j = nb_pos;
    Index_Set t_save = save; 

    for( ; j-- ; ) {  // outgoing arcs from node
     Index arc = *(t_save++);
     Index t2 = Endn[ arc ];
     FNumber f2 = Dfct[ t2 ];

     if( GTZ( -f2 , EpsDfct ) ) {  // decrease the total deficit by
      quit = true;                 // decreasing flow of arc
      FNumber f = X[ arc ];
      FNumber dx = ( deficit < -f2 ? deficit : -f2 );
      if( f < dx )
       dx = f;

      deficit -= dx;
      Dfct[ t2 ] += dx;

      if( ! queue[ t2 ] ) {
       queue[ prvnde ] = t2;
       queue[ t2 ] = node;
       prvnde = t2;
       }

      X[ arc ] -= dx;
      U[ arc ] += dx;

      if( ETZ( deficit , EpsDfct ) )
       break;
      }
     }  // end for( j )

    for( j = m - nb_neg , t_save = save + m ; j-- ; ) {  // incoming arcs
     Index arc = *(--t_save);                            // into node
     Index t2 = Startn[ arc ];
     FNumber f2 = Dfct[ t2 ];

     if( GTZ( -f2 , EpsDfct ) ) {  // decrease the total deficit by
      quit = true;                 // increasing flow of arc
      FNumber f = U[ arc ];
      FNumber dx = ( deficit < -f2 ? deficit : -f2 );
      if( f < dx )
       dx = f;

      deficit -= dx;
      Dfct[ t2 ] += dx;

      if( ! queue[ t2 ] ) {
       queue[ prvnde ] = t2;
       queue[ t2 ] = node;
       prvnde = t2;
       }

      X[ arc ] += dx;
      U[ arc ] -= dx;

      if( ETZ( deficit , EpsDfct ) )
       break;
      }
     }  // end for( j )
    }  // end if( cond1 )

   Dfct[ node ] = deficit;

   // reconstruct the linked list of balanced arcs incident to this node:
   // for each adjacent node, add any newly balanced arcs to the list, but do
   // not bother removing formerly balanced ones (they will be removed the
   // next time each adjacent node is scanned)

   if( pchange )
    relist( node );

   }  // end if( posit ) - single node iteration for deficit > 0 node- - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else {
   // attempt a single node iteration from node with negative deficit - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   bool pchange = false;
   deficit = - deficit;
   FNumber indef = deficit;
   FNumber delx = svblncdarcs( node , tfstin , tnxtin , tfstou , tnxtou );

   // end of initial node scan - - - - - - - - - - - - - - - - - - - - - - -

   bool cond1;

   if( GEZ( delx - deficit , EpsDfct ) ) {  // if no price change is
    quit = ( deficit < indef );             // possible, skip do loop
    cond1 = true;
    }
   else {
    // search along the ascent direction for the best price by checking the
    // slope of the dual cost at successive breakpoints;  first, compute the
    // the distance to the next breakpoint

    CNumber delprc = nxtbrkpt( FOu + node , NxtOu , FIn + node , NxtIn );

    for(;;) {
     if( GTZ( deficit - delx , EpsDfct ) && ( delprc == Inf<CNumber>() ) ) {
      error_node = node;
      error_info = 6;
      status = MCFClass::kUnfeasible;
      return;
      }

     if( ! ETZ( delx , EpsDfct ) ) {
      // skip flow adjustment if there is no flow to modify
      // adjust the flow on the balanced arcs incident to node to maintain
      // complementary slackness after the price change

      Index j = nb_pos;
      Index_Set t_save = save;

      for( ; j-- ; ) {
       Index arc = *(t_save++);
       Index t2 = Startn[ arc ];
       FNumber f = X[ arc ];
       Dfct[ t2 ] -= f;

       if( ! queue[ t2 ] ) {
        queue[ prvnde ] = t2;
        queue[ t2 ] = node;
        prvnde = t2;
        }

       U[ arc ] += f;
       X[ arc ] = 0;
       }

      for( j = m - nb_neg , t_save = save + m ; j-- ; ) {
       Index arc = *(--t_save);
       Index t2 = Endn[ arc ];
       FNumber f = U[ arc ];
       Dfct[ t2 ] -= f;

       if( ! queue[ t2 ] ) {
        queue[ prvnde ] = t2;
        queue[ t2 ] = node;
        prvnde = t2;
        }

       X[ arc ] += f;
       U[ arc ] = 0;
       }

      deficit -= delx;
      }

     if( delprc == Inf<CNumber>() ) {
      quit = true;
      cond1 = false;
      break;
      }

     // node corresponds to a dual ascent direction: increase the price of
     // node by delprc and compute the stepsize to the next breakpoint in
     // the dual cost

     pchange = true;
     delx = dascnt( node , delprc , FIn , NxtIn , FOu , NxtOu );

     if( GEZ( delx - deficit , EpsDfct ) ) {  // if no price change is
      quit = ( deficit < indef );             // possible, exit do loop
      cond1 = true;
      break;
      }
     }  // end for( ever )
    }  // end else( delx != deficit )

   // perform flow augmentation at node- - - - - - - - - - - - - - - - - - -

   if( cond1 ) {
    Index j = nb_pos;
    Index_Set t_save = save;

    for( ; j-- ; ) {  // incoming arcs into node
     Index arc = *(t_save++);
     Index t2 = Startn[ arc ];
     FNumber f2 = Dfct[ t2 ];

     if( GTZ( f2 , EpsDfct ) ) {  // decrease the total deficit by
      quit = true;                // decreasing flow of arc
      FNumber f = X[ arc ];
      FNumber dx = ( deficit < f2 ? deficit : f2 );
      if( f < dx )
       dx = f;

      deficit -= dx;
      Dfct[ t2 ] -= dx;

      if( ! queue[ t2 ] ) {
       queue[ prvnde ] = t2;
       queue[ t2 ] = node;
       prvnde = t2;
       }

      X[ arc ] -= dx;
      U[ arc ] += dx;

      if( ETZ( deficit , EpsDfct ) )
       break;
      }
     }

    for( j = m - nb_neg , t_save = save + m ; j-- ; ) {  // outgoing arcs
     Index arc = *(--t_save);                            // from node
     Index t2 = Endn[ arc ];
     FNumber f2 = Dfct[ t2 ];

     if( GTZ( f2 , EpsDfct ) ) {  // decrease the total deficit by
      quit = true;                // increasing flow on arc
      FNumber f = U[ arc ];
      FNumber dx = ( deficit < f2 ? deficit : f2 );
      if( f < dx )
       dx = f;

      deficit -= dx;
      Dfct[ t2 ] -= dx;

      if( ! queue[ t2 ] ) {
       queue[ prvnde ] = t2;
       queue[ t2 ] = node;
       prvnde = t2;
       }

      X[ arc ] += dx;
      U[ arc ] -= dx;

      if( ETZ( deficit , EpsDfct ) )
       break;
      }
     }
    }

   Dfct[ node ] = -deficit;

   // reconstruct the linked list of balanced arcs incident to this node:
   // for each adjacent node, add any newly balanced arcs to the list, but do
   // not bother removing formerly balanced ones (they will be removed the
   // next time each adjacent node is scanned)

   if( pchange )
    relist( node );

   }  // end else( ! posit ) - single node iteration for deficit < 0 node - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ! quit ) && ( npass >= npasslim ) ) {
   // multinode iteration from node - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   pp( nmultinode );

   Switch = ( ndfct < tp );  // if number of nonzero deficit nodes is small,
                             // continue labelling until a flow augmentation
                             // is done
   bool nzerodfct = true;    // if the deficit of node has been zeroed

   // unmark nodes labeled earlier- - - - - - - - - - - - - - - - - - - - - -

   for( Index_Set tl = label + nlabel ; tl > label ; ) {
    Index i = *(--tl);
    mark[ i ] = scan[ i ] = false;
    }

   // initialize labelling- - - - - - - - - - - - - - - - - - - - - - - - - -

   nlabel = 1;
   mark[ *label = node ] = true;
   Prdcsr[ node ] = 0;

   // scan starting node- - - - - - - - - - - - - - - - - - - - - - - - - - -

   scan[ node ] = true;
   Index nscan = 1;
   FNumber dm = Dfct[ node ];
   FNumber delx = 0;

   Index j = nb_pos;
   Index_Set t_save = save; 

   for( ; j-- ; ) {
    Index arc = *(t_save++);
    Index t2 = ( posit ? Endn[ arc ] : Startn[ arc ] );

    if( ! mark[ t2 ] ) {
     Prdcsr[ t2 ] = arc;
     label[ nlabel++ ] = t2;
     mark[ t2 ] = true;
     delx += X[ arc ];
     }
    }

   for( j = m - nb_neg , t_save = save + m ; j-- ; ) {
    Index arc = *(--t_save);
    Index t2 = ( posit ? Startn[ arc ] : Endn[ arc ] );

    if( ! mark[ t2 ] ) {
     Prdcsr[ t2 ] = -arc;
     label[ nlabel++ ] = t2;
     mark[ t2 ] = true;
     delx += U[ arc ];
     }
    }

   // start scanning a labeled but unscanned node - - - - - - - - - - - - - -

   bool continua;

   do {
    Index i = label[ nscan++ ];

    // check to see if Switch needs to be set to true so to continue
    // scanning even after a price change

    Switch = Switch || ( ( nscan > ts ) && ( ndfct < ts ) );

    /* Scanning will continue until either an overestimate of the residual
       capacity across the cut corresponding to the scanned set of nodes
       (called delx) exceeds the absolute value of the total deficit of the
       scanned nodes (called dm), or an augmenting path is found. Arcs that
       are in the tree but are not balced are removed as part of the
       scanning process. */

    Index naugnod = 0;
    scan[ i ] = true;

    if( posit ) {  // scanning node i in case of positive deficit - - - - - -
     Index prvarc = 0;
     Index arc = tfstou[ i ];

     while( arc ) {  // arc is an outgoing arc from node i
      if( ETZ( RC[ arc ] , EpsCst ) ) {
       Index t2 = Endn[ arc ];
       if( ! mark[ t2 ] )  // t2 is not labeled
       {
        if( GTZ( X[ arc ], EpsFlw ) ) {
         if( LTZ( Dfct[ t2 ], EpsDfct ) )
          save[ naugnod++ ] = t2;

         Prdcsr[ t2 ] = arc;
         label[ nlabel++ ] = t2;
         mark[ t2 ] = true;
         delx += X[ arc ];
        } else {
         Prdcsr[ t2 ] = 0;
        }
       }

       prvarc = arc;
       arc = tnxtou[ arc ];
       }
      else {
       Index tmparc = arc;
       arc = tnxtou[ arc ];
       tnxtou[ tmparc ] = tmparc;

       if( prvarc )
        tnxtou[ prvarc ] = arc;
       else
        tfstou[ i ] = arc;
       }
      }

     prvarc = 0;
     arc = tfstin[ i ];

     while( arc ) {  // arc is an incoming arc into node i
      if( ETZ( RC[ arc ] , EpsCst ) ) {
       Index t2 = Startn[ arc ];
       if( ! mark[ t2 ] )  // t2 is not labeled
       {
        if( GTZ( U[ arc ], EpsFlw ) ) {
         if( LTZ( Dfct[ t2 ], EpsDfct ) )
          save[ naugnod++ ] = t2;

         Prdcsr[ t2 ] = -arc;
         label[ nlabel++ ] = t2;
         mark[ t2 ] = true;
         delx += U[ arc ];
        } else {
         Prdcsr[ t2 ] = 0;
        }
       }

       prvarc = arc;
       arc = tnxtin[ arc ];
       }
      else {
       Index tmparc = arc;
       arc = tnxtin[ arc ];
       tnxtin[ tmparc ] = tmparc;

       if( prvarc )
        tnxtin[ prvarc ] = arc;
       else
        tfstin[ i ] = arc;
       }
      }
     }
    else {  // scanning node i in case of negative deficit- - - - - - - - - -
     Index prvarc = 0;
     Index arc = tfstin[ i ];

     while( arc ) {  // arc is an incoming arc into node i
      if( ETZ( RC[ arc ] , EpsCst ) ) {
       Index t2 = Startn[ arc ];
       if( ! mark[ t2 ] )  // t2 is not labelled
       {
        if( GTZ( X[ arc ], EpsFlw ) ) {
         if( GTZ( Dfct[ t2 ], EpsDfct ) )
          save[ naugnod++ ] = t2;

         Prdcsr[ t2 ] = arc;
         label[ nlabel++ ] = t2;
         mark[ t2 ] = true;
         delx += X[ arc ];
        } else {
         Prdcsr[ t2 ] = 0;
        }
       }

       prvarc = arc;
       arc = tnxtin[ arc ];
       }
      else {
       Index tmparc = arc;
       arc = tnxtin[ arc ];
       tnxtin[ tmparc ] = tmparc;

       if( prvarc )
        tnxtin[ prvarc ] = arc;
       else
        tfstin[ i ] = arc;
       }
      }

     prvarc = 0;
     arc = tfstou[ i ];

     while( arc ) {  // arc is an outgoing arc from node i
      if( ETZ( RC[ arc ] , EpsCst ) ) {
       Index t2 = Endn[ arc ];
       if( ! mark[ t2 ] )  // t2 is not albelled
       {
        if( GTZ( U[ arc ], EpsFlw ) ) {
         if( GTZ( Dfct[ t2 ], EpsDfct ) )
          save[ naugnod++ ] = t2;

         Prdcsr[ t2 ] = -arc;
         label[ nlabel++ ] = t2;
         mark[ t2 ] = true;
         delx += U[ arc ];
        } else {
         Prdcsr[ t2 ] = 0;
        }
       }

       prvarc = arc;
       arc = tnxtou[ arc ];
       }
      else {
       Index tmparc = arc;
       arc = tnxtou[ arc ];
       tnxtou[ tmparc ] = tmparc;

       if( prvarc )
        tnxtou[ prvarc ] = arc;
       else
        tfstou[ i ] = arc;
       }
      }
     }  // end of scanning of node i for negative deficit case - - - - - - -

    SIndex arc = Prdcsr[ i ];

    if( arc > 0 )
     delx -= X[ arc ];
    else
     delx -= U[ -arc ];

    dm += Dfct[ i ];  // add deficit of node scanned to dm

    // check if the set of scanned nodes corresponds to a dual ascent
    // direction; if yes, performa a price adjustment step, otherwise
    // continue labelling

    continua = false;
    bool do_ascnt = true;

    if( nscan < nlabel )
     if( Switch || ( ( delx >= dm ) && ( delx >= -dm ) ) )
      do_ascnt = false;

    if( do_ascnt ) {
     // try a price change: since delx - abs( dm ) is an overestimate of the
     // ascent slope, we may occasionally try a direction that is not of
     // ascent - in this case, Ascnt() returns with quit = false and we
     // continue labeling nodes

     FNumber sdm;
     Index_Set Term1;
     Index_Set Term2;
     Index_Set F1;
     Index_Set F2;
     Index_Set Nxt1;
     Index_Set Nxt2;

     if( posit ) {
      sdm = dm;
      Term1 = Startn;
      Term2 = Endn;
      F1 = FOu;
      Nxt1 = NxtOu;
      F2 = FIn;
      Nxt2 = NxtIn;
      }
     else {
      sdm = -dm;
      Term1 = Endn;
      Term2 = Startn;
      F1 = FIn;
      Nxt1 = NxtIn;
      F2 = FOu;
      Nxt2 = NxtOu;
      }

     if( ! Ascnt( sdm , delx , nlabel , Switch , nscan , node , Term1 ,
                  Term2 , F1 , Nxt1 , F2 , Nxt2 ) ) {
      error_node = node;
      error_info = 7;
      status = MCFClass::kUnfeasible;
      return;
      }

     if( ! Switch )
      break;

     // store those newly labeled nodes to which - - - - - - - - - - - - - -
     // flow augmentation is possible- - - - - - - - - - - - - - - - - - - -

     Index_Set t_label = label + nscan;
     Index h = nlabel - nscan;
     naugnod = 0;

     if( posit )
      for( ; h-- ; ) {
       Index t2 = *(t_label++);
       if( LTZ( Dfct[ t2 ] , EpsDfct ) )
        save[ naugnod++ ] = t2;
       }
     else
      for( ; h-- ; ) {
       Index t2 = *(t_label++);
       if( GTZ( Dfct[ t2 ] , EpsDfct ) )
        save[ naugnod++ ] = t2;
       }

     }  // end if( do_ascnt )

    if( naugnod > 0 ) {  // if flow augmentation is possible- - - - - - - - -
                         // do the augmentation - - - - - - - - - - - - - - -
     Index_Set tsv = save;
     for( Index h = naugnod ; h-- ; ) {
      num_augm++;
      Index augnod = *(tsv++);
      Index node_n;
      Index node_p;
      Index_Set Term1;
      Index_Set Term2;

      if( posit ) {
       node_n = node;
       node_p = augnod;
       Term1 = Startn;
       Term2 = Endn;
       }
      else {
       node_p = node;
       node_n = augnod;
       Term1 = Endn;
       Term2 = Startn;
       }

      AugFlow( augnod , node , node_p , node_n , Term1 , Term2 );

      if( ETZ( Dfct[ node ] , EpsDfct ) ) {
       nzerodfct = false;
       break;
       }

      if( ! ETZ( Dfct[ augnod ] , EpsDfct ) )
       Switch = false;
      }
     }
    else
     continua = true;

    // if node has still nonzero deficit and all newly labeled nodes have
    // the same sign for their deficit as node, we continue labeling
    // only when flow augmentation is done relatively infrequently

    } while( nzerodfct && 
             ( continua || ( Switch && ( iter > it_aug * num_augm ) ) ) );
   }

  if( MaxIter && ( iter > MaxIter ) ) {  // iterations limit
   status = kStopped;
   break;
   }

  if( MCFt && MaxTime && ( MCFt->Read() > MaxTime ) ) {  // time limit
   status = kStopped;
   break;
   }
  }  // end for( ever ) - main loop ends here - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MCFt )
  MCFt->Stop();

 }  // end( RelaxIV::SolveMCF )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void RelaxIV::MCFGetX( FRow F , Index_Set nms , cIndex strt , Index stp )
{
 if( stp > m )
   stp = m;

 if( nms ) {
  cFRow tX = X + strt;
  for( Index i = strt ; i < stp ; i++ ) {
   cFNumber ttX = *(++tX);
   if( GTZ( ttX , EpsFlw ) ) {
    *(F++) = ttX;
    *(nms++) = i;
    }

   *nms = Inf<Index>();
   }
  }
 else {
  cFRow tX = X + stp;
  for( F += stp - strt ; tX > X + strt ; )
   *(--F) = *(tX--);
  }
 }

/*--------------------------------------------------------------------------*/

MCFClass::cFRow RelaxIV::MCFGetX( void )
{
 return( X + 1 );
 }

/*--------------------------------------------------------------------------*/

void RelaxIV::MCFGetRC( CRow CR , cIndex_Set nms , cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(CR++) = RC[ ++h ];
  }
 else {
  if( stp > m )
   stp = m;

  CRow tRC = RC + stp;
  for( CR += stp - strt ; tRC > RC + strt ; )
   *(--CR) = *(tRC--);
  }
 }

/*--------------------------------------------------------------------------*/

MCFClass::cCRow RelaxIV::MCFGetRC( void )
{
 return( RC + 1 );
 }

/*--------------------------------------------------------------------------*/

MCFClass::CNumber RelaxIV::MCFGetRC( cIndex i )
{
 return( RC[ i + 1 ] );
 }

/*--------------------------------------------------------------------------*/

void RelaxIV::MCFGetPi( CRow P , cIndex_Set nms , cIndex strt , Index stp )
{
 #if( P_ALLOC )
  if( PiOwnr != this )
   cmptprices();

  if( nms ) {
   while( *nms < strt )
    nms++;

   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(P++) = Pi[ ++h ];
   }
  else {
   if( stp > n )
    stp = n;

   CRow tPi = Pi + stp;
   for( P += stp - strt ; tPi > Pi + strt ; )
    *(--P) = *(tPi--);
   }
 #else
  Pi = P - 1;

  cmptprices();

  Pi++;

  if( nms ) {
   while( *nms < strt )
    nms++;

   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(P++) = Pi[ h ];
   }
  else
   if( strt ) {
    if( stp > n )
     stp = n;

    CRow tP = P;
    for( stp -= strt , P += strt ; stp-- ; )
     *(tP++) = *(P++);
    }
 #endif

 }  // end( RelaxIV::MCFGetPi( CRow ) )

/*--------------------------------------------------------------------------*/

MCFClass::cCRow RelaxIV::MCFGetPi( void )
{
 #if( P_ALLOC )
  if( PiOwnr != this )
   cmptprices();

  return( Pi + 1 );
 #else
  return( NULL );
 #endif

 }  // end( RelaxIV::MCFGetPi( void ) )

/*--------------------------------------------------------------------------*/

MCFClass::FONumber RelaxIV::MCFGetFO( void )
{
 if( status == MCFClass::kUnfeasible )
  return( Inf<FONumber>() );
 else {
  if( FO == Inf<FONumber>() ) {
   cFRow tX = X + m;
   cCRow tC = C + m;
   for( FO = 0 ; tX > X ; )
    FO += *(tX--) * (*(tC--));
   }

  return( FO );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------ SPECIALIZED INTERFACE ---------------------------*/
/*--------------------------------------------------------------------------*/

MCFClass::MCFStatePtr RelaxIV::MCFGetState( void )
{
 RIVState *S = new RIVState( m );

 CRow tRC = RC + m;
 CRow tSR = (S->RedCost) + m;
 FRow tSF = (S->Flow) + m;

 for( FRow tX = X + m ; tX > X ; ) {
  *(--tSF) = *(tX--);
  *(--tSR) = *(tRC--);
  }

 return( S );

 }  // end( MCFGetState )

/*--------------------------------------------------------------------------*/

void RelaxIV::MCFPutState( MCFClass::MCFStatePtr S )
{
 RIVState *RS = dynamic_cast<RIVState*>( S );
 if( ! RS )
  return;

 // complementary slackness conditions and bounds must be verified- - - - - -

 cFRow tCap = Cap + m;
 cFRow tSF = (RS->Flow) + m;
 cCRow tSRC = (RS->RedCost) + m;
 for( ; tCap > Cap ; ) {
  FNumber ttCap = *(tCap--);
  cFNumber ttSF = *(--tSF);
  cCNumber ttSRC = *(--tSRC);
  if( ttSRC == Inf<CNumber>() )
   continue;

  ttCap -= ttSF;

  if( LTZ( ttCap , EpsFlw ) || LTZ( ttSF , EpsFlw ) )
   return;

  if( GTZ( ttSRC , EpsCst ) && GTZ( ttSF , EpsFlw ) )
   return;

  if( LTZ( ttSRC , EpsCst ) && GTZ( ttCap , EpsFlw ) )
   return;
  }

 #if( P_ALLOC )
 {
  // check that RC[ i , j ] = C[ i , j ] + P[ i ] - P[ j ]- - - - - - - - - -

  CRow tRC = RC;         // save Reduced Costs pointer
  RC = RS->RedCost - 1;  // temporarily use the new presumed RC
  cmptprices();          // compute Pi[] with the new RCs

  cCRow SRC = RC;
  RC = tRC;             // restore the current prices

  cCRow tPi = Pi + n;
  for( Index_Set tou = FOu + n ; tou > FOu ; ) {
   cCNumber Pii = *(tPi--);
   Index arc = *(tou--);

   while( arc ) {
    cCNumber Dlt = Pi[ Endn[ arc ] ] - Pii - C[ arc ] + SRC[ arc ];
    if( ! ETZ( Dlt , EpsCst ) )
     return;

    arc = NxtOu[ arc ];
    }
   }
  }
 #endif

 // correct the internal state of RelaxIV - - - - - - - - - - - - - - - - - -

 FRow tDfct = Dfct + n;
 for( FRow tB = B + n ; tB > B ; )
  *(tDfct--) = *(tB--);

 FRow tX = X + m;
 FRow tU = U + m;
 CRow tRC = RC + m;
 Index_Set tEn = Endn + m;
 Index_Set tSn = Startn + m;
 for( tSF = RS->Flow + m , tSRC = RS->RedCost + m , tCap = Cap + m ;
      tU > U ; ) {
  *(tRC--) = *(--tSRC);

  FNumber tXi = *(tX--) = *(--tSF);
  *(tU--) = *(tCap--) - tXi;
  Dfct[ *(tSn--) ] += tXi;
  Dfct[ *(tEn--) ] -= tXi;
  }

 init_tree();

 }  // end( MCFPutState )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void RelaxIV::MCFArcs( Index_Set Startv , Index_Set Endv , cIndex_Set nms ,
                       cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;
  
  if( Startv && Endv )
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    *(Startv++) = Startn[ ++h ] - USENAME0;
    *(Endv++) = Endn[ h ] - USENAME0;
    }
  else {
   if( Startv )
    for( Index h ; ( h = *(nms++) ) < stp ; )
     *(Startv++) = Startn[ ++h ] - USENAME0;

   if( Endv )
    for( Index h ; ( h = *(nms++) ) < stp ; )
     *(Endv++) = Endn[ ++h ] - USENAME0;
   }
  }
 else {
  if( stp > m )
   stp = m;

  if( Endv ) {
   Index_Set tE = Endn + stp;
   for( Endv += stp - strt ; tE > Endn + strt ; )
    *(--Endv) = *(tE--) - USENAME0;
   }

  if( Startv ) {
   Index_Set tS = Startn + stp;
   for( Startv += stp - strt ; tS > Startn + strt ; )
    *(--Startv) = *(tS--) - USENAME0;
   }
  }
 }  // end( RelaxIV::MCFArcs )

/*--------------------------------------------------------------------------*/

void RelaxIV::MCFCosts( CRow Costv , cIndex_Set nms , cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(Costv++) = C[ ++h ];
  }
 else {
  if( stp > m )
   stp = m;

  CRow tC = C + stp;
  for( Costv += stp - strt ; tC > C + strt ; )
   *(--Costv) = *(tC--);
  }
 }  // end( RelaxIV::MCFCosts )

/*--------------------------------------------------------------------------*/

void RelaxIV::MCFUCaps( FRow UCapv , cIndex_Set nms , cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(UCapv++) = Cap[ ++h ];
  }
 else {
  if( stp > m )
   stp = m;

  FRow tCap = Cap + stp;
  for( UCapv += stp - strt ; tCap > Cap + strt ; )
   *(--UCapv) = *(tCap--);
  }
 }  // end( RelaxIV::MCFUCaps )
  
/*--------------------------------------------------------------------------*/

void RelaxIV::MCFDfcts( FRow Dfctv , cIndex_Set nms , cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(Dfctv++) = B[ ++h ];
  }
 else {
  if( stp > n )
   stp = n;

  FRow tB = B + stp;
  for( Dfctv += stp - strt ; tB > B + strt ; )
   *(--Dfctv) = *(tB--);
  }
 }  // end( RelaxIV::MCFDfcts )

/*--------------------------------------------------------------------------*/

void RelaxIV::WriteMCF( ostream &oStrm , int frmt )
{
 #if( ( Ctype == REAL_TYPE ) || ( Ftype == REAL_TYPE ) )
  oStrm.precision( 12 );
 #endif

 switch( frmt ) {
  case( kCLP ):  // LP format - - - - - - - - - - - - - - - - - - - - - - - -
                 //-  - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  {   // writing the objective function - - - - - - - - - - - - - - - - - - -

    cCRow tC = C;
    cCRow tRC = RC;
    short chrcnt = 0;
    for( Index i = 0 ; i < m ; i++ ) {
     CNumber ttC = *(++tC);
     if( *(++tRC) == Inf<CNumber>() )
      continue;

     if( chrcnt > 230 ) {
      oStrm << endl;
      chrcnt = 20;
      }
     else
      chrcnt += 20;
 
     if( ttC >= 0 )
      oStrm << '+';

     oStrm << ttC << 'x' << i;
     }
    }

   oStrm << endl << "subject to" << endl;

   {  // writing the flow conservation constraints- - - - - - - - - - - - - -

    cFRow tB = B;
    cIndex_Set tFOu = FOu;
    cIndex_Set tFIn = FIn;
    for( Index i = n ; i-- ; ) {
     short chrcnt = 0;
     Index arc = *(++tFOu);

     while( arc ) {
      if( chrcnt > 240 ) {
       oStrm << endl;
       chrcnt = 11;
       }
      else
       chrcnt += 11;

      oStrm << "-x" << arc - 1;
      arc = NxtOu[ arc ];
      }

     arc = *(++tFIn);

     while( arc ) {
      if( chrcnt > 240 ) {
       oStrm << endl;
       chrcnt = 11;
       }
      else
       chrcnt += 11;

      oStrm << "+x" << arc - 1;
      arc = NxtIn[ arc ];
      }

     if( chrcnt > 240 )
      oStrm << endl;

     oStrm << '=' << *(++tB) << endl;  
     }
    }

   oStrm << "bounds" << endl;

   {  // writing the bounds - - - - - - - - - - - - - - - - - - - - - - - - -

    cCRow tRC = RC;
    FRow tCap = Cap;
    for( Index i = 0 ; i < m ; tCap++ , i++ ) {
     cFNumber ttC = *(++tCap);
     if( *(++tRC) == Inf<CNumber>() )
      oStrm << "0 <= x" << i << " <= " << ttC << endl;
     }
    }

   oStrm << "end" << endl;
   break;

  case( kRIV ):  // RelaxIV format- - - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   oStrm << n << "\t" << m << endl << endl;

   Index i;
   for( i = 0 ; i++ < m ; ) {
    oStrm << Startn[ i ] << "\t" << Endn[ i ] << "\t" << U[ i ] << "\t";

    if( RC[ i ] < Inf<CNumber>() )
     oStrm << RC[ i ] << endl;
    else
     oStrm << "+INF" << endl;
    }

   oStrm << endl;

   for( i = 0 ; i < n ; )
    oStrm << Dfct[ ++i ] << endl;

   break;

  default:       // any other format- - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   MCFClass::WriteMCF( oStrm , frmt );
  }
 }  // end( MCFClass::WriteMCF )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/*----- Changing the costs, deficits and upper capacities of the (MCF) -----*/
/*--------------------------------------------------------------------------*/

void RelaxIV::ChgCosts( cCRow NCost , cIndex_Set nms , cIndex strt ,
			Index stp )
{
 if( nms )
  while( *nms < strt ) {
   nms++;
   NCost++;
   }

 if( stp > m )
  stp = m;

 if( status || ( ! Senstv ) ) {
  CRow tC = C;
  if( nms ) {
   Index h;
   for( tC++ ; ( h = *(nms++) ) < stp ; )
    tC[ h ] = *(NCost++);
   }
  else
   for( tC += stp , NCost += stp - strt ; tC > C + strt ; )
    *(tC--) = *(--NCost);

  status = MCFClass::kUnSolved;
  }
 else
  if( nms )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    chgcsti( ++h , *(NCost++) );
  else
   for( Index h = strt ; h < stp ; )
    chgcsti( ++h , *(NCost++) );

 }  // end( RelaxIV::ChgCosts( some / all ) )

/*--------------------------------------------------------------------------*/

void RelaxIV::ChgCost( Index arc , cCNumber NCost )
{
 if( status || ( ! Senstv ) ) {
  C[ ++arc ] = NCost;
  status = MCFClass::kUnSolved;
  }
 else
  chgcsti( ++arc , NCost );

 }  // end( RelaxIV::ChgCost( one ) )

/*--------------------------------------------------------------------------*/

void RelaxIV::ChgDfcts( cFRow NDfct , cIndex_Set nms , cIndex strt ,
			Index stp )
{
 if( nms )
  while( *nms < strt ) {
   nms++;
   NDfct++;
   }

 if( stp > n )
  stp = n;

 FRow tB = B + 1;

 if( status || ( ! Senstv ) ) {
  if( nms )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    tB[ h ] = *(NDfct++);
  else
   for( tB += stp , NDfct += stp - strt ; tB-- > B + strt ; )
    *tB = *(--NDfct);

  status = MCFClass::kUnSolved;
  }
 else {
  FRow tDfct = Dfct + 1;
  if( nms )
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    FNumber NDh = *(NDfct++);
    tDfct[ h ] += NDh - tB[ h ];
    tB[ h ] = NDh;
    }
  else
   for( tDfct += stp , NDfct += stp - strt , tB += stp ; tB-- > B + strt ; ) {
    *(--tDfct) += *(--NDfct) - *tB;
    *tB = *NDfct;
    }
  }
 }  // end( RelaxIV::ChgDfcts( some / all ) )

/*--------------------------------------------------------------------------*/

void RelaxIV::ChgDfct( Index nod , cFNumber NDfct )
{
 nod++;
 if( status || ( ! Senstv ) ) {
  B[ nod ] = NDfct;
  status = MCFClass::kUnSolved;
  }
 else {   
  Dfct[ nod ] += NDfct - B[ nod ];
  B[ nod ] = NDfct;
  }
 }  // end( RelaxIV::ChgDfcts( one ) )

/*--------------------------------------------------------------------------*/

void RelaxIV::ChgUCaps( cFRow NCap , cIndex_Set nms , cIndex strt , Index stp )
{
 if( nms )
  while( *nms < strt ) {
   nms++;
   NCap++;
   }

 if( stp > m )
  stp = m;

 if( status || ( ! Senstv ) ) {
  FRow tCap = Cap;
  if( nms ) {
   Index h;
   for( tCap++ ; ( h = *(nms++) ) < stp ; )
    tCap[ h ] = *(NCap++);
   }
  else
   for( tCap += stp , NCap += stp - strt ; tCap > Cap + strt ; )
    *(tCap--) = *(--NCap);

  status = MCFClass::kUnSolved;
  }
 else
  if( nms )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    chgcapi( ++h , *(NCap++) );
  else
   for( Index h = strt ; h < stp ; )
    chgcapi( ++h , *(NCap++) );

 }  // end( RelaxIV::ChgUCaps( some / all ) )

/*--------------------------------------------------------------------------*/

void RelaxIV::ChgUCap( Index arc , cFNumber NCap )
{
 if( status || ( ! Senstv ) ) {
  Cap[ ++arc ] = NCap;
  status = MCFClass::kUnSolved;
  }
 else   
  chgcapi( ++arc , NCap );

 }  // end( RelaxIV::ChgUCaps( one ) )

/*--------------------------------------------------------------------------*/
/*------------------ Modifying the structure of the graph ------------------*/ 
/*--------------------------------------------------------------------------*/

void RelaxIV::CloseArc( cIndex name )
{
 #if( DYNMC_MCF_RIV )
  delarci( name + 1 );
 #else
  throw(
   MCFException( "RelaxIV::CloseArc() not implemented if DYNMC_MCF_RIV < 1"
                 ) );
 #endif

 }  // end( RelaxIV::CloseArc )

/*--------------------------------------------------------------------------*/

void RelaxIV::DelNode( cIndex name )
{
 #if( DYNMC_MCF_RIV )
  Index node = name + USENAME0;
  Index arc = FOu[ node ];

  while( arc ) {
   CloseArc( arc );
   arc = FOu[ node ];
   }

  arc = FIn[ node ];

  while( arc ) {
   CloseArc( arc );
   arc = FIn[ node ];
   }
  
  Dfct[ node ] = 0;

  if( node == n )
   do
    n--;
   while( FOu[ n ] == FIn[ n ] );

  status = MCFClass::kUnSolved;
 #else
  throw(
   MCFException( "RelaxIV::DelNode() not implemented if DYNMC_MCF_RIV < 1"
                 ) );
 #endif

 }  // end( RelaxIV::DelNode )

/*--------------------------------------------------------------------------*/

void RelaxIV::OpenArc( cIndex name )
{
 #if( DYNMC_MCF_RIV > 1 )
  addarci( name + 1 );
 #else
  throw(
   MCFException( "RelaxIV::OpenArc() not implemented if DYNMC_MCF_RIV < 2"
                 ) );
 #endif

 }  // end( RelaxIV::UnsetArcFree )

/*--------------------------------------------------------------------------*/

MCFClass::Index RelaxIV::AddNode( cFNumber aDfct )
{
 #if( DYNMC_MCF_RIV > 1 )
  if( n == nmax )
   return( Inf<Index>() );

  n++;

  B[ n ] = aDfct;
  FOu[ n ] = FIn[ n ] = 0;

  if( status || ( ! Senstv ) )
   status = MCFClass::kUnSolved;
  else {
   Dfct[ n ] = aDfct;
   tfstou[ n ] = tfstin[ n ] = 0;

   #if( P_ALLOC )
    if( PiOwnr == this )
     Pi[ n ] = 0;
   #endif
   }

  return( n - USENAME0 );
 #else
  throw(
   MCFException( "RelaxIV::AddNode() not implemented if DYNMC_MCF_RIV < 2"
                 ) );

  return( Inf<Index>() );
 #endif

 }  // end( RelaxIV::AddNode )

/*--------------------------------------------------------------------------*/

void RelaxIV::ChangeArc( cIndex name , cIndex nSN , cIndex nEN )
{
 #if( DYNMC_MCF_RIV > 2 )
  Index arc = name + 1;
  if( RC[ arc ] < Inf<CNumber>() )  // the arc is currently open- - - - - - -
   if( status || ( ! Senstv ) ) {  // no need to reoptimize - - - - - - - - -
    delarci( arc );

    if( nSN < Inf<Index>() )
     Startn[ arc ] = nSN + USENAME0;

    if( nEN < Inf<Index>() )
     Endn[ arc ] =  nEN + USENAME0;

    addarci( arc );
    }
   else {  // reoptimization is required- - - - - - - - - - - - - - - - - - -
    if( PiOwnr != this )  // compute the dual prices
     cmptprices();

    cCNumber RCa = RC[ arc ] = C[ arc ] + Pi[ nSN ] - Pi[ nEN ];
    cFNumber Ua = Cap[ arc ];

    if( nSN < Inf<Index>() ) {  // the start node changes
     Dfct[ Startn[ arc ] ] -= X[ arc ];  // update the deficit

     Index arc1 = FOu[ Startn[ arc ] ];  // update the FS
     if( arc1 == arc )
      FOu[ Startn[ arc ] ] = NxtOu[ arc1 ];
     else {
      Index arc2;
      do {
       arc2 = arc1;
       arc1 = NxtOu[ arc1 ];
       } while( arc1 != arc );

      NxtOu[ arc2 ] = NxtOu[ arc ];
      }

     if( tnxtou[ arc ] != arc ) {  // update the "restricted" FS
      if( ( arc1 = tfstou[ Startn[ arc ] ] ) == arc )
       tfstou[ Startn[ arc ] ] = tnxtou[ arc ];
      else {
       Index arc2;
       do {
        arc2 = arc1;
        arc1 = tnxtou[ arc1 ];
        } while( arc1 != arc );

       tnxtou[ arc2 ] = tnxtou[ arc ];
       }

      tnxtou[ arc ] = arc;
      }

     cIndex sn = nSN + USENAME0;  // update Startn[]
     Startn[ arc ] = sn;
     NxtOu[ arc ] = FOu[ sn ];
     FOu[ sn ] = arc;

     if( LTZ( RCa , EpsCst ) )
      Dfct[ sn ] += Ua;

     if( ETZ( RCa , EpsCst ) ) {
      tnxtou[ arc ] = tfstou[ sn ];
      tfstou[ sn ] = arc;
      }
     }

    if( nEN < Inf<Index>() ) {  // the end node changes
     Dfct[ Endn[ arc ] ] += X[ arc ];  // update the deficit

     Index arc1 = FIn[ Endn[ arc ] ];
     if( arc1 == arc )  // update the BS
      FIn[ Endn[ arc ] ] = NxtIn[ arc1 ];
     else {
      Index arc2;
      do {
       arc2 = arc1;
       arc1 = NxtIn[ arc1 ];
       } while( arc1 != arc );

      NxtIn[ arc2 ] = NxtIn[ arc ];
      }

     if( tnxtin[ arc ] != arc ) {  // update the "restricted" BS
      if( ( arc1 = tfstin[ Endn[ arc ] ] ) == arc )
       tfstin[ Endn[ arc ] ] = tnxtin[ arc ];
      else {
       Index arc2;
       do {
        arc2 = arc1;
        arc1 = tnxtin[ arc1 ];
        } while( arc1 != arc );

       tnxtin[ arc2 ] = tnxtin[ arc ];
       }

      tnxtin[ arc ] = arc;
      }

     cIndex en = nEN + USENAME0;  // update Endn[]
     Endn[ arc ] = en;
     NxtIn[ arc ] = FIn[ en ];
     FIn[ en ] = arc;

     if( LTZ( RCa , EpsCst ) )
      Dfct[ en ] -= Ua;

     if( ETZ( RCa , EpsCst ) ) {
      tnxtin[ arc ] = tfstin[ en ];
      tfstin[ en ] = arc;
      }
     }

   if( LTZ( RCa , EpsCst ) ) {
    X[ arc ] = Ua;
    U[ arc ] = 0;
    }
   else {
    X[ arc ] = 0;
    U[ arc ] = Ua;
    }
   }      // end( else( reoptimization is required ) )- - - - - - - - - - - -
  else {  // the arc is currently closed- - - - - - - - - - - - - - - - - - -
          //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( nSN < Inf<Index>() )
    Startn[ arc ] = nSN + USENAME0;

   if( nEN < Inf<Index>() )
    Endn[ arc ] =  nEN + USENAME0;
   }
 #else
  throw(
   MCFException( "RelaxIV::ChangeArc() not implemented if DYNMC_MCF_RIV < 3"
                 ) );
 #endif

 }  // end( RelaxIV::ChangeArc )

/*--------------------------------------------------------------------------*/

void RelaxIV::DelArc( cIndex name )
{
 #if( DYNMC_MCF_RIV > 2 )
  Index arc = name;
  delarci( ++arc );

  if( arc == m ) {
   while( arc && ( RC[ arc ] == Inf<CNumber>() ) )
    arc--;
 
   m = arc;
   }

  Startn[ arc ] = Inf<Index>();
  Endn[ arc ] = ffp;
  ffp = arc;
 #else
  throw(
   MCFException( "RelaxIV::DelArc() not implemented if DYNMC_MCF_RIV < 3"
                 ) );
 #endif

 }  // end( RelaxIV::DelArc )

/*--------------------------------------------------------------------------*/

MCFClass::Index RelaxIV::AddArc( cIndex Start , cIndex End , cFNumber aU ,
				 cCNumber aC )
{
 #if( DYNMC_MCF_RIV > 2 )
  // select position - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  Index arc = ffp;
  while( arc ) {   
   ffp = Endn[ ffp ];

   if( arc <= m )
    break;

   arc = ffp;
   }

  if( ! arc ) {
   if( m < mmax ) {
    arc = ++m;
   } else {
    return ( Inf< Index >() );
   }
  }

 // insert new arc in position arc - - - - - - - - - - - - - - - - - - - - -

  C[ arc ] = aC;
  Cap[ arc ] = aU;
  Endn[ arc ] = End + USENAME0;
  Startn[ arc ] = Start + USENAME0;

  addarci( arc );

  return( arc - 1 );
 #else
  throw(
   MCFException( "RelaxIV::AddArc() not implemented if DYNMC_MCF_RIV < 3"
                 ) );

  return( Inf<Index>() );
 #endif

 }  // end( RelaxIV::AddArc )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

RelaxIV::~RelaxIV()
{
 if( ! --InstCntr ) {  // deallocating static members - - - - - - - - - - - -
  maxnmax = maxmmax = 0;

  #if( AUCTION )
   delete[] ++NxtpushB;
   delete[] ++NxtpushF;
  #endif

  delete[] save;
  save = NULL;

  #if( AUCTION )
   delete[] ++SB_arc;
   delete[] ++FpushB;
   delete[] ++extend_arc;
   delete[] ++SB_level;
  #endif

  #if( P_ALLOC )
   delete[] ++Pi;
  #endif

  delete[] ++DDNeg;
  delete[] ++DDPos;
  delete[] ++queue;
  delete[] ++scan;
  delete[] ++mark;
  mark = NULL;
  delete[] label;
  delete[] ++Prdcsr;
  }

 MemDeAlloc();  // deallocate all the rest- - - - - - - - - - - - - - - - - -

 }  // end( ~RelaxIV )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void RelaxIV::init_tree( void )
{
 Index_Set tou = tfstou + n;
 Index_Set tin = tfstin + n;
 for( ; tin > tfstin ; )
  *(tou--) = *(tin--) = 0;

 Index arc = 0;
 cCRow tRC = RC;
 for( tou = tnxtou , tin = tnxtin ; arc++ < m ; ) {
  cCNumber RCarc = *(++tRC);  // beware of ETZ()
  if( ETZ( RCarc , EpsCst ) ) {
   Index i = Startn[ arc ];
   *(++tou) = tfstou[ i ];
   tfstou[ i ] = arc;

   *(++tin) = tfstin[ i = Endn[ arc ] ];
   tfstin[ i ] = arc;
   continue;
   }

  *(++tin) = *(++tou) = arc;
  }
 }  // end( init_tree )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::init_standard( void )
{
 // reset B and initialize the directional derivative for each coordinate - -
 // (note that DDPos and DDNeg are only used inside init_standard())- - - - -

 FRow tDfct = Dfct + n;
 FRow tDDPos = DDPos + n;
 FRow tDDNeg = DDNeg + n;
 for( cFRow tB = B + n ; tB > B ; ) {
  cFNumber ttB = *(tB--);
  *(tDfct--) = ttB;
  *(tDDPos--) = ttB;
  *(tDDNeg--) = - ttB;
  }

 // initialize the arc flows X and the reduced capacities U to satisfy- - - -
 // complementary slackness with the prices; meanwhile, compute the - - - - -
 // directional derivative for each coordinate and the actual deficits- - - -

 FRow tX = X + m;
 FRow tU = U + m;
 cCRow tC = C + m;
 CRow tRC = RC + m;
 cFRow tCap = Cap + m;
 for( ; tC > C ; tX-- , tU-- , tRC-- , tCap-- , tC-- ) {
  #if( ( ! SAME_GRPH_RIV ) || DYNMC_MCF_RIV )
   if( *tRC == Inf<CNumber>() )
    continue;
  #endif

  cFNumber f = *tCap;
  cCNumber RCi = *tRC = *tC;

  if( LEZ( RCi , EpsCst ) ) {
   cIndex i = tX - X;
   cIndex si = Startn[ i ];
   cIndex ei = Endn[ i ];
   DDPos[ si ] += f;
   DDNeg[ ei ] += f;

   if( LTZ( RCi , EpsCst ) ) {
    *tU = 0;
    *tX = f;
    Dfct[ si ] += f;
    Dfct[ ei ] -= f;
    DDNeg[ si ] -= f;
    DDPos[ ei ] -= f;
    continue;
    }
   }

  *tU = f;
  *tX = 0;
  }

 // make 2 or 3 (depending on the density of the network) passes through- - -
 // all nodes, performing only single node relaxation iterations- - - - - - -

 for( char npass = ( m > n * maxdns ? 2 : 3 ) ; npass-- ; ) {
  Index_Set tFOu = FOu + n;
  Index_Set tFIn = FIn + n;
  for( tDfct = Dfct + n , tDDPos = DDPos + n , tDDNeg = DDNeg + n ;
       tDfct > Dfct ; tDfct-- , tDDPos-- , tDDNeg-- , tFOu-- , tFIn-- )
   if( ! ETZ( *tDfct , EpsDfct ) ) {
    if( LEZ( *tDDPos , EpsDfct ) ) {
     // compute delprc, the stepsize to the next breakpoint in the dual cost
     // as the price of node is increased: since the reduced cost of all
     // outgoing (resp. incoming) arcs will decrease (resp. increase) as the
     // price of node is increased, the next breakpoint is the minimum of the
     // positive reduced cost on outgoing arcs and of the negative reduced
     // cost on incoming arcs

     CNumber delprc = nxtbrkpt( tFOu , NxtOu , tFIn , NxtIn );

     if( delprc == Inf<CNumber>() ) {           // if no breakpoint is left
      if( ETZ( *tDDPos , EpsDfct ) )  // and dual ascent is still possible
       continue;                       // try another dual ascent

      status = MCFClass::kUnfeasible;  // else the problem is unfeasible
      error_node = tDfct - Dfct;
      error_info = 3;
      return;
      }

     for(;;) {
      // delprc is the stepsize to next breakpoint: increase price of node
      // by delprc and compute the stepsize to the next breakpoint

      CNumber nxtbrk = Inf<CNumber>();

      Index arc = *tFOu;
      while( arc ) {  // look at all arcs out of the current node
       CNumber trc = mvflw1( arc , tDfct , tDDNeg , Endn , U , X );

       decrsRC( arc , trc , delprc , nxtbrk , tDDPos , DDNeg , Endn );

       arc = NxtOu[ arc ];
       }

      arc = *tFIn;
      while( arc ) {  // look at all arcs into the current node
       CNumber trc = mvflw1( arc , tDfct , tDDNeg , Startn , X , U );

       incrsRC( arc , trc , delprc , nxtbrk , tDDPos , DDNeg , Startn );

       arc = NxtIn[ arc ];
       }

      // if price of current node can be increased further without decreasing
      // the dual cost (even if dual cost doesn't increase ), do another
      // iteration to increase the price

      if( LEZ( *tDDPos , EpsDfct ) && ( nxtbrk < Inf<CNumber>() ) )
       delprc = nxtbrk;
      else
       break;
      }         // end for( ever )
     }
    else {
     if( LEZ( *tDDNeg , EpsDfct ) ) {
      // Compute delprc, the stepsize to the next breakpoint in the dual cost
      // as the price of node is decreased. Since the reduced cost of all
      // outgoing (resp. incoming) arcs will increase (resp. decrease) as the
      // price  of node is decreased, the next breakpoint is the minimum of
      // the negative reduced cost on outgoing arcs and of the positive
      // reduced cost on incoming arcs

      CNumber delprc = nxtbrkpt( tFIn , NxtIn , tFOu , NxtOu );

      if( delprc == Inf<CNumber>() ) {           // if no breakpoint is left
       if( ETZ( *tDDNeg , EpsDfct ) )  // and dual ascent is still possible
        continue;                       // try another dual ascent

       status = MCFClass::kUnfeasible;  // else the problem is unfeasible
       error_node = tDfct - Dfct;
       error_info = 4;
       return;
       }

      for(;;) {
       // delprc is the stepsize to next breakpoint: decrease price of node
       // by delprc and compute the stepsize to the next breakpoint

       CNumber nxtbrk = Inf<CNumber>();

       Index arc = *tFOu;
       while( arc ) {  // look at all arcs out of the current node
        CNumber trc = mvflw2( arc , tDfct , tDDPos , Endn , X , U );

        incrsRC( arc , trc , delprc , nxtbrk , tDDNeg , DDPos , Endn );

        arc = NxtOu[ arc ];
        }

       arc = *tFIn;
       while( arc ) {  // look at all arcs into the current node
        CNumber trc = mvflw2( arc , tDfct , tDDPos , Startn , U , X );

        decrsRC( arc , trc , delprc , nxtbrk , tDDNeg , DDPos, Startn );

        arc = NxtIn[ arc ];
        }

       // if price of current node can be decreased further without
       // decreasing the dual cost (even if dual cost doesn't increase),
       // do another iteration to decrease the price

       if( LEZ( *tDDNeg , EpsDfct ) && ( nxtbrk < Inf<CNumber>() ) )
        delprc = nxtbrk;
       else
        break;
       }  // end for( ever )
      }
    }
   }
  }  // end( for( npass ) )
 }  // end( init_standard )

/*--------------------------------------------------------------------------*/

inline MCFClass::FNumber RelaxIV::svblncdarcs( cIndex node ,
                                     cIndex_Set tfst1 , cIndex_Set tnxt1 ,
                                     cIndex_Set tfst2 , cIndex_Set tnxt2 )
{
 FNumber delx = 0;
 nb_pos = 0;
 nb_neg = m;

 // check first (probably) balanced arcs incident to node: first arcs are
 // outgoing arcs in case of positive deficit node, incoming arcs otherwise

 Index arc = tfst1[ node ];

 while( arc ) {
  if( ETZ( RC[ arc ] , EpsCst ) && GTZ( X[ arc ] , EpsFlw ) )
   delx += X[ save[ nb_pos++ ] = arc ];

  arc = tnxt1[ arc ];
  }

 // check second (probably) balanced arcs incident to node: second arcs are
 // incoming arcs in case of positive deficit node, outgoing arcs otherwise

 arc = tfst2[ node ];

 while( arc ) {
  if( ETZ( RC[ arc ] , EpsCst ) && GTZ( U[ arc ] , EpsFlw ) )
   delx += U[ save[ --nb_neg ] = arc ];

  arc = tnxt2[ arc ];
  }

 return( delx );

 }  // end( svblncdarcs )

/*--------------------------------------------------------------------------*/

inline MCFClass::FNumber RelaxIV::dascnt( cIndex node , CNumber &delprc ,
					  cIndex_Set F1 , cIndex_Set Nxt1 ,
					  cIndex_Set F2 , cIndex_Set Nxt2 )
{
 nb_pos = 0;
 nb_neg = m;
 CNumber dp = delprc;
 delprc = Inf<CNumber>();
 FNumber delx = 0;
 Index arc = F1[ node ];

 while( arc ) {
  CNumber rdcost = ( RC[ arc ] += dp );

  if( ETZ( rdcost , EpsCst ) ) {
   save[ nb_pos++ ] = arc;
   delx += X[ arc ];
   }

  if( LTZ( rdcost , EpsCst ) && ( -rdcost < delprc ) )
   delprc = -rdcost;

  arc = Nxt1[ arc ];
  }

 for( arc = F2[ node ] ; arc ; ) {
  CNumber rdcost = ( RC[ arc ] -= dp );

  if( ETZ( rdcost , EpsCst ) ) {
   save[ --nb_neg ] = arc;
   delx += U[ arc ];
   }

  if( GTZ( rdcost , EpsCst ) && ( rdcost < delprc ) )
   delprc = rdcost;

  arc = Nxt2[ arc ];
  }

 return( delx );

 }  // end( dascnt )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::relist( cIndex node )
{
 // first clear the list of balanced arcs in/out of node

 Index arc = tfstou[ node ];
 tfstou[ node ] = 0;

 while( arc ) {
  Index nxtarc = tnxtou[ arc ];
  tnxtou[ arc ] = arc;
  arc = nxtarc;
  }

 arc = tfstin[ node ];
 tfstin[ node ] = 0;

 while( arc ) {
  Index nxtarc = tnxtin[ arc ];
  tnxtin[ arc ] = arc;
  arc = nxtarc;
  }

 // now add the currently balanced arcs to the list for this node and the
 // appropriate adjacent ones

 Index_Set t_save = save;
 Index j = nb_pos;   

 for( ; j-- ; ) {
  arc = *(t_save++);

  if( tnxtou[ arc ] == arc ) {
   tnxtou[ arc ] = tfstou[ Startn[ arc ] ];
   tfstou[ Startn[ arc ] ] = arc;
   }

  if( tnxtin[ arc ] == arc ) {
   tnxtin[ arc ] = tfstin[ Endn[ arc ] ];
   tfstin[ Endn[ arc ] ] = arc;
   }
  }

 for( j = m - nb_neg , t_save = save + m ; j-- ; ) {
  arc = *(--t_save);

  if( tnxtou[ arc ] == arc ) {
   tnxtou[ arc ] = tfstou[ Startn[ arc ] ];
   tfstou[ Startn[ arc ] ] = arc;
   }

  if( tnxtin[ arc ] == arc ) {
   tnxtin[ arc ] = tfstin[ Endn[ arc ] ];
   tfstin[ Endn[ arc ] ] = arc;
   }
  }
 }  // end( relist )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::AugFlow( cIndex augnod , cIndex root ,
                              cIndex node_p , cIndex node_m ,
                              cIndex_Set Term1 , cIndex_Set Term2 )

/* This subroutine execute the flow augmenting step: the augmenting path has
   been already found in the scanning step. The flow on the arcs of the path
   is modified: for the positive (negative) oriented arcs the flow is reduced
   (augmented) to reduce the total deficit. */
{
 Index ib = augnod;
 FNumber dx = - Dfct[ node_p ];

 while( ib != root ) {
  SIndex arc = Prdcsr[ ib ];

  if( arc > 0 ) {
   if( X[ arc ] < dx )
    dx = X[ arc ];

   ib = Term1[ arc ];
   }
  else {
   if( U[ arc = - arc ] < dx )
    dx = U[ arc ];

   ib = Term2[ arc ];
   }
  } 

 if( Dfct[ node_m ] < dx )
  dx = Dfct[ node_m ];

 if( GTZ( dx , EpsDfct ) ) {
  // increase (decrease) the flow of all forward (backward) arcs in the- - -
  // flow augmenting path, adjusting node deficits accordingly - - - - - - -

  if( ! queue[ augnod ] ) {
   queue[ prvnde ] = augnod;
   queue[ augnod ] = root;
   prvnde = augnod;
   }

  Dfct[ node_p ] += dx;
  Dfct[ node_m ] -= dx;

  for( ib = augnod ; ib != root ; ) {
   SIndex arc = Prdcsr[ ib ];

   if( arc > 0 ) {
    X[ arc ] -= dx;
    U[ arc ] += dx;

    ib = Term1[ arc ];
    }
   else {
    X[ arc = - arc ] += dx;
    U[ arc ] -= dx;

    ib = Term2[ arc ];
    }
   }
  }
 }  // end( AugFlow )

/*--------------------------------------------------------------------------*/

inline bool RelaxIV::Ascnt( cFNumber sdm , FNumber delx , Index &nlabel ,
                            bool &Switch , Index &nscan , Index &curnode ,
                            cIndex_Set Term1 , cIndex_Set Term2 ,
                            cIndex_Set F1 , cIndex_Set Nxt1 ,
                            cIndex_Set F2 , cIndex_Set Nxt2 )

/* This soubroutine performs the multinode price adjustment step. If the
   scanned nodes have positive deficit, it first checks if decreasing the
   price of the scanned nodes increases the dual cost. If yes, then it
   decreases the price of all scanned nodes.
   There are two possibilities for price decrease:

   - if Switch == true, then the set of scanned nodes corresponds to an
     elementary direction of maximal rate of ascent, hence the price of all
     scanned nodes is decreased until the next breakpoint in the dual cost is
     encountered: at this point, some arc becomes balanced, some nodes are
     added to the labeled set and the subroutine is exited;

   - if Switch == false, then the price of all scanned nodes are decreased
     until the rate of ascent becomes negative (this corresponds to the
     price adjustment step in which both the line search and the degenerate
     ascent iteration are implemented). */
{
 // Store the arcs between the set of scanned nodes and its complement and- -
 // compute delprc, the stepsize to the next breakpoint in the dual cost in -
 // the direction of decreasing prices of the scanned nodes. The arcs are - -
 // stored into save (the positive ones at the beginning and the negative - -
 // ones at the end) by looking at the arcs incident to either the set of - -
 // scanned nodes or its complement, whichever is smaller - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 pp( num_ascnt );

 FNumber dlx = 0;
 Index nsave_pos = 0;
 Index nsave_neg = m;
 CNumber delprc = Inf<CNumber>();

 if( nscan <= n / 2 ) {  // - - - - - - - - - - - - - - - - - - - - - - - - -
  Index_Set t_label = label;
  for( Index i = nscan ; i-- ; ) {
   Index node = *(t_label++);
   Index arc = F1[ node ];

   while( arc ) {  // arc points from a scanned node to an unscanned node
    Index node2 = Term2[ arc ];

    if( ! scan[ node2 ] ) {
     save[ nsave_pos++ ] = arc;

     CNumber rdcost = -RC[ arc ];
     if( ETZ( rdcost , EpsCst ) && ( Prdcsr[ node2 ] != SIndex( arc ) ) )
      dlx += X[ arc ];

     if( GTZ( rdcost , EpsCst ) && ( rdcost < delprc ) )
      delprc = rdcost;
     }

    arc = Nxt1[ arc ];
    }

   arc = F2[ node ];
   while( arc ) {  // arc points from an unscanned node to a scanned node
    Index node2 = Term1[ arc ];

    if( ! scan[ node2 ] ) {
     save[ --nsave_neg ] = arc;

     CNumber rdcost = RC[ arc ];
     if( ETZ( rdcost , EpsCst ) && ( Prdcsr[ node2 ] != - SIndex( arc ) ) )
      dlx += U[ arc ];

     if( GTZ( rdcost , EpsCst ) && ( rdcost < delprc ) )
      delprc = rdcost;
     }

    arc = Nxt2[ arc ];
    }
   }
  }
 else {  // nscan > n / 2 - - - - - - - - - - - - - - - - - - - - - - - - - -
  Bool_Vec tscan = scan;
  for( Index i = 0 ; i++ < n ; )
   if( ! *(++tscan) ) {
    SIndex Prdi = Prdcsr[ i ];
    Index arc = F2[ i ];

    while( arc ) {
     Index node2 = Term1[ arc ];

     if( scan[ node2 ] ) {
      save[ nsave_pos++ ] = arc;

      CNumber rdcost = -RC[ arc ];
      if( ETZ( rdcost , EpsCst ) && ( Prdi != RelaxIV::SIndex( arc ) ) )
       dlx += X[ arc ];

      if( GTZ( rdcost , EpsCst ) && ( rdcost < delprc ) )
       delprc = rdcost;
      }

     arc = Nxt2[ arc ];
     }

    for( arc = F1[ i ] ; arc ; ) {
     Index node2 = Term2[ arc ];

     if( scan[ node2 ] ) {
      save[ --nsave_neg ] = arc;

      CNumber rdcost = RC[ arc ];
      if( ETZ( rdcost , EpsCst ) && ( Prdi != -RelaxIV::SIndex( arc ) ) )
       dlx += U[ arc ];

      if( GTZ( rdcost , EpsCst ) && ( rdcost < delprc ) )
       delprc = rdcost;
      }

     arc = Nxt1[ arc ];
     }
    }

  }  // end else( nscan > n / 2 ) - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! Switch ) && ( GEZ( delx + dlx - sdm , EpsDfct ) ) ) {
  // check if the set of scanned nodes truly corresponds to a dual ascent
  // direction (here delx is the exact sum of the flow on arcs from the
  // scanned set to the unscanned set plus the residual capacity on arcs
  // from the unscanned set to scanned set): if this is *not* the case,
  // set Switch to true and exit subroutine

  Switch = true;
  return( true );
  }

 delx += dlx;

 // check that the problem is feasible- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  if( delprc == Inf<CNumber>() )  // can increase the dual cost indefinitely
   return( false );               // the primal problem is infeasible

  if( Switch ) {  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // decrease (increase when deficit is negative) the prices of the scanned
   // nodes, add more nodes to the labeled set and check if a newly labeled
   // node has negative deficit

   Index i = nsave_pos;
   Index_Set t_save = save;
   for( ; i-- ; ) {
    Index arc = *(t_save++);
    RC[ arc ] += delprc;

    if( ETZ( RC[ arc ] , EpsCst ) ) {
     Index node2 = Term2[ arc ];

     if( tnxtou[ arc ] == arc ) {
      tnxtou[ arc ] = tfstou[ Startn[ arc ] ];
      tfstou[ Startn[ arc ] ] = arc;
      }

     if( tnxtin[ arc ] == arc ) {
      tnxtin[ arc ] = tfstin[ Endn[ arc ] ];
      tfstin[ Endn[ arc ] ] = arc;
      }

     if( ! mark[ node2 ] ) {
      Prdcsr[ node2 ] = arc;
      label[ nlabel++ ] = node2;
      mark[ node2 ] = true;
      }
     }
    }

   for( i = m - nsave_neg , t_save = save + m ; i-- ; ) {
    Index arc = *(--t_save);
    RC[ arc ] -= delprc;

    if( ETZ( RC[ arc ] , EpsCst ) ) {
     Index node2 = Term1[ arc ];

     if( tnxtou[ arc ] == arc ) {
      tnxtou[ arc ] = tfstou[ Startn[ arc ] ];
      tfstou[ Startn[ arc ] ] = arc;
      }

     if( tnxtin[ arc ] == arc ) {
      tnxtin[ arc ] = tfstin[ Endn[ arc ] ];
      tfstin[ Endn[ arc ] ] = arc;
      }

     if( ! mark[ node2 ] ) {
      Prdcsr[ node2 ] = -arc;
      label[ nlabel++ ] = node2;
      mark[ node2 ] = true;
      }
     }
    }

   break;  // force termination
   }

  // decrease the prices of the scanned nodes by delprc and adjust the flow
  // to maintain complementary slackness with the prices

  Index h = nsave_pos;
  Index_Set tsv = save;
  for( ; h-- ; ) {
   Index arc = *(tsv++);
   if( ETZ( RC[ arc ] , EpsCst ) ) {
    FNumber t2 = X[ arc ];
    Index t3 = Startn[ arc ];
    Dfct[ t3 ] -= t2;

    if( ! queue[ t3 ] ) {
     queue[ prvnde ] = t3;
     queue[ t3 ] = curnode;
     prvnde = t3;
     }

    Dfct[ t3 = Endn[ arc ] ] += t2;

    if( ! queue[ t3 ] ) {
     queue[ prvnde ] = t3;
     queue[ t3 ] = curnode;
     prvnde = t3;
     }

    U[ arc ] += t2;
    X[ arc ] = 0;
    }

   RC[ arc ] += delprc;
   if( ETZ( RC[ arc ] , EpsCst ) ) {
    delx += X[ arc ];

    if( tnxtin[ arc ] == arc ) {
     Index j = Endn[ arc ];
     tnxtin[ arc ] = tfstin[ j ];
     tfstin[ j ] = arc;
     }

    if( tnxtou[ arc ] == arc ) {
     Index j = Startn[ arc ];
     tnxtou[ arc ] = tfstou[ j ];
     tfstou[ j ] = arc;
     }
    }
   }

  for( h = m - nsave_neg , tsv = save + m ; h-- ; ) {
   Index arc = *(--tsv);
   if( ETZ( RC[ arc ] , EpsCst ) ) {
    FNumber t2 = U[ arc ];
    Index t3 = Startn[ arc ];
    Dfct[ t3 ] += t2;

    if( ! queue[ t3 ] ) {
     queue[ prvnde ] = t3;
     queue[ t3 ] = curnode;
     prvnde = t3;
     }

    Dfct[ t3 = Endn[ arc ] ] -= t2;

    if( ! queue[ t3 ] ) {
     queue[ prvnde ] = t3;
     queue[ t3 ] = curnode;
     prvnde = t3;
     }

    X[ arc ] += t2;
    U[ arc ] = 0;
    }

   RC[ arc ] -= delprc;
   if( ETZ( RC[ arc ] , EpsCst ) ) {
    delx += U[ arc ];

    if( tnxtin[ arc ] == arc ) {
     Index j = Endn[ arc ];
     tnxtin[ arc ] = tfstin[ j ];
     tfstin[ j ] = arc;
     }

    if( tnxtou[ arc ] == arc ) {
     Index j = Startn[ arc ];
     tnxtou[ arc ] = tfstou[ j ];
     tfstou[ j ] = arc;
     }
    }
   }

  if( GEZ( sdm - delx , EpsDfct ) ) {
   // the set of scanned nodes still correspond to a dual (possibly
   // degenerate) ascent direction: compute the stepsize delprc to the
   // next breakpoint in the dual cost function

   Index i = nsave_pos;
   Index_Set t_save = save;
   for( delprc = Inf<CNumber>() ; i-- ; ) {
    CNumber rdcost = -RC[ *(t_save++) ];

    if( GTZ( rdcost , EpsCst ) && ( rdcost < delprc ) )
     delprc = rdcost;
    }

   for( i = m - nsave_neg , t_save = save + m ; i-- ; ) {
    CNumber rdcost = RC[ *(--t_save) ];

    if( GTZ( rdcost , EpsCst ) && ( rdcost < delprc ) )
     delprc = rdcost;
    }
   }
  } while( GEZ( sdm - delx , EpsDfct ) &&
           ( ( delprc < Inf<CNumber>() ) || GTZ( sdm - delx , EpsDfct ) ) );

 return( true );

 }  // end( Ascnt )

/*--------------------------------------------------------------------------*/

inline MCFClass::CNumber RelaxIV::nxtbrkpt( cIndex_Set tSt1 , cIndex_Set NSt1 ,
					    cIndex_Set tSt2 , cIndex_Set NSt2 )
{
 CNumber delprc = Inf<CNumber>();
 Index arc = *tSt1;
 while( arc ) {
  CNumber trc = RC[ arc ];
  if( GTZ( trc , EpsCst ) && ( trc < delprc ) )
   delprc = trc;

  arc = NSt1[ arc ];
  }

 for( arc = *tSt2 ; arc ; ) {
  CNumber trc = -RC[ arc ];
  if( GTZ( trc , EpsCst ) && ( trc < delprc ) )
   delprc = trc;

  arc = NSt2[ arc ];
  }

 return( delprc );

 }  // end( nxtbrkpt )

/*--------------------------------------------------------------------------*/

inline MCFClass::CNumber RelaxIV::mvflw1( cIndex arc , FRow tDfct ,
					  FRow tDDNeg , cIndex_Set Term ,
					  FRow Flow1 , FRow Flow2 )
{
 CNumber trc = RC[ arc ];

 if( ETZ( trc , EpsCst ) ) {
  Index t1 = Term[ arc ];
  FNumber f = Flow1[ arc ];

  if( GTZ( f , EpsFlw ) ) {
   (*tDfct) += f;
   Dfct[ t1 ] -= f;
   Flow2[ arc ] = f;
   Flow1[ arc ] = 0;
   }
  else
   f = Flow2[ arc ];

  (*tDDNeg) -= f;
  DDPos[ t1 ] -= f;
  }

 return( trc );

 }  // end( mvflw1 )

/*--------------------------------------------------------------------------*/

inline MCFClass::CNumber RelaxIV::mvflw2( cIndex arc , FRow tDfct ,
					  FRow tDDPos , cIndex_Set Term ,
					  FRow Flow1 , FRow Flow2 )
{
 CNumber trc = RC[ arc ];

 if( ETZ( trc , EpsCst ) ) {
  Index t1 = Term[ arc ];
  FNumber f = Flow1[ arc ];

  if( GTZ( f , EpsFlw ) ) {
   (*tDfct) -= f;
   Dfct[ t1 ] += f;
   Flow2[ arc ] = f;
   Flow1[ arc ] = 0;
   }
  else
   f = Flow2[ arc ];

  (*tDDPos) -= f;
  DDNeg[ t1 ] -= f;
  }

 return( trc );

 }  // end( mvflw2 )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::decrsRC( cIndex arc , CNumber trc , cCNumber delprc ,
                              CNumber &nxtbrk , FRow tDD1 , FRow DD2 ,
                              cIndex_Set Term )
{
 trc -= delprc;

 if( GTZ( trc , EpsCst ) && ( trc < nxtbrk ) )
  nxtbrk = trc;
 else
  if( ETZ( trc , EpsCst ) ) {  // arc goes from inactive to balanced: update
   *tDD1 += U[ arc ];          // the rate of dual ascent at current node and
   DD2[ Term[ arc ] ] += U[ arc ];  // at its neighbours
   }

 RC[ arc ] = trc;

 }  // end( dcrsRC )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::incrsRC( cIndex arc , CNumber trc , cCNumber delprc ,
                              CNumber &nxtbrk , FRow tDD1 , FRow DD2 ,
                              cIndex_Set Term )
{ 
 trc += delprc;

 if( LTZ( trc , EpsCst ) && ( -trc < nxtbrk ) )
  nxtbrk = -trc;
 else
  if( ETZ( trc , EpsCst ) ) {  // arc goes from active to balanced: update
   (*tDD1) += X[ arc ];        // the rate of dual ascent at current node and
   DD2[ Term[ arc ] ] += X[ arc ];   // at its neighbours
   }

 RC[ arc ] = trc;

 }  // end( incrsRC )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::chgcsti( cIndex i , CNumber NCost )
{
 CNumber RCi = RC[ i ];
 cCNumber DCst = NCost - C[ i ];
 C[ i ] = NCost;

 #if( ( ! SAME_GRPH_RIV ) || DYNMC_MCF_RIV )
  if( RCi < Inf<CNumber>() )
 #endif
  {
   RC[ i ] = ( RCi += DCst );

   if( GTZ( RCi , EpsCst ) ) {  // RC[ i ] > 0 - - - - - - - - - - - - - - -
    cFNumber Xi = X[ i ];
    if( GTZ( Xi , EpsFlw ) ) {
     X[ i ] = 0;
 
     Dfct[ Startn[ i ] ] -= Xi;
     Dfct[ Endn[ i ] ] += Xi;
     U[ i ] += Xi;
     }
    }
   else
    if( LTZ( RCi , EpsCst ) ) {  // RC[ i ] < 0- - - - - - - - - - - - - - -
     cFNumber Ui = U[ i ];
     if( GTZ( Ui , EpsFlw ) ) {
      U[ i ] = 0;

      Dfct[ Startn[ i ] ] += Ui;
      Dfct[ Endn[ i ] ] -= Ui;
      X[ i ] += Ui;
      }
     }
    else {  // RC[ i ] == 0 - - - - - - - - - - - - - - - - - - - - - - - - -
     if( tnxtou[ i ] == i ) {
      Index node = Startn[ i ];
      tnxtou[ i ] = tfstou[ node ];
      tfstou[ node ] = i;
      }

     if( tnxtin[ i ] == i ) {
      Index node = Endn[ i ];
      tnxtin[ i ] = tfstin[ node ];
      tfstin[ node ] = i;
      }
     }
   }  // end( if( Cap[ i ] ) )
 }  // end( chgcsti )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::chgcapi( cIndex i , cFNumber NCap )
{
 Cap[ i ] = NCap;
 cFNumber diffX = NCap - X[ i ];

 if( LTZ( RC[ i ] , EpsCst ) || LTZ( diffX , EpsFlw ) ) {
  Dfct[ Startn[ i ] ] += diffX;
  Dfct[ Endn[ i ] ] -= diffX;
  X[ i ] = NCap;
  U[ i ] = 0;
  }
 else
  U[ i ] = diffX;

 }  // end( chgcapi )

/*--------------------------------------------------------------------------*/

#if( DYNMC_MCF_RIV )

inline void RelaxIV::delarci( cIndex arc )
{
 Index arc1 = FOu[ Startn[ arc ] ];
 if( arc1 == arc )
  FOu[ Startn[ arc ] ] = NxtOu[ arc1 ];
 else {
  Index arc2;
  do {
   arc2 = arc1;
   arc1 = NxtOu[ arc1 ];
   } while( arc1 != arc );

  NxtOu[ arc2 ] = NxtOu[ arc ];
  }

 if( ( arc1 = FIn[ Endn[ arc ] ] ) == arc )
  FIn[ Endn[ arc ] ] = NxtIn[ arc1 ];
 else {
  Index arc2;
  do {
   arc2 = arc1;
   arc1 = NxtIn[ arc1 ];
   } while( arc1 != arc );

  NxtIn[ arc2 ] = NxtIn[ arc ];
  }

 if( status || ( ! Senstv ) )
  status = MCFClass::kUnSolved;
 else {
  Dfct[ Startn[ arc ] ] -= X[ arc ];
  Dfct[ Endn[ arc ] ] += X[ arc ];

  if( tnxtou[ arc ] != arc ) {
   if( ( arc1 = tfstou[ Startn[ arc ] ] ) == arc )
    tfstou[ Startn[ arc ] ] = tnxtou[ arc ];
   else {
    Index arc2;
    do {
     arc2 = arc1;
     arc1 = tnxtou[ arc1 ];
     } while( arc1 != arc );

    tnxtou[ arc2 ] = tnxtou[ arc ];
    }

   tnxtou[ arc ] = arc;
   }

  if( tnxtin[ arc ] != arc ) {
   if( ( arc1 = tfstin[ Endn[ arc ] ] ) == arc )
    tfstin[ Endn[ arc ] ] = tnxtin[ arc ];
   else {
    Index arc2;
    do {
     arc2 = arc1;
     arc1 = tnxtin[ arc1 ];
     } while( arc1 != arc );

    tnxtin[ arc2 ] = tnxtin[ arc ];
    }

   tnxtin[ arc ] = arc;
   }
  }

 X[ arc ] = 0;
 RC[ arc ] = Inf<CNumber>();

 }  // end( delarci )

/*--------------------------------------------------------------------------*/

#if( DYNMC_MCF_RIV > 1 )

inline void RelaxIV::addarci( cIndex arc )
{
 RC[ arc ] = C[ arc ];
 cIndex sn = Startn[ arc ]; 
 cIndex en = Endn[ arc ];

 if( status || ( ! Senstv ) )
  status = MCFClass::kUnSolved;
 else {
  if( PiOwnr != this )  // compute the dual prices
   cmptprices();

  cFNumber Ua = Cap[ arc ];
  cCNumber RCa = ( RC[ arc ] += ( Pi[ sn ] - Pi[ en ] ) );

  if( LTZ( RCa , EpsCst ) ) {
   Dfct[ sn ] += Ua;
   Dfct[ en ] -= Ua;
   X[ arc ] = Ua;
   U[ arc ] = 0;
   }
  else {
   X[ arc ] = 0;
   U[ arc ] = Ua;
   }

  if( ETZ( RCa , EpsCst ) ) {
   tnxtou[ arc ] = tfstou[ sn ];
   tfstou[ sn ] = arc;
   tnxtin[ arc ] = tfstin[ en ];
   tfstin[ en ] = arc;
   }
  }

 // update FS & BS: this *must* be done *after* the call to cmptprices()- - -

 NxtOu[ arc ] = FOu[ sn ];
 FOu[ sn ] = arc;
 NxtIn[ arc ] = FIn[ en ];
 FIn[ en ] = arc;

 }  // end( addarci )

#endif  // DYNMC_MCF_RIV > 1
#endif  // DYNMC_MCF_RIV

/*--------------------------------------------------------------------------*/

inline void RelaxIV::cmptprices( void )
{
 CRow tPi = Pi + n;
 for( ; tPi > Pi ; )
  *(tPi--) = Inf<CNumber>();        // reset all potentials to +INF

 Index pcnt = n;           // # of potentials still to compute
 do {                      // outer loop: for all connected components- - - -
  do
   ++tPi;
  while( *tPi < Inf<CNumber>() );   // search "root" of unvisited component

  *tPi = 0;                // set any initial potential
  if( ! --pcnt )
   break;

  CNumber Pstart = 0;
  Index_Set tq = queue;
  for( Index start = tPi - Pi ;; ) {  // inner loop: visit this
   Index arc = FOu[ start ];          // connected component
   while( arc ) {                     // scan FS( start )
    #if( ( ! SAME_GRPH_RIV ) || DYNMC_MCF_RIV )
     if( RC[ arc ] < Inf<CNumber>() )
    #endif
     {
      cIndex end = Endn[ arc ];
      if( Pi[ end ] == Inf<CNumber>() ) {
       Pi[ *(++tq) = end ] = C[ arc ] - RC[ arc ] + Pstart;
       pcnt--;
       }
      }

    arc = NxtOu[ arc ];
    }

   if( ! pcnt )
    break;

   for( arc = FIn[ start ] ; arc ; ) {  // scan BS( start )
    #if( ( ! SAME_GRPH_RIV ) || DYNMC_MCF_RIV )
     if( RC[ arc ] < Inf<CNumber>() )
    #endif
     {
      cIndex end = Startn[ arc ];
      if( Pi[ end ] == Inf<CNumber>() ) {
       Pi[ *(++tq) = end ] = RC[ arc ] + Pstart - C[ arc ];
       pcnt--;
       }
      }

    arc = NxtIn[ arc ];
    }

   if( pcnt && ( tq > queue ) )
    Pstart = Pi[ start = *(tq--) ];
   else
    break;
   }                // end inner loop
  } while( pcnt );  // end outer loop - - - - - - - - - - - - - - - - - - - -

 #if( P_ALLOC )
  PiOwnr = this;  // "sign" these potentials
 #endif

 }  // end( cmptprices )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( AUCTION )

void RelaxIV::auction( void )
{
 // this method uses a version of the Auction/Shortest Paths algorithm for
 // Min Cost Flow problems to compute initial flow and prices for the
 // Relaxation algorithm

 // reset B - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FRow tDfct = Dfct + n;
 for( cFRow tB = B + n ; tB > B ; )
  *(tDfct--) = *(tB--);

 // compute maxcost and mincost - - - - - - - - - - - - - - - - - - - - - - -

 CNumber maxcost = - C_LARGE;
 CNumber mincost = C_LARGE;
 {
  // meanwhile, initialize the arc flows X and the reduced capacities U to
  // satisfy the complementary slackness with the reduced costs RC

  FRow tX = X + m;
  FRow tU = U + m;
  cCRow tC = C + m;
  CRow tRC = RC + m;
  cFRow tCap = Cap + m;
  for( ; tC > C ; tX-- , tU-- , tC-- , tRC-- , tCap-- ) {
   #if( ( ! SAME_GRPH_RIV ) || DYNMC_MCF_RIV )
    if( *tRC == Inf<CNumber>() )
     continue;
   #endif

   cCNumber RCi = *tRC = *tC;
   if( maxcost < RCi )
    maxcost = RCi;
   if( mincost > RCi )
    mincost = RCi;

   cFNumber f = *(tCap--);
   if( LTZ( RCi , EpsCst ) ) {
    *tU = 0;
    *tX = f;
    cIndex i = tX - X;
    Dfct[ Startn[ i ] ] += f;
    Dfct[ Endn[ i ] ] -= f;
    }
   else {
    *tU = f;
    *tX = 0;
    }
   }
  }

 // set initial eps - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CNumber eps = int( ( maxcost - mincost ) / maxdf );
 if( eps < 1 )
  eps = 1;

 // set initial prices to zero- - - - - - - - - - - - - - - - - - - - - - - -
  
 for( CRow tP = Pi + n ; tP > Pi ; )
  *(tP--) = 0;

 // start of the first scaling phase- - - - - - - - - - - - - - - - - - - - -

 #if( RELAXIV_STATISTICS )
  nsp = 0;
 #endif
 int pass = 0;
 FNumber thresh_dfct = 0;
 bool cond100;

 do {   /*100*/
  cond100 = false;
  if( ( ++pass == npassauct ) || ( eps == 1 ) )
   crash = false;

  Index nolist = 0;

  /*------------ Construct list of  positive surplus nodes ------------*/
  /*--------------- and queue of negative surplus nodes ---------------*/
  {
   SIndex_Set tPrdcsr = Prdcsr;
   SIndex_Set textend_arc = extend_arc;
   Bool_Vec tmark = mark;
   Index_Set tqueue = queue;
   CRow tSB_level = SB_level;
   FRow tDfct = Dfct;
   for( Index node = 0 ; node++ < n ; ) {
    *(++tPrdcsr) = 0;
    *(++tmark) = false;
    *(++textend_arc) = 0;
    *(++tSB_level) = -C_LARGE;
    *(++tqueue) = node + 1;
    if( GTZ( *(++tDfct) , EpsDfct ) )
     save[ nolist++ ] = node;
    }
   }

  Index root = queue[ n ] = 1;
  prvnde = lastq = n;

  /*- Initialization with down iterations  for negative surplus nodes -*/
  {
   Index i = nolist;
   Index_Set t_save = save;
   for( ; i-- ; ) {
    Index node = *(t_save++);
    pp( nsp );

    /*-------- Build the list arcs w/ room for pushing flow -------*/
    /*--------- and find proper prices  for down iteration ---------*/

    CNumber bstlevel = -C_LARGE;
    SIndex extarc;
    FpushF[ node ] = 0;
    Index arc = FOu[ node ];
    Index last;
    while( arc ) {
     if( GTZ( U[ arc ] , EpsFlw ) )
      if( ! FpushF[ node ] ) {
       FpushF[ node ] = arc;
       NxtpushF[ arc ] = 0;
       last = arc;
       }
      else {
       NxtpushF[ last ] = arc;
       NxtpushF[ arc ] = 0;
       last = arc;
       }

     if( GTZ( X[ arc ] , EpsFlw ) ) {
      CNumber new_level = Pi[ Endn[ arc ] ] + RC[ arc ];
      if( new_level > bstlevel ) {
       bstlevel = new_level;
       extarc = arc;
       }
      }

     arc = NxtOu[ arc ];
     }

    FpushB[ node ] = 0;
    arc = FIn[ node ];
    while( arc ) {
     if( GTZ( X[ arc ] , EpsFlw ) )
      if( ! FpushB[ node ] ) {
       FpushB[ node ] = arc;
       NxtpushB[ arc ] = 0;
       last = arc;
       }
      else {
       NxtpushB[ last ] = arc;
       NxtpushB[ arc ] = 0;
       last = arc;
       }

     if( GTZ( U[ arc ] , EpsFlw ) ) {
      CNumber new_level = Pi[ Startn[ arc ] ] - RC[ arc ];
      if( new_level > bstlevel ) {
       bstlevel = new_level;
       extarc = -arc;
       }
      }

     arc = NxtIn[ arc ];
     }

    extend_arc[ node ] = extarc;    
    Pi[ node ] = bstlevel - eps;
    }
   }

  /*-------- Start the augmention cycles of the scaling phase. --------*/

  bool cond200;
  do {
   cond200 = false;

   if( GTZ( thresh_dfct - Dfct[ root ] , EpsDfct ) ) {
    Index term = root;
    mark[ root ] = true;

    /*-------- Main forward algorithm with root as origin. --------*/

    bool cond500;
    do {
     cond500 = false;
     int salto = 550;
     CNumber bstlevel , seclevel;

     /*---------- Start of a new forward iteration ------------*/

     CNumber pterm = Pi[term];
     SIndex extarc = extend_arc[ term ];

     if( ! extarc ) {
      /*- Build the list of arcs w/ room for pushing flow --*/

      FpushF[ term ] = 0;
      Index arc = FOu[ term ];
      Index last;
      while( arc ) {
       if( GTZ( U[ arc ] , EpsFlw ) )
        if( ! FpushF[ term ] ) {
         FpushF[ term ] = arc;
         NxtpushF[ arc ] = 0;
         last = arc;
         }
        else {
         NxtpushF[ last ] = arc;
         NxtpushF[ arc ] = 0;
         last = arc;
         }

       arc = NxtOu[ arc ];
       }

      FpushB[ term ] = 0;
      arc = FIn[ term ];
      while( arc ) {
       if( GTZ( X[ arc ] , EpsFlw ) )
        if( ! FpushB[ term ] ) {
         FpushB[ term ] = arc;
         NxtpushB[ arc ] = 0;
         last = arc;
         }
        else {
         NxtpushB[ last ] = arc;
         NxtpushB[ arc ] = 0;
         last = arc;
         }

       arc = NxtIn[ arc ];
       }

      salto = 600;
      }
     else {
      /*------- Speculative path  extension attempt -------------------------
        Note: extarc > 0 means that extarc is oriented from the root to the
        destinations, extarc < 0 means that extarc is oriented from the
        destinations to the root, extarc = 0 or Prdarc = 0 means the
        extension arc or the predecessor arc, respectively, has not been
        established. */

      if( extarc > 0 ) {
       if( ETZ( U[ extarc ] , EpsFlw ) ) {
        seclevel = SB_level[ term ];
        salto = 580;
        }
       else {
        Index end = Endn[ extarc ];
        bstlevel = Pi[ end ] + RC[ extarc ];
        if( GEZ( pterm - bstlevel , EpsCst ) ) {
         if( mark[ end ] )
          salto = 1200;
         else {
          term = end;
          Prdcsr[ term ] = extarc;
          mark[ term ] = true;

          // if negative surplus node is found, do an augmentation

          if( GTZ( Dfct[ term ] , EpsDfct ) )
           salto = 2000;
          else
           cond500 = true;  // return for another iteration
          }
         }
        }
       }
      else {
       extarc = -extarc;
       if( ETZ( X[ extarc ] , EpsFlw ) ) {
        seclevel = SB_level[ term ];
        salto = 580;
        }
       else {
        Index start = Startn[ extarc ];
        bstlevel = Pi[ start ] - RC[ extarc ];
        if( GEZ( pterm - bstlevel , EpsCst ) ) {
         if( mark[ start ] )
          salto = 1200;
         else {
          term = start;
          Prdcsr[ term ] = -extarc;
          mark[ term ] = true;

          // if negative surplus node is found, do an augmentation

          if( GTZ( Dfct[ term ] , EpsDfct ) )
           salto = 2000;
          else
           cond500 = true;  // return for another iteration
          }
         }
        }
       }
      }

     if( cond500 )
      continue;

     do {
      // second best logic test applied to save a full node scan: if old
      // best level continues to be best go for another contraction

      if( salto == 550 ) {
       seclevel = SB_level[ term ];
       if( LEZ( bstlevel - seclevel , EpsCst ) )
        salto = 800;
       else
        salto = 580;
       }

      // if second best can be used, either do a contraction or start over
      // with a speculative extension

      if( salto == 580 )
       if( seclevel > -C_LARGE ) {
        extarc = SB_arc[ term ];
        if( extarc > 0 )
         if( ETZ( U[ extarc ] , EpsFlw ) )
          salto = 600;
         else
          bstlevel = Pi[ Endn[ extarc ] ] + RC[ extarc ];
        else
         if( ETZ( X[ -extarc ] , EpsFlw ) )
          salto = 600;
         else
          bstlevel = Pi[ Startn[ -extarc ] ] - RC[ -extarc ];

        if( salto != 600 )
         if( ETZ( bstlevel - seclevel , EpsCst ) ) {
          SB_level[ term ] = -C_LARGE;
          extend_arc[ term ] = extarc;
          salto = 800;
          }
         else
          salto = 600;
        }
       else
        salto = 600;

      // extention/contraction attempt was unsuccessful, so scan terminal
      // node

      if( salto == 600 ) {
       pp( nsp );
       SIndex secarc;
       bstlevel = seclevel = C_LARGE;

       Index arc = FpushF[ term ];
       while( arc ) {
        CNumber new_level = Pi[ Endn[ arc ] ] + RC[ arc ];
        if( GTZ( seclevel-new_level , EpsCst ) )
         if( GTZ( bstlevel - new_level , EpsCst ) ) {
          seclevel = bstlevel;
          bstlevel = new_level;
          secarc = extarc;
          extarc = arc;
          }
         else {
          seclevel = new_level;
          secarc = arc;
          }

        arc = NxtpushF[ arc ];
        }

       for( arc = FpushB[ term ] ; arc ; ) {
        CNumber new_level = Pi[ Startn[ arc ] ] - RC[ arc ];
        if( GTZ( seclevel - new_level , EpsCst ) )
         if( GTZ( bstlevel - new_level , EpsCst ) ) {
          seclevel = bstlevel;
          bstlevel = new_level;
          secarc = extarc;
          extarc = -arc;
          }
         else {
          seclevel = new_level;
          secarc = -arc;
          }

        arc = NxtpushB[ arc ];
        }

       SB_level[ term ] = seclevel;
       SB_arc[ term ] = secarc;
       extend_arc[ term ] = extarc;
       salto = 800;
       }

      // end of node scan: if the terminal node is the root, adjust its
      // price and change root

      if( salto == 800 ) {
       if( term == root ) {
        Pi[ term ] = bstlevel + eps;
        if( pterm >= C_LARGE ) {
         status = MCFClass::kUnfeasible;  // the problem is unfeasible
         error_node = root;
         error_info = 8;
         return;
         }

        mark[ root ] = false;
        prvnde = root;
        root = queue[ root ];
        cond200 = true;
        break;
        }

       // check whether extension or contraction

       SIndex prd = Prdcsr[ term ];
       CNumber prevlevel; 
       Index pr_term;   
       if( prd > 0 ) {
        pr_term = Startn[ prd ];
        prevlevel = Pi[ pr_term ] - RC[ prd ];
        }
       else {
        pr_term = Endn[ -prd ];
        prevlevel = Pi[ pr_term ] + RC[ -prd ];
        }

       if( GTZ( prevlevel - bstlevel , EpsCst ) ) {  /* Path extension */
        if( GEZ( prevlevel - ( bstlevel + eps ) , EpsCst ) )
         Pi[ term ] = bstlevel + eps;
        else
         Pi[ term ] = prevlevel;

        Index nde = ( extarc > 0 ? Endn[ extarc ] : Startn[ - extarc ] );
        if( mark[ nde ] )
         salto = 1200;
        else
         term = nde;

        if( salto != 1200 ) {
         Prdcsr[ term ] = extarc;
         mark[ term ] = true;

         // if negative surplus node is found, do an augmentation

         if( GTZ( Dfct[ term ] , EpsDfct ) )
          salto = 2000;
         else {
          cond500 = true;
          break;
          }
         }
        }
       else {  /* Path contraction. */
        Pi[ term ] = bstlevel + eps;
        mark[ term ] = false;
        term = pr_term;
        if( pr_term != root )
         if( LEZ( bstlevel - ( pterm + eps ) , EpsCst ) )
          salto = 2000;

        if( salto != 2000 ) {
         pterm = Pi[ term ];
         extarc = prd;
         if( prd > 0 )
          bstlevel +=eps + RC[ prd ];
         else
          bstlevel += eps - RC[ -prd ];

         // do a second best test and, if that fails, do a full node scan

         salto = 550;
         }
        }
       }

      // a cycle is about to form; do a retreat sequence

      if( salto == 1200 ) {
       Index pr_term;
       Index node = term;
       while( node != root ) {
        mark[ node ] = false;
        SIndex prd = Prdcsr[ node ];
        if( prd > 0 ) {
         pr_term = Startn[ prd ];
         CNumber prdRC = Pi[ pr_term ] - Pi[ node ] - RC[ prd ] - eps;
         if( ETZ( prdRC , EpsCst ) )  // beware of ETZ()
          node = pr_term;
         else
          break;
         }
        else {
         pr_term = Endn[ -prd ];
         CNumber prdRC = Pi[ pr_term ] - Pi[ node ] - RC[ -prd ] - eps;
         if( ETZ( prdRC , EpsCst ) )  // beware of ETZ()
          node = pr_term;
         else
          break;
         }
        }

       // do a full scan and price rise at pr_term

       if( node != root ) {
        pp( nsp );
        SIndex secarc;
        bstlevel = seclevel = C_LARGE;
        Index arc = FpushF[ pr_term ];
        while( arc ) {
         CNumber new_level = Pi[ Endn[ arc ] ] + RC[ arc ];
         if( GTZ( seclevel - new_level , EpsCst ) )
          if( GTZ( bstlevel - new_level , EpsCst ) ) {
           seclevel = bstlevel;
           bstlevel = new_level;
           secarc = extarc;
           extarc = arc;
           }
          else {
           seclevel = new_level;
           secarc = arc;
           }

         arc = NxtpushF[ arc ];
         }

        for( arc = FpushB[ pr_term ] ; arc ; ) {
         CNumber new_level = Pi[ Startn[ arc ] ] - RC[ arc ];
         if( GTZ( seclevel - new_level , EpsCst ) )
          if( GTZ( bstlevel - new_level , EpsCst ) ) {
           seclevel = bstlevel;
           bstlevel = new_level;
           secarc = extarc;
           extarc = -arc;
           }
          else {
           seclevel = new_level;
           secarc = -arc;
           }

         arc = NxtpushB[ arc ];
         }

        SB_level[ pr_term ] = seclevel;
        SB_arc[ pr_term ] = secarc;
        extend_arc[ pr_term ] = extarc;
        Pi[ pr_term ] = bstlevel + eps;

        if( pr_term == root ) {
         mark[ prvnde = root ] = false;
         root = queue[ root ];
         cond200 = true;
         break;
         }

        mark[ pr_term ] = false;
        SIndex prd = Prdcsr[ pr_term ];
        if( prd > 0 )
         term = Startn[ prd ];
        else
         term = Endn[ -prd ];

        if( term == root ) {
         mark[ prvnde = root ] = false;
         root = queue[ root ];
         cond200 = true;
         break;
         }
        }

       salto = 2000;
       }

      /*------ End of auction/shortest path routine. ------*/
      /* Do augmentation from root and correct the push lists */

      if( salto == 2000 ) {
       FNumber incr = -Dfct[ root ];
       Index node = root;

       do {
        extarc = extend_arc[ node ];
        mark[ node ] = false;
        if( extarc > 0 ) {
         node = Endn[ extarc ];
         if( GTZ( incr - U[ extarc ] , EpsFlw ) )
          incr = U[ extarc ];
         }
        else {
         node = Startn[ -extarc ];
         if( GTZ( incr - X[ -extarc ] , EpsFlw ) )
          incr = X[ -extarc ];
         }
        } while( node != term );

       mark[ term ] = false;
       if( GTZ( Dfct[ term ] , EpsDfct ) )
        if( GTZ( incr - Dfct[ term ] , EpsDfct ) )
         incr = Dfct[ term ];

       node = root;
       do {
        extarc = extend_arc[ node ];
        if( extarc > 0 ) {
         Index end = Endn[ extarc ];

         /*----- Add arc to the reduced graph -----*/

         if( ETZ( X[ extarc ] , EpsFlw ) ) {
          NxtpushB[ extarc ] = FpushB[ end ];
          FpushB[ end ] = extarc;
          CNumber new_level = Pi[ node ] - RC[ extarc ];
          if( GTZ( SB_level[ end ] - new_level , EpsCst ) ) {
           SB_level[ end ] = new_level;
           SB_arc[ end ] = -extarc;
           }
          }

         X[ extarc ] += incr;
         U[ extarc ] -= incr;

         /*-- Remove arc from the reduced graph --*/

         if( ETZ( U[ extarc ] , EpsFlw ) ) {
          SIndex arc = FpushF[ node ];
          if( arc == extarc )
           FpushF[ node ] = NxtpushF[ arc ];
          else {
           SIndex prevarc = arc;
           for( arc = NxtpushF[ arc ] ; arc ; ) {
            if( arc == extarc ) {
             NxtpushF[ prevarc ] = NxtpushF[ arc ];
             break;
             }
            prevarc = arc;
            arc = NxtpushF[ arc ];
            }
           }
          }

         node = end;
         }
        else {
         extarc = -extarc;
         Index start = Startn[ extarc ];

         /*---- Add arc to the reduced graph -----*/

         if( ETZ( U[ extarc ] , EpsFlw ) ) {
          NxtpushF[ extarc ] = FpushF[ start ];
          FpushF[ start ] = extarc;
          CNumber new_level = Pi[ node ] + RC[ extarc ];
          if( GTZ( SB_level[ start ] - new_level , EpsCst ) ) {
           SB_level[ start ] = new_level;
           SB_arc[ start ] = extarc;
           }
          }

         U[ extarc ] += incr;
         X[ extarc ] -= incr;

         /*-- Remove arc from the reduced graph ---*/

         if( ETZ( X[ extarc ] , EpsFlw ) ) {
          SIndex arc = FpushB[ node ];
          if( arc == extarc )
           FpushB[ node ] = NxtpushB[ arc ];
          else {
           SIndex prevarc = arc;
           for( arc = NxtpushB[ arc ] ; arc ; ) {
            if( arc == extarc ) {
             NxtpushB[ prevarc ] = NxtpushB[ arc ];
             break;
             }
            prevarc = arc;
            arc = NxtpushB[ arc ];
            }
           }
          }

         node = start;
         }
        } while( node != term );

       Dfct[ term ] -= incr;
       Dfct[ root ] += incr;

       /*------- Insert term in the queue if it --------*/
       /*--------- has a large enough surplus ----------*/

       if( GTZ( thresh_dfct - Dfct[ term ] , EpsDfct ) )
        if( ! queue[ term ] ) {
         Index nxtnode = queue[ root ];
         if( LEZ( Pi[ nxtnode ] - Pi[ term ] , EpsCst )
             && ( root != nxtnode ) ) {
          queue[ root ] = term;
          queue[ term ] = nxtnode;
          }
         else {
          queue[ prvnde ] = term;
          queue[ term ] = root;
          prvnde = term;
          }
         }

       /* If root has a large enough surplus, keep it in */
       /*---- the queue and return for another iteration */

       if( GTZ( thresh_dfct - Dfct[ root ] , EpsDfct ) ) {
        prvnde = root;
        root = queue[ root ];
        cond200 = true;
        break;
        }
       else
        salto = 3000;
       }
      } while( salto != 3000 );
     } while( cond500 );
    }

   if( cond200 )
    continue;

   /* End of augmenting cycle */
   /* Check for termination of scaling phase. If scaling phase is not */
   /* finished, advance the queue  and return to take another node. */

   Index nxtnode = queue[ root ];
   if( root != nxtnode ) {
    queue[ root ] = 0;
    queue[ prvnde ] = root = nxtnode;
    cond200 = true;
    }
   } while( cond200 );

  /*---------- End of subproblem (scaling phase) ----------*/
  /*-------------------- Reduce epsilon -------------------*/

  eps = eps / factor;
  if( eps < 1 )
   eps = 1;

  if( eps == 1 )
   thresh_dfct = 0;
  else
   thresh_dfct = thresh_dfct / FNumber( factor );

  /*----- If another auction scaling phase remains, reset -----*/
  /*----- the flows & the push lists else reset arc flows -----*/
  /*--------- to satisfy cs and compute reduced costs ---------*/

  if( crash ) {
   #if( DYNMC_MCF_RIV > 1 )
    cFRow tCap = Cap;
   #endif
   FRow tU = U;
   FRow tX = X;
   cCRow tRC = RC;
   cIndex_Set tEndn = Endn;
   cIndex_Set tStartn = Startn;
   for( Index i = m ; i-- ; ) {
    tRC++; tU++; tX++;
    tStartn++; tEndn++;
    #if( DYNMC_MCF_RIV > 1 )
     if( ! *(++tCap) )
      continue;
    #endif

    Index end = *tEndn;
    Index start = *tStartn;
    CNumber dRC = Pi[ start ] - Pi[ end ] - (*tRC);
    if( GTZ( dRC - eps , EpsCst ) ) {
     FNumber resid = *tU;
     if( GTZ( resid , EpsFlw ) ) {
      Dfct[ start ] += resid;
      Dfct[ end ] -= resid;
      *tX += resid;
      *tU = 0;
      }
     }
    else
     if( LTZ( dRC + eps , EpsCst ) ) {
      FNumber flow = *tX;
      if( GTZ( flow , EpsFlw ) ) {
       Dfct[ start ] -= flow;
       Dfct[ end ] += flow;
       *tU += flow;
       *tX = 0;
       }
      }
    }

   cond100 = true;
   }
  } while( cond100 );

 crash = true;

 FRow tX = X + m;
 FRow tU = U + m;
 CRow tRC = RC + m;
 Index_Set tEn = Endn + m;
 Index_Set tSn = Startn + m;
 for( ; tX > X ; tX-- , tU-- , tRC-- , tEn-- , tSn-- ) {
  #if( ( ! SAME_GRPH_RIV ) || DYNMC_MCF_RIV )
   if( *tRC == Inf<CNumber>() )
    continue;
  #endif

  cIndex en = *tEn;
  cIndex sn = *tSn;
  CNumber ttRC = ( *tRC += Pi[ en ] - Pi[ sn ] );
  if( LTZ( ttRC , EpsCst ) ) {
   FNumber ttU = *tU;
   if( GTZ( ttU , EpsFlw ) ) {
    Dfct[ sn ] += ttU;
    Dfct[ en ] -= ttU;
    *tX += ttU;
    *tU = 0;
    }
   }
  else
   if( GTZ( ttRC , EpsCst ) ) {
    FNumber ttX = *tX;
    if( GTZ( ttX , EpsFlw ) ) {
     Dfct[ sn ] -= ttX;
     Dfct[ en ] += ttX;
     *tU += ttX;
     *tX = 0;
     }
    }
  }

 }  // end( auction )

#endif // ( AUCTION )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::MemAlloc( void )
{
 #if( SAME_GRPH_RIV )
  if( ! Startn )  // allocating data structures for graph topology- - - - - -
 #endif
  {
   Startn = new Index[ mmax ];
   Endn   = new Index[ mmax ];

   #if( SAME_GRPH_RIV )
    *Startn = 0;  // 0 is not a feasible value for a Startn[ 1 ]; this tells
                  // that Startn and Endn have not been initialized yet
   #endif
   Startn--;
   Endn--;
   }

 #if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
  if( ! FIn )     // allocating data structures for FS / BS - - - - - - - - -
 #endif
  {
   FIn   = new Index[ nmax ];
   NxtIn = new Index[ mmax ];

   FOu   = new Index[ nmax ];
   NxtOu = new Index[ mmax ];

   #if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
    *FIn = Inf<Index>();  // Inf<Index>() is not a feasible value for FIn[ 1 ]; this tells
                   // that FS and BS have not been initialized yet
   #endif
   FIn--;
   NxtIn--;
   FOu--;
   NxtOu--;
   }

 if( maxnmax < nmax ) {  // allocating node-wise temporaries- - - - - - - - -
  maxnmax = nmax;

  if( mark ) {
   #if( AUCTION )
    delete[] ++SB_arc;
    delete[] ++FpushB;
    delete[] ++extend_arc;
    delete[] ++SB_level;
   #endif

   #if( P_ALLOC )
    delete[] ++Pi;
   #endif

   delete[] ++DDNeg;
   delete[] ++DDPos;
   delete[] ++queue;
   delete[] ++scan;
   delete[] ++mark;
   delete[] label;
   delete[] ++Prdcsr;
   }

  Prdcsr = new SIndex[ nmax ];
  label  = new Index[ nmax ];
  mark   = new bool[ nmax ];
  scan   = new bool[ nmax ];
  queue  = new Index[ nmax ];
  DDPos  = new FNumber[ nmax ];
  DDNeg  = new FNumber[ nmax ];

  Prdcsr--;
  mark--;
  scan--;
  queue--;
  DDPos--;
  DDNeg--;

  #if( P_ALLOC )
   Pi = new CNumber[ nmax ];
   Pi--;
  #endif

  #if( AUCTION )
   SB_level   = new CNumber[ nmax ];
   extend_arc = new SIndex[ nmax ];

   FpushF = label;

   FpushB = new Index[ nmax ];
   SB_arc = new SIndex[ nmax ];

   SB_level--;
   extend_arc--;
   FpushF--;
   FpushB--;
   SB_arc--;
  #endif

  }  // end if( allocating node-wise temporaries )

 if( maxmmax < mmax ) {  // allocating arc-wise temporaries - - - - - - - - -
  maxmmax = mmax;

  if( save ) {
   #if( AUCTION )
    delete[] ++NxtpushB;
    delete[] ++NxtpushF;
   #endif

   delete[] save;
   }

  save = new Index[ mmax ];

  #if( AUCTION )
   NxtpushF = new Index[ mmax ];
   NxtpushB = new Index[ mmax ];

   NxtpushF--;
   NxtpushB--;
  #endif

  }  // end if( allocating arc-wise temporaries )

 // allocating flows, reduced costs, potentials etc - - - - - - - - - - - - -

 X    = new FNumber[ mmax ];
 U    = new FNumber[ mmax ];
 Cap  = new FNumber[ mmax ];
 RC   = new CNumber[ mmax ];
 C    = new CNumber[ mmax ];
 Dfct = new FNumber[ nmax ];
 B    = new FNumber[ nmax ];

 U--;
 Cap--;
 RC--;
 C--;
 X--;
 Dfct--;
 B--;

 // allocating "restricted" (to balanecd arcs) BS and FS- - - - - - - - - - -

 tfstin = new Index[ nmax ];
 tnxtin = new Index[ mmax ];

 tfstou = new Index[ nmax ];
 tnxtou = new Index[ mmax ];

 tfstin--;
 tnxtin--;
 tfstou--;
 tnxtou--;

 }  // end( MemAlloc )

/*--------------------------------------------------------------------------*/

inline void RelaxIV::MemDeAlloc( void )
{
 // deallocating "local" vectors- - - - - - - - - - - - - - - - - - - - - - -

 if( mmax && nmax ) {
  delete[] ++tnxtou;
  delete[] ++tfstou;

  delete[] ++tnxtin;
  delete[] ++tfstin;

  delete[] ++B;
  delete[] ++Dfct;
  delete[] ++C;
  delete[] ++RC;
  delete[] ++Cap;
  delete[] ++U;
  delete[] ++X;
  }

 #if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
  if( ! InstCntr )  // deallocating FS and BS data structures - - - - - - - -
 #else
  if( mmax && nmax )
 #endif
  {
   delete[] ++NxtOu;
   delete[] ++FOu;

   delete[] ++NxtIn;
   delete[] ++FIn;
   #if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
    FIn = NULL;
   #endif
   }

 #if( SAME_GRPH_RIV )
  if( ! InstCntr )  // deallocating graph topology description- - - - - - - -
 #else
  if( mmax && nmax )
 #endif
  {
   delete[] ++Endn;
   delete[] ++Startn;

   #if( SAME_GRPH_RIV )
    Startn = NULL;
   #endif
   }

 }  // end( MemDeAlloc )

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC MEMBERS -------------------------------*/
/*--------------------------------------------------------------------------*/

MCFClass::Index RelaxIV::InstCntr = 0;
MCFClass::Index RelaxIV::maxnmax = 0;
MCFClass::Index RelaxIV::maxmmax = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#if( SAME_GRPH_RIV )
 MCFClass::Index_Set RelaxIV::Startn = NULL;
 MCFClass::Index_Set RelaxIV::Endn = NULL;
#endif

#if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
 MCFClass::Index_Set RelaxIV::NxtOu = NULL;
 MCFClass::Index_Set RelaxIV::NxtIn = NULL;
 MCFClass::Index_Set RelaxIV::FOu = NULL;
 MCFClass::Index_Set RelaxIV::FIn = NULL;
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#if( P_ALLOC )
 RelaxIV *RelaxIV::PiOwnr = NULL;
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

MCFClass::CRow       RelaxIV::Pi = NULL;

RelaxIV::Bool_Vec    RelaxIV::mark = NULL;
MCFClass::Index_Set  RelaxIV::save = NULL;
MCFClass::Index_Set  RelaxIV::label = NULL;
RelaxIV::SIndex_Set  RelaxIV::Prdcsr = NULL;

RelaxIV::Bool_Vec    RelaxIV::scan = NULL;

MCFClass::Index_Set  RelaxIV::queue = NULL;
MCFClass::Index      RelaxIV::lastq = 0;
MCFClass::Index      RelaxIV::prvnde = 0;

MCFClass::FRow       RelaxIV::DDPos = NULL;
MCFClass::FRow       RelaxIV::DDNeg = NULL;
#if( AUCTION )
 MCFClass::CRow      RelaxIV::SB_level = NULL;
 RelaxIV::SIndex_Set RelaxIV::extend_arc = NULL;
 RelaxIV::SIndex_Set RelaxIV::SB_arc = NULL;
 MCFClass::Index_Set RelaxIV::FpushF = NULL;
 MCFClass::Index_Set RelaxIV::NxtpushF = NULL;
 MCFClass::Index_Set RelaxIV::FpushB = NULL;
 MCFClass::Index_Set RelaxIV::NxtpushB = NULL;
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- End File RelaxIV.C -----------------------------*/
/*--------------------------------------------------------------------------*/
