/*--------------------------------------------------------------------------*/
/*---------------------------- File MCFSimplex.C ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Linear and Quadratic Min Cost Flow problems solver based on the      --*/
/*-- (primal and dual) simplex algorithm. Conforms to the standard MCF    --*/
/*-- interface defined in MCFClass.h.                                     --*/
/*--                                                                      --*/
/*--                            VERSION 1.00                              --*/
/*--                           29 - 08 - 2011                             --*/
/*--                                                                      --*/
/*--                           Implementation:                            --*/
/*--                                                                      --*/
/*--                         Alessandro Bertolini                         --*/
/*--                          Antonio Frangioni                           --*/
/*--                                                                      --*/
/*--                       Operations Research Group                      --*/
/*--                      Dipartimento di Informatica                     --*/
/*--                         Universita' di Pisa                          --*/
/*--                                                                      --*/
/*-- Copyright (C) 2008 - 2011 by Alessandro Bertolini, Antonio Frangioni --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFSimplex.h"
#include <algorithm>
#include <iostream>

#include <cstdlib>
#include <ctime>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace MCFClass_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LIMITATEPRECISION 1

/* If LIMITATEPRECISION is 1, in the quadratic case the Primal Simplex accepts
   entering arc in base only if the decrease of the o.f. value is bigger than 
   a fixed thresold (EpsOpt * oldFOValue / n). Otherwise, any strict decrease
   in the o.f. value is accepted. */
   
#define UNIPI_PRIMAL_INITIAL_SHOW 0

/* If UNIPI_PRIMAL_INITIAL_SHOW = 1, Primal Simplex shows the initial condition
   (arcs and nodes) of the network. */

#define UNIPI_PRIMAL_ITER_SHOW 0

/* If UNIPI_PRIMAL_FINAL_SHOW = x with x > 0, Primal Simplex shows the condition
   (arcs and nodes) of the network every x iterations. */

#define UNIPI_PRIMAL_FINAL_SHOW 0

/* If UNIPI_PRIMAL_FINAL_SHOW = 1, Primal Simplex shows the final condition
   (arcs and nodes) of the network. */

#define UNIPI_DUAL_INITIAL_SHOW 0

/* If UNIPI_DUAL_INITIAL_SHOW = 1, Dual Simplex shows the initial condition
   (arcs and nodes) of the network. */

#define UNIPI_DUAL_ITER_SHOW 0

/* If UNIPI_DUAL_FINAL_SHOW = x with x > 0, Dual Simplex shows the condition
   (arcs and nodes) of the network every x iterations. */

#define UNIPI_DUAL_FINAL_SHOW 0

/* If UNIPI_DUAL_FINAL_SHOW = 1, Dual Simplex shows the final condition
   (arcs and nodes) of the network. */

#define UNIPI_VIS_DUMMY_ARCS 1

/* If UNIPI_VIS_DUMMY_ARCS = 1, Primal Simplex or Dual Simplex shows the
   conditions of the dummy arcs. */

#define UNIPI_VIS_ARC_UPPER 1
#define UNIPI_VIS_ARC_COST 1
#define UNIPI_VIS_ARC_Q_COST 1
#define UNIPI_VIS_ARC_REDUCT_COST 1
#define UNIPI_VIS_ARC_STATE 1
#define UNIPI_VIS_NODE_BASIC_ARC 1

/* When Primal Simplex or Dual Simplex shows the conditions of the network, 
   for every arcs the algorithm shows the flow, for every nodes it shows
   balance and potential. These 6 flags decide if the algorithm shows a
   particular value of the arcs/nodes;  for example if
   UNIPI_VIS_ARC_UPPER == 1, the algorithm shows the upper bounds of all
   arcs. */

#define FOSHOW 0

/* If FOSHOW is 1, the algorithm shows the f.o. value every x iterations 
   (x = UNIPI_PRIMAL_ITER_SHOW or x = UNIPI_DUAL_ITER_SHOW). */

#define OPTQUADRATIC 0

/* If OPTQUADRATIC is 1 the Primal Simplex, in the quadratic case, tries to
   optimize the update of the potential.
   Unfortunately this doesn't work well: for this reason it is set to 0. */

/*--------------------------------------------------------------------------*/
/*--------------------------- FUNCTIONS ------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline T ABS( const T x )
{
 return( x >= T( 0 ) ? x : - x );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline void Swap( T &v1 , T &v2 )
{
 T temp = v1;
 v1 = v2;
 v2 = temp;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- CONSTANTS ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( QUADRATICCOST == 0 )
 static const char DELETED  = -2;  // ident for deleted arcs
 static const char CLOSED   = -1;  // ident for closed arcs
 static const char BASIC    =  0;  // ident for basis arcs
 static const char AT_LOWER =  1;  // ident for arcs in L
 static const char AT_UPPER =  2;  // ident for arcs in U
#endif

/* These macros will be used by method MemAllocCandidateList() to set the
   values of numCandidateList and hotListSize. There are different macros,
   according to:

   - the used Simplex
   - the size of the network 
   - (obviously) the different variables

   This set of values tries to improve the performance of the two algorithms
   according to diversified situations. */

static const int PRIMAL_LOW_NUM_CANDIDATE_LIST =  30;
static const int PRIMAL_MEDIUM_NUM_CANDIDATE_LIST =  50;
static const int PRIMAL_HIGH_NUM_CANDIDATE_LIST =  200;
static const int PRIMAL_LOW_HOT_LIST_SIZE =  5;
static const int PRIMAL_MEDIUM_HOT_LIST_SIZE =  10;
static const int PRIMAL_HIGH_HOT_LIST_SIZE =  20;
static const int DUAL_LOW_NUM_CANDIDATE_LIST =  6;
static const int DUAL_HIGH_NUM_CANDIDATE_LIST =  10;
static const int DUAL_LOW_HOT_LIST_SIZE =  1;
static const int DUAL_HIGH_HOT_LIST_SIZE =  2;

/*--------------------------------------------------------------------------*/
/*--------------------------- COSTRUCTOR -----------------------------------*/
/*--------------------------------------------------------------------------*/

MCFSimplex::MCFSimplex( cIndex nmx , cIndex mmx )
            :
            MCFClass( nmx , mmx )
{
 #if( QUADRATICCOST )
  if( numeric_limits<FNumber>::is_integer )
   throw( MCFException( "FNumber must be float if QUADRATICCOST == 1" ) );

  if( numeric_limits<CNumber>::is_integer )
   throw( MCFException( "CNumber must be float if QUADRATICCOST == 1" ) );
 #endif

 newSession = true;
 if( nmax && mmax )
  MemAlloc();
 else
  nmax = mmax = 0;

 #if( QUADRATICCOST )
  recomputeFOLimits = 100;
  // recomputeFOLimits represents the limit of the iteration in which 
  // quadratic Primal Simplex computes "manually" the f.o. value
  EpsOpt = 1e-13;
  // EpsOpt is the fixed precision of the quadratic Primal Simplex
 #endif

 pricingRule = kCandidateListPivot;
 forcedNumCandidateList = 0;
 forcedHotListSize = 0;
 usePrimalSimplex = true;
 nodesP = NULL;
 nodesD = NULL;
 arcsP = NULL;
 arcsD = NULL;
 candP = NULL;
 candD = NULL;

 if( numeric_limits<CNumber>::is_integer )
  MAX_ART_COST = CNumber( 1e7 );
 else
  MAX_ART_COST = CNumber( 1e10 );

 }  // end( MCFSimplex )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void MCFSimplex::LoadNet( cIndex nmx , cIndex mmx , cIndex pn , cIndex pm ,
                          cFRow pU , cCRow pC , cFRow pDfct ,
                          cIndex_Set pSn , cIndex_Set pEn )
{
 MemDeAllocCandidateList();

 if( ( nmx != nmax ) || ( mmx != mmax ) ) {
  // if the size of the allocated memory changes

  if( nmax && mmax )  {  // if the memory was already allocated
   MemDeAlloc(true);         // deallocate the Primal 
   MemDeAlloc(false);        // and the Dual data structures
   nmax = mmax = 0;
   }

  if( nmx && mmx ) {   // if the new dimensions of the memory are correct
   nmax = nmx;
   mmax = mmx;
   MemAlloc();
   }
  }

 if( ( ! nmax ) || ( ! mmax ) )
  // if one of the two dimension of the memory isn't correct
  nmax = mmax = 0;

 if( nmax ) {  // if the new dimensions of the memory are correct
  n = pn;
  m = pm;

  if( usePrimalSimplex ) {
   stopNodesP = nodesP + n;
   dummyRootP = nodesP + nmax;
   for( nodePType *node = nodesP ; node != stopNodesP ; node++ )
    node->balance = pDfct[ node - nodesP ];  // initialize nodes

   stopArcsP = arcsP + m;
   dummyArcsP = arcsP + mmax;
   stopDummyP = dummyArcsP + n;
   for( arcPType *arc = arcsP ; arc != stopArcsP ; arc++ ) {
    // initialize real arcs
    arc->cost = pC[ arc - arcsP ];
    #if( QUADRATICCOST )
     arc->quadraticCost = 0; 
    #endif
    arc->upper = pU[ arc - arcsP ];
    arc->tail = nodesP + pSn[ arc - arcsP ] - 1;
    arc->head = nodesP + pEn[ arc - arcsP ] - 1;
    }
   }
  else {
   stopNodesD = nodesD + n;
   dummyRootD = nodesD + nmax;
   for( nodeDType *node = nodesD ; node != stopNodesD ; node++ )
    node->balance = pDfct[ node - nodesD ];  // initialize nodes

   stopArcsD = arcsD + m;
   dummyArcsD = arcsD + mmax;
   stopDummyD = dummyArcsD + n;
   for( arcDType *arc = arcsD ; arc != stopArcsD ; arc++ ) {
    // initialize real arcs
    arc->cost = pC[ arc - arcsD ];
    #if( QUADRATICCOST )
     arc->quadraticCost = 0; 
    #endif
    arc->upper = pU[ arc - arcsD ];
    arc->tail = nodesD + pSn[ arc - arcsD ] - 1;
    arc->head = nodesD + pEn[ arc - arcsD ] - 1;
    }

   CreateAdditionalDualStructures();
   }

  if( pricingRule == kCandidateListPivot )
   MemAllocCandidateList();

  status = kUnSolved;
  }
 }  // end( MCFSimplex::LoadNet )
                
/*-------------------------------------------------------------------------*/

void MCFSimplex::SetAlg( bool UsPrml , char WhchPrc )
{
 bool newUsePrimalSimplex = UsPrml;
 bool oldUsePrimalSimplex = usePrimalSimplex;
 char newPricingRule = WhchPrc;
 char oldPricingRule = pricingRule;

 usePrimalSimplex = newUsePrimalSimplex;
 pricingRule = newPricingRule;

 if( ( ! usePrimalSimplex ) && ( pricingRule == kDantzig) ) {
  pricingRule = kFirstEligibleArc;
  newPricingRule = pricingRule;
  }

 if( ( newUsePrimalSimplex != oldUsePrimalSimplex ) ||
     ( newPricingRule != oldPricingRule ) ) {
  MemDeAllocCandidateList();

  if( newUsePrimalSimplex != oldUsePrimalSimplex ) {
   #if( QUADRATICCOST )
     throw(
      MCFException( "Primal Simplex is the only option if QUADRATICCOST == 1"
                    ) );
   }
   #else
    MemAlloc();
    nodePType *nP = nodesP;
    nodeDType *nD = nodesD;
    arcPType *aP = arcsP;
    arcDType *aD = arcsD;

    if( newUsePrimalSimplex ) { // from Dual to Primal
     if( nodesD == NULL )
      return;

     if( newSession )
      CreateInitialDualBase();

     dummyRootP = nodesP + nmax;
     stopNodesP = nodesP + n;
     dummyArcsP = arcsP + mmax;
     stopArcsP = arcsP + m;
     stopDummyP = dummyArcsP + n;
     // Copy the old Dual data structure in a new Primal data structure
     while( nD != stopNodesD ) {
      nP->prevInT = nodesP + ( nD->prevInT - nodesD );
      nP->nextInT = nodesP + ( nD->nextInT - nodesD );
      nP->enteringTArc = arcsP + ( nD->enteringTArc - arcsD );
      nP->balance = nD->balance;
      nP->potential = nD->potential;
      nP->subTreeLevel = nD->subTreeLevel;
      nP++;
      nD++;
      }

     dummyRootP->prevInT = NULL;
     dummyRootP->nextInT = nodesP + ( dummyRootD->nextInT - nodesD );
     dummyRootP->enteringTArc = arcsP + ( dummyRootD->enteringTArc - arcsD );
     dummyRootP->balance = dummyRootD->balance;
     dummyRootP->potential = dummyRootD->potential;
     dummyRootP->subTreeLevel = dummyRootD->subTreeLevel;
     while( aD != stopArcsD ) {
      aP->tail = nodesP + ( aD->tail - nodesD );
      aP->head = nodesP + ( aD->head - nodesD );
      aP->flow = aD->flow;
      aP->cost = aD->cost;
      aP->ident = aD->ident;
      aP->upper = aD->upper;
      aP++;
      aD++;
      }

     aP = dummyArcsP;
     aD = dummyArcsD;
     while( aD != stopDummyD ) {
      aP->tail = nodesP + ( aD->tail - nodesD );
      aP->head = nodesP + ( aD->head - nodesD );
      aP->flow = aD->flow;
      aP->cost = aD->cost;
      aP->ident = aD->ident;
      if( ( ETZ(aP->flow, EpsFlw) ) && (aP->ident == AT_UPPER) )
       aP->ident = AT_LOWER;

      aP->upper = Inf<FNumber>();
      aP++;
      aD++;
      }

     MemDeAlloc(false);
     if( Senstv && ( status != kUnSolved ) ) {
      nodePType *node = dummyRootP;
      for( int i = 0 ; i < n ; i++ )
       node = node->nextInT;

      node->nextInT = NULL;
      dummyRootP->prevInT = NULL;
      dummyRootP->enteringTArc = NULL;
      // balance the flow
      CreateInitialPModifiedBalanceVector();
      PostPVisit( dummyRootP );
      // restore the primal admissibility
      BalanceFlow( dummyRootP );
      ComputePotential( dummyRootP );
      }
     else
      status = kUnSolved;
     }
   else {  // from Primal to Dual
    if( nodesP == NULL )
     return;

    if( newSession )
     CreateInitialPrimalBase();

    dummyRootD = nodesD + nmax;
    stopNodesD = nodesD + n;
    dummyArcsD = arcsD + mmax;
    stopArcsD = arcsD + m;
    stopDummyD = dummyArcsD + n;
    // Copy the old Primal data structure in a new Dual data structure
    while( nP != stopNodesP ) {
     nD->prevInT = nodesD + ( nP->prevInT - nodesP );
     nD->nextInT = nodesD + ( nP->nextInT - nodesP );
     nD->enteringTArc = arcsD + ( nP->enteringTArc - arcsP );
     nD->balance = nP->balance;
     nD->potential = nP->potential;
     nD->subTreeLevel = nP->subTreeLevel;
     nP++;
     nD++;
     }

    dummyRootD->prevInT = NULL;
    dummyRootD->nextInT = nodesD + ( dummyRootP->nextInT - nodesP );
    dummyRootD->enteringTArc = NULL;
    dummyRootD->balance = dummyRootP->balance;
    dummyRootD->potential = dummyRootP->potential;
    dummyRootD->subTreeLevel = dummyRootP->subTreeLevel;
    while( aP != stopArcsP ) {
     aD->tail = nodesD + ( aP->tail - nodesP );
     aD->head = nodesD + ( aP->head - nodesP );
     aD->flow = aP->flow;
     aD->cost = aP->cost;
     aD->ident = aP->ident;
     aD->upper = aP->upper;
     aP++;
     aD++;
     }

    aP = dummyArcsP;
    aD = dummyArcsD;
    while( aP != stopDummyP ) {
     aD->tail = nodesD + ( aP->tail - nodesP );
     aD->head = nodesD + ( aP->head - nodesP );
     aD->flow = aP->flow;
     aD->cost = aP->cost;
     aD->ident = aP->ident;
     aD->upper = 0;
     aP++;
     aD++;
     }

    CreateAdditionalDualStructures();
    MemDeAlloc(true);
    nodeDType *node = dummyRootD;
    for( int i = 0 ; i < n ; i++ )
     node = node->nextInT;

    node->nextInT = NULL;
    dummyRootD->enteringTArc = NULL;
    dummyRootD->prevInT = NULL;
    if( Senstv && ( status != kUnSolved ) ) {
     // fix every flow arc according to its reduct cost
     for( arcDType *arc = arcsD ; arc != stopArcsD ; arc++ ) {
      if( (arc->ident == AT_LOWER) && LTZ( ReductCost( arc ) , EpsCst ) ) {
       arc->flow = arc->upper;
       arc->ident = AT_UPPER;
       }

      if( ( arc->ident == AT_UPPER ) && GTZ( ReductCost( arc ) , EpsCst ) ) {
       arc->flow = 0;
       arc->ident = AT_LOWER;
       }
      }

     // balance the flow
     CreateInitialDModifiedBalanceVector();
     PostDVisit( dummyRootD );
     }
    else
     status = kUnSolved;
    }
   //#endif
   }
#endif
 if( newPricingRule == kFirstEligibleArc )
  if( newUsePrimalSimplex )
   arcToStartP = arcsP;
  else
   arcToStartD = arcsD;

 if( ( nmax && mmax ) && ( newPricingRule == kCandidateListPivot ) )
  MemAllocCandidateList();
  }
 }  // end( SetAlg )

/*-------------------------------------------------------------------------*/

void MCFSimplex::SetPar( int par, int val )
{
 switch( par ) {
 case kAlgPrimal:
  if( val == kYes )
   SetAlg( true , pricingRule);

  if( val == kNo )
   SetAlg( false, pricingRule );

  break;

 case kAlgPricing:
 
  if( ( val == kDantzig ) || ( val == kFirstEligibleArc ) ||
      ( val == kCandidateListPivot ) )
   SetAlg( usePrimalSimplex , val );

  break;

 case kNumCandList:

  MemDeAllocCandidateList();
  forcedNumCandidateList = val;
  MemAllocCandidateList();
  forcedNumCandidateList = 0;
  forcedHotListSize = 0;
  break;

 case kHotListSize:

  MemDeAllocCandidateList();
  forcedHotListSize = val;
  MemAllocCandidateList();
  forcedNumCandidateList = 0;
  forcedHotListSize = 0;
  break;

 case kRecomputeFOLimits:

  recomputeFOLimits = val;
  break;

 default:

  MCFClass::SetPar(par, val);
 }
 }  // end( SetPar( int ) )

/*-------------------------------------------------------------------------*/

void MCFSimplex::SetPar( int par , double val )
{
 switch( par ) {
 case kEpsOpt:

  EpsOpt = val;
  break;

 default:
  MCFClass::SetPar( par , val );
 }
 }  // end( SetPar( double )

/*-------------------------------------------------------------------------*/
/*--------------- METHODS FOR SOLVING THE PROBLEM -------------------------*/
/*-------------------------------------------------------------------------*/

void MCFSimplex::SolveMCF( void )
{
 if( MCFt )
  MCFt->Start();

 if( status == kUnSolved )
  if(usePrimalSimplex )
   CreateInitialPrimalBase();
  else
   CreateInitialDualBase();

 newSession = false;
 if( usePrimalSimplex )
  PrimalSimplex();
 else
  DualSimplex();

 if( MCFt )
  MCFt->Stop();

 }  // end( MCFSimplex::SolveMCF )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void MCFSimplex::MCFGetX( FRow F , Index_Set nms , cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  if( usePrimalSimplex )
   for( Index i = strt ; i < stp ; i++ ) {
    FNumber tXi = ( arcsP + i )->flow;
    if( GTZ( tXi , EpsFlw ) ) {
     *(F++) = tXi;
     *(nms++) = i;
     }
    }
  else
   for( Index i = strt ; i < stp ; i++ ) {
    FNumber tXi = ( arcsD + i )->flow;
    if( GTZ( tXi , EpsFlw ) ) {
     *(F++) = tXi;
     *(nms++) = i;
     }
    }

  *nms = Inf<Index>();
  }        
 else
  if( usePrimalSimplex )
   for( Index i = strt; i < stp; i++ )
    *(F++) = ( arcsP + i )->flow;
  else
   for( Index i = strt; i < stp; i++ )
    *(F++) = ( arcsD + i )->flow;

 }  // end( MCFSimplex::MCFGetX( some ) )

/*--------------------------------------------------------------------------*/

void MCFSimplex::MCFGetRC( CRow CR , cIndex_Set nms , cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  if( usePrimalSimplex )
   for( Index h ; ( h = *(nms++) ) < stp ; )
   *(CR++) = CNumber( ReductCost( arcsP + h ) );
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(CR++) = ReductCost( arcsD + h );
  }
 else {
  if( stp > m )
   stp = m;

  if( usePrimalSimplex )
   for( arcPType* arc = arcsP + strt ; arc < arcsP + stp ; arc++ )
    *(CR++) = CNumber( ReductCost( arc ) );
  else
   for( arcDType* arc = arcsD + strt ; arc < arcsD + stp ; arc++ )
    *(CR++) = ReductCost( arc );
  }
 }  // end( MCFSimplex::MCFGetRC( some ) )

/*--------------------------------------------------------------------------*/

MCFSimplex::CNumber MCFSimplex::MCFGetRC( cIndex i )
{
 if( usePrimalSimplex )
  return CNumber( ReductCost( arcsP + i ) );
 else
  return( ReductCost( arcsD + i ) );

 }  // end( MCFSimplex::MCFGetRC( i ) )

/*--------------------------------------------------------------------------*/

void MCFSimplex::MCFGetPi( CRow P , cIndex_Set nms , cIndex strt , Index stp )
{
 if( stp > n )
  stp = n;

 if( nms ) {
  while( *nms < strt )
   nms++;

  if( usePrimalSimplex )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(P++) = CNumber( (nodesP + h)->potential );
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(P++) = (nodesD + h )->potential;
  }
 else
  if( usePrimalSimplex )
   for( nodePType *node = nodesP + strt ; node < ( nodesP + stp ) ; node++ )
    *(P++) = CNumber( node->potential );
  else
   for( nodeDType *node = nodesD + strt ; node++ < ( nodesD + stp ) ; node++ )
    *(P++) = node->potential;

 }  // end(  MCFSimplex::MCFGetPi( some ) )

/*--------------------------------------------------------------------------*/

MCFSimplex::FONumber MCFSimplex::MCFGetFO( void )
{
 if( status == kOK )
  return( (FONumber) GetFO() );
 else
  if( status == kUnbounded ) 
   return( - Inf<FONumber>() );
  else
   return( Inf<FONumber>() );

 }  // end( MCFSimplex::MCFGetFO )

/*-------------------------------------------------------------------------*/
/*----------METHODS FOR READING THE DATA OF THE PROBLEM--------------------*/
/*-------------------------------------------------------------------------*/

void MCFSimplex::MCFArcs( Index_Set Startv , Index_Set Endv ,
			  cIndex_Set nms , cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  while( *nms < strt )
   nms++;

  if( usePrimalSimplex )
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    if( Startv )
     *(Startv++) = Index( (arcsP + h)->tail - nodesP) - USENAME0;
    if( Endv )
     *(Endv++) = Index( (arcsP + h)->head - nodesP ) - USENAME0;
    }
  else
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    if( Startv )
     *(Startv++) = Index( (arcsD + h)->tail - nodesD) - USENAME0;
    if( Endv )
     *(Endv++) = Index( (arcsD + h)->head - nodesD ) - USENAME0;
    }
  }
 else
  if( usePrimalSimplex )
   for( arcPType* arc = arcsP + strt ; arc < (arcsP + stp) ; arc++ ) {
    if( Startv )
     *(Startv++) = Index( arc->tail - nodesP ) - USENAME0;
    if( Endv )
     *(Endv++) = Index( arc->head - nodesP ) - USENAME0;
    }
  else
   for( arcDType* arc = arcsD + strt ; arc < (arcsD + stp) ; arc++ ) {
    if( Startv )
     *(Startv++) = Index( arc->tail - nodesD ) - USENAME0;
    if( Endv )
     *(Endv++) = Index( arc->head - nodesD ) - USENAME0;
    }

 }  // end( MCFSimplex::MCFArcs )

/*-------------------------------------------------------------------------*/

void MCFSimplex::MCFCosts( CRow Costv , cIndex_Set nms ,
			   cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  while( *nms < strt )
   nms++;

  if( usePrimalSimplex )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(Costv++) = (arcsP + h)->cost;
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(Costv++) = (arcsD + h)->cost;
  }
 else
  if( usePrimalSimplex )
   for( arcPType* arc = arcsP + strt ; arc < (arcsP + stp) ; arc++ )
    *(Costv++) = arc->cost;           
  else
   for( arcDType* arc = arcsD + strt ; arc < (arcsD + stp) ; arc++ )
    *(Costv++) = arc->cost;           

 }  // end( MCFSimplex::MCFCosts ( some ) )

/*-------------------------------------------------------------------------*/

void MCFSimplex::MCFQCoef( CRow Qv , cIndex_Set nms  ,
			   cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  while( *nms < strt )
   nms++;

  #if( QUADRATICCOST )
   if( usePrimalSimplex )
    for( Index h ; ( h = *(nms++) ) < stp ; )
     *(Qv++) = (arcsP + h)->quadraticCost;
   else
    for( Index h ; ( h = *(nms++) ) < stp ; )
     *(Qv++) = (arcsD + h)->quadraticCost;
  #else
    for( Index h ; ( h = *(nms++) ) < stp ; )
     *(Qv++) = 0;
  #endif
  }
 else
  #if( QUADRATICCOST )
   if( usePrimalSimplex )
    for( arcPType* arc = arcsP + strt ; arc < ( arcsP + stp ) ; arc++ )
     *(Qv++) = arc->quadraticCost;
   else
    for( arcDType* arc = arcsD + strt ; arc < ( arcsD + stp ) ; arc++ )
     *(Qv++) = arc->quadraticCost;
  #else
   for( Index h = strt ; h++ < stp ; )
    *(Qv++) = 0;           
  #endif

 }  // end( MCFSimplex::MCFQCoef ( some ) )

/*-------------------------------------------------------------------------*/

void MCFSimplex::MCFUCaps( FRow UCapv , cIndex_Set nms ,
			   cIndex strt , Index stp ) 
{
 if( stp > m )
  stp = m;

 if( nms ) {
  while( *nms < strt )
   nms++;

  if( usePrimalSimplex )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(UCapv++) = (arcsP + h)->upper;
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(UCapv++) = (arcsD + h)->upper;
  }
 else
  if( usePrimalSimplex )
   for( arcPType* arc = arcsP + strt ; arc <  (arcsP + stp ) ; arc++ )
    *(UCapv++) = arc->upper;
  else
   for( arcDType* arc = arcsD + strt ; arc < ( arcsD + stp ) ; arc++ )
    *(UCapv++) = arc->upper;

 }  // end( MCFSimplex::MCFUCaps ( some ) )
 
/*-------------------------------------------------------------------------*/

void MCFSimplex::MCFDfcts( FRow Dfctv , cIndex_Set nms ,
			   cIndex strt , Index stp )
{
 if( stp > n )
  stp = n;

 if( nms ) {
  while( *nms < strt )
   nms++;

  if( usePrimalSimplex )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(Dfctv++) = ( nodesP + h )->balance;
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(Dfctv++) = (nodesD + h )->balance;
  }
 else
  if( usePrimalSimplex )
   for( nodePType* node = nodesP + strt ; node < ( nodesP + stp ) ; node++ )
    *(Dfctv++) = node->balance;
  else
   for( nodeDType* node = nodesD + strt ; node < ( nodesD + stp ) ; node++ )
    *(Dfctv++) = node->balance;

 }  // end( MCFSimplex::MCFDfcts )

/*-------------------------------------------------------------------------*/
/*--------- METHODS FOR ADDING / REMOVING / CHANGING DATA -----------------*/
/*-------------------------------------------------------------------------*/

void MCFSimplex::ChgCosts( cCRow NCost , cIndex_Set nms ,
			   cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NCost++;
   }

  cIndex_Set tnms = nms;  // nms may be needed below
  if( usePrimalSimplex )
   for( Index h ; ( h = *(tnms++) ) < stp ; )
    arcsP[ h ].cost = *(NCost++);
  else
   for( Index h ; ( h = *(tnms++) ) < stp ; )
    arcsD[ h ].cost = *(NCost++);
  }
 else
  if( usePrimalSimplex )
   for( arcPType *arc = arcsP + strt ; arc < (arcsP + stp) ; arc++ )
    arc->cost = *(NCost++); 
  else
   for( arcDType *arc = arcsD + strt ; arc < (arcsD + stp) ; arc++ )
    arc->cost = *(NCost++); 

 if( Senstv && ( status != kUnSolved ) )
  if( usePrimalSimplex )
   ComputePotential( dummyRootP );
  else {
   #if( QUADRATICCOST == 0 )
    for( arcDType *arc = arcsD ; arc != stopArcsD ; arc++ )
     if( arc->ident > BASIC )
      if( GTZ( ReductCost( arc ) , EpsCst ) ) {
       arc->flow = 0;
       arc->ident = AT_LOWER;
       }
      else {
       arc->flow = arc->upper; 
       arc->ident = AT_UPPER;
       }

    CreateInitialDModifiedBalanceVector();
    PostDVisit( dummyRootD );
   #endif
   }
 else
  status = kUnSolved;

 }  // end( MCFSimplex::ChgCosts )

/*-------------------------------------------------------------------------*/

void MCFSimplex::ChgCost( Index arc , cCNumber NCost )
{
 if( arc >= m )
  return;

 if( usePrimalSimplex )
  ( arcsP + arc )->cost = NCost;
 else
  ( arcsD + arc )->cost = NCost;

 if( Senstv && ( status != kUnSolved ) ) {
  if( usePrimalSimplex ) {
   nodePType *node = ( arcsP + arc )->tail;
   if( ( ( arcsP + arc )->head)->subTreeLevel < node->subTreeLevel )
    node = ( arcsP + arc )->head;

   ComputePotential( dummyRootP );
   }
  else {
   #if( QUADRATICCOST == 0 )
    nodeDType *node = ( arcsD + arc )->tail;
    if( ( ( arcsD + arc )->head )->subTreeLevel < node->subTreeLevel )
     node = ( arcsD + arc )->head;

    ComputePotential( dummyRootD );
    for( arcDType *a = arcsD ; a != stopArcsD ; a++)
     if( a->ident > BASIC )
      if( GTZ( ReductCost( a ) , EpsCst ) ) {
       a->flow = 0;
       a->ident = AT_LOWER;
       }
      else {
       a->flow = a->upper; 
       a->ident = AT_UPPER;
       }

    CreateInitialDModifiedBalanceVector();
    PostDVisit( dummyRootD );
   #endif
   }
  }
 else
  status = kUnSolved;

 }  // end( MCFSimplex::ChgCost )

/*-------------------------------------------------------------------------*/

void MCFSimplex::ChgQCoef( cCRow NQCoef , cIndex_Set nms ,
			   cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 #if( QUADRATICCOST )
  if( nms ) {
   while( *nms < strt ) {
    nms++;
    NQCoef++;
    }

   cIndex_Set tnms = nms;  // nms may be needed below
   if( usePrimalSimplex )
    for( Index h ; ( h = *(tnms++) ) < stp ; )
     arcsP[ h ].quadraticCost = *(NQCoef++);
   else
    for( Index h ; ( h = *(tnms++) ) < stp ; )
     arcsD[ h ].quadraticCost = *(NQCoef++);
   }
  else
   if( usePrimalSimplex )
    for( arcPType *arc = arcsP + strt ; arc < ( arcsP + stp ) ; arc++ )
     arc->quadraticCost = *(NQCoef++); 
   else
    for( arcDType *arc = arcsD + strt ; arc < ( arcsD + stp ) ; arc++ )
     arc->quadraticCost = *(NQCoef++); 

  if( Senstv && (status != kUnSolved ) )
   ComputePotential( dummyRootP );
  else
   status = kUnSolved;
 #else
  if( NQCoef )
   throw( MCFException( "ChgQCoef: nonzero coefficients not allowed" ) );
 #endif
}  // end( MCFSimplex::ChgQCoef )

/*-------------------------------------------------------------------------*/

void MCFSimplex::ChgQCoef( Index arc , cCNumber NQCoef )
{
 #if( QUADRATICCOST )
  if( arc >= m )
   return;

  if( usePrimalSimplex )
   ( arcsP + arc )->quadraticCost = NQCoef;
  else
   ( arcsD + arc )->quadraticCost = NQCoef;

  if( Senstv && ( status != kUnSolved ) ) {
   nodePType *node = ( arcsP + arc )->tail;
   if( ( ( arcsP + arc )->head )->subTreeLevel < node->subTreeLevel )
    node = ( arcsP + arc )->head;

   ComputePotential( node );
   }
  else
   status = kUnSolved;
 #else
  if( NQCoef )
   throw( MCFException( "ChgQCoef: nonzero coefficients not allowed" ) );
 #endif
 }  // end( MCFSimplex::ChgQCoef )

/*-------------------------------------------------------------------------*/
    
void MCFSimplex::ChgDfcts( cFRow NDfct , cIndex_Set nms ,
			   cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NDfct++;
   }

  cIndex_Set tnms = nms;  // nms may be needed below
  if( usePrimalSimplex )
   for( Index h ; ( h = *(tnms++) ) < stp ; )
    nodesP[ h ].balance = *(NDfct++);
  else
   for( Index h ; ( h = *(tnms++) ) < stp ; )
    nodesD[ h ].balance = *(NDfct++);
  }
 else
  if( usePrimalSimplex )
   for( nodePType *node = nodesP + strt ; node < ( nodesP + stp ) ; node++ )
    node->balance = *(NDfct++);
  else
   for( nodeDType *node = nodesD + strt ; node < ( nodesD + stp ) ; node++ )
    node->balance = *(NDfct++);

 if( Senstv && (status != kUnSolved ) )
  if( usePrimalSimplex ) {
   CreateInitialPModifiedBalanceVector();
   PostPVisit( dummyRootP );
   BalanceFlow( dummyRootP );
   ComputePotential( dummyRootP );
   }
  else {
   CreateInitialDModifiedBalanceVector();
   PostDVisit( dummyRootD );
   }
 else
  status = kUnSolved;

 }  // end( MCFSimplex::ChgDfcts )

/*-------------------------------------------------------------------------*/

void MCFSimplex::ChgDfct( Index nod , cFNumber NDfct )
{ 
 if( nod > n )
  return;

 if( usePrimalSimplex )
  ( nodesP + nod - 1 )->balance = NDfct;
 else
  ( nodesD + nod - 1 )->balance = NDfct;

 if( Senstv && (status != kUnSolved ) )
  if( usePrimalSimplex ) {
   CreateInitialPModifiedBalanceVector();
   PostPVisit( dummyRootP );
   BalanceFlow( dummyRootP );
   ComputePotential( dummyRootP );
   }
  else {
   CreateInitialDModifiedBalanceVector();
   PostDVisit( dummyRootD );
   }
 else
  status = kUnSolved;

 }  // end( MCFSimplex::ChgDfct )

/*-------------------------------------------------------------------------*/

void MCFSimplex::ChgUCaps( cFRow NCap , cIndex_Set nms ,
			   cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NCap++;
   }

  cIndex_Set tnms = nms;  // nms may be needed below
  if( usePrimalSimplex )
   for( Index h ; ( h = *(tnms++) ) < stp ; )
    arcsP[ h ].upper = *(NCap++);
  else
   for( Index h ; ( h = *(tnms++) ) < stp ; )
    arcsD[ h ].upper = *(NCap++);
  }
 else
  if( usePrimalSimplex )
   for( arcPType *arc = arcsP + strt ; arc < ( arcsP + stp ) ; arc++ )
    arc->upper = *(NCap++);
  else
   for( arcDType *arc = arcsD + strt ; arc < ( arcsD + stp ) ; arc++ )
    arc->upper = *(NCap++); 

 if( Senstv && (status != kUnSolved ) ) {
  if( usePrimalSimplex ) {
   for( arcPType *arc = arcsP ; arc != stopArcsP ; arc++)
    #if( QUADRATICCOST )
     if( GT( arc->flow , arc->upper , EpsFlw ) ) 
      arc->flow = arc->upper;
    #else
     if( GT(arc->flow , arc->upper , EpsFlw ) || 
	 ( ( arc->ident == AT_UPPER ) &&
	   ( ! ETZ( arc->flow - arc->upper , EpsFlw ) ) ) )
      arc->flow = arc->upper;
    #endif

   CreateInitialPModifiedBalanceVector();
   PostPVisit( dummyRootP );
   BalanceFlow( dummyRootP );
   ComputePotential( dummyRootP );
   }
  else {
   #if( QUADRATICCOST == 0 )
    for( arcDType *arc = arcsD ; arc != stopArcsD ; arc++ )
     if( ( GT( arc->flow , arc->upper , EpsFlw ) && ( arc->ident != BASIC ) ) ||
	 ( ( arc->ident == AT_UPPER ) &&
	   ( ! ETZ( arc->flow - arc->upper , EpsFlw ) ) ) )
      arc->flow = arc->upper;

    CreateInitialDModifiedBalanceVector();
    PostDVisit( dummyRootD );
    ComputePotential( dummyRootD );
   #endif
   }
  }
 else
  status = kUnSolved;

 }  // end( MCFSimplex::ChgUCaps )

/*-------------------------------------------------------------------------*/

void MCFSimplex::ChgUCap( Index arc , cFNumber NCap )
{
 if( arc >= m )
  return;

 if( usePrimalSimplex )
  ( arcsP + arc )->upper = NCap;
 else
  ( arcsD + arc )->upper = NCap;

 if( Senstv && (status != kUnSolved ) ) {
  if( usePrimalSimplex ) {
   #if( QUADRATICCOST )
    if( GT( ( arcsP + arc )->flow , ( arcsP + arc )->upper , EpsFlw ) ) 
     ( arcsP + arc )->flow = ( arcsP + arc )->upper;
   #else
    if( GT( ( arcsP + arc )->flow , ( arcsP + arc )->upper , EpsFlw ) ||
	( ( ( arcsP + arc )->ident == AT_UPPER ) &&
	  ( ! ETZ( ( arcsP + arc )->flow - ( arcsP + arc )->upper , EpsFlw ) ) ) )
     ( arcsP + arc )->flow = ( arcsP + arc )->upper;
   #endif

   CreateInitialPModifiedBalanceVector();
   PostPVisit( dummyRootP );
   BalanceFlow( dummyRootP );
   ComputePotential( dummyRootP );
   }
  else {
   #if( QUADRATICCOST == 0 )
    if( ( GT( ( arcsD + arc )->flow , ( arcsD + arc )->upper , EpsFlw ) &&
	  ( ( ( arcsD + arc )->ident != BASIC ) ) ) || 
	( ( ( arcsD + arc )->ident == AT_UPPER ) &&
	  ( ! ETZ( ( arcsD + arc )->flow - ( arcsD + arc )->upper , EpsFlw ) ) ) ) {
     ( arcsD + arc )->flow = ( arcsD + arc )->upper;
     ( arcsD + arc )->ident = AT_UPPER;
     }

    CreateInitialDModifiedBalanceVector();
    PostDVisit( dummyRootD );
    ComputePotential( dummyRootD );
   #endif
   }
  }
 else
  status = kUnSolved;

 }  // end( MCFSimplex::ChgUCap )

/*-------------------------------------------------------------------------*/

bool MCFSimplex::IsClosedArc( cIndex name )
{
 if( name >= m )
  return( false );

 #if( QUADRATICCOST )
  return( ( arcsP + name )->cost == Inf<CNumber>() );
 #else
  if( usePrimalSimplex )
   return( ( ( arcsP + name )->ident < BASIC) );
  else
   return( ( ( arcsD + name )->ident < BASIC) );
 #endif
 }

/*-------------------------------------------------------------------------*/

void MCFSimplex::CloseArc( cIndex name )
{
 if( name >= m )
  return;

 if( usePrimalSimplex ) {
  arcPType *arc = arcsP + name;
  #if( QUADRATICCOST )
   if( arc->cost == Inf<CNumber>() )
    return;

   arc->cost = Inf<CNumber>();
  #else
   if( arc->ident < BASIC )
    return;

   arc->ident = CLOSED;
  #endif

  arc->flow = 0;

  if( Senstv && ( status != kUnSolved ) ) {
   nodePType *node = NULL;
   if( (arc->tail)->enteringTArc == arc )
    node = arc->tail;

   if( (arc->head)->enteringTArc == arc )
    node = arc->head;

   if( node ) {
    node->enteringTArc = dummyArcsP + ( node - nodesP );
    nodePType *last = CutAndUpdateSubtree( node , -node->subTreeLevel + 1 );
    PasteSubtree( node , last , dummyRootP );
    node->enteringTArc = dummyArcsP + ( node - nodesP );
    }

   CreateInitialPModifiedBalanceVector();
   PostPVisit(dummyRootP);
   BalanceFlow(dummyRootP);                
   ComputePotential(dummyRootP);
   }
  else
   status = kUnSolved;
  }
 else {
  #if( QUADRATICCOST == 0 )
   arcDType *arc = arcsD + name;
   if( arc->ident < BASIC )
    return;

   arc->flow = 0;
   arc->ident = CLOSED;

   if( Senstv && ( status != kUnSolved ) ) {
    nodeDType *node = NULL;
    if( ( arc->tail )->enteringTArc == arc)
     node = arc->tail;

    if( ( arc->head )->enteringTArc == arc )
     node = arc->head;

    if( node ) {
     node->enteringTArc = dummyArcsD + ( node - nodesD );
     nodeDType *last = CutAndUpdateSubtree( node , -node->subTreeLevel + 1 );
     PasteSubtree( node , last , dummyRootD );
     node->enteringTArc = dummyArcsD + ( node - nodesD );
     ComputePotential( dummyRootD );

     for( arcDType *a = arcsD ; a != stopArcsD ; a++ )
      if( a->ident > BASIC )
       if( GTZ( ReductCost( a ) , EpsCst ) ) {
	a->flow = 0;
	a->ident = AT_LOWER;
        }
       else {
	a->flow = a->upper; 
	a->ident = AT_UPPER;
        }
     }

    CreateInitialDModifiedBalanceVector();
    PostDVisit(dummyRootD);
    ComputePotential(dummyRootD);
    }
   else
    status = kUnSolved;
  #endif
  }
 }  // end( MCFSimplex::CloseArc )

/*--------------------------------------------------------------------------*/

void MCFSimplex::DelNode( cIndex name )
{
 if( name >= n )
  return;

 if( usePrimalSimplex ) {
  nodePType *node = nodesP + name;
  nodePType *last = CutAndUpdateSubtree(node, -node->subTreeLevel);
  nodePType *n = node->nextInT;
  while( n ) {
   if( n->subTreeLevel == 1 ) 
    n->enteringTArc = dummyArcsP + ( n - nodesP );

   n = n->nextInT;
   }

  PasteSubtree( node , last , dummyRootP );
  n = node->nextInT;
  dummyRootP->nextInT = n;
  n->prevInT = dummyRootP;
  
  for( arcPType *arc = arcsP ; arc != stopArcsP ; arc++ ) {
   if( ( ( arc->tail ) == node) || ( ( arc->head ) == node ) ) {
    arc->flow = 0;
    #if( QUADRATICCOST )
     arc->cost = Inf<CNumber>();
    #else
     arc->ident = CLOSED;
    #endif
    }
   }

  for( arcPType *arc = dummyArcsP ; arc != stopDummyP ; arc++ ) {
   if( ( ( arc->tail ) == node ) || ( ( arc->head ) == node ) ) {
    arc->flow = 0;
    #if( QUADRATICCOST )
     arc->cost = Inf<CNumber>();
    #else
     arc->ident = CLOSED;
    #endif
    }
   }

  CreateInitialPModifiedBalanceVector();
  PostPVisit( dummyRootP );
  BalanceFlow( dummyRootP );
  ComputePotential( dummyRootP );
  }
 else {
  #if( QUADRATICCOST == 0 )
   nodeDType *node = nodesD + name;
   nodeDType *last = CutAndUpdateSubtree( node , -node->subTreeLevel );
   nodeDType *n = node->nextInT;
   while( n ) {
    if( n->subTreeLevel == 1 )
     n->enteringTArc = dummyArcsD + ( n - nodesD );

    n = n->nextInT;
    }

   PasteSubtree( node , last , dummyRootD );
   n = node->nextInT;
   dummyRootD->nextInT = n;
   n->prevInT = dummyRootD;

   for( arcDType *arc = arcsD ; arc != stopArcsD ; arc++ )
    if( ( ( arc->tail ) == node) || ( ( arc->head ) == node ) ) {
     arc->flow = 0;
     arc->ident = CLOSED;
     }

   for( arcDType *arc = dummyArcsD ; arc != stopDummyD ; arc++ )
    if( ( ( arc->tail ) == node ) || ( ( arc->head ) == node ) ) {
     arc->flow = 0;
     arc->ident = CLOSED;
     }

   CreateInitialDModifiedBalanceVector();
   PostDVisit( dummyRootD );
   ComputePotential( dummyRootD );
  #endif
  }
 }  // end( MCFSimplex::DelNode )

/*--------------------------------------------------------------------------*/

void MCFSimplex::OpenArc( cIndex name )
{
 if( name >= m )
  return;

 if( usePrimalSimplex ) {
  /* Quadratic case is not implemented for a theory bug.
     Infact a closed arc in the quadratic case has its cost fixed to infinity, and
     it's impossible to restore the old value. */

  #if( QUADRATICCOST == 0 )
   arcPType *arc = arcsP + name;
   if( arc->ident == CLOSED ) {
    arc->ident = AT_LOWER;
    arc->flow = 0; 
    }
  #endif
  }
 else {
  #if( QUADRATICCOST == 0 )
   arcDType *arc = arcsD + name;
   if( arc->ident == CLOSED ) {
    if( GTZ( ReductCost( arc ) , EpsCst ) )
     arc->ident = AT_LOWER;
    else {
     arc->ident = AT_UPPER;
     arc->flow = arc->upper;
     if( Senstv && ( status != kUnSolved ) ) {
      CreateInitialDModifiedBalanceVector();
      PostDVisit( dummyRootD );
      }
     else
      status = kUnSolved;
     }
    }
  #endif
  }
 }  // end( MCFSimplex:OpenArc )

/*--------------------------------------------------------------------------*/

MCFSimplex::Index MCFSimplex::AddNode( cFNumber aDfct )
{
 if( n >= nmax )
  return( Inf<Index>() );        

 n++;
 if( usePrimalSimplex ) {
  nodePType *newNode = nodesP + n - 1;
  stopArcsP->tail = newNode;
  stopArcsP->head = dummyRootP;
  stopArcsP->upper = Inf<FNumber>();
  stopArcsP->flow = 0;
  stopArcsP->cost = Inf<CNumber>();
  #if( QUADRATICCOST )
   stopArcsP->quadraticCost = 0;
  #else
   stopArcsP->ident = BASIC;
  #endif
  stopArcsP++;
  newNode->balance = aDfct;
  newNode->prevInT = dummyRootP;
  newNode->nextInT = dummyRootP->nextInT;
  (dummyRootP->nextInT)->prevInT = newNode;
  dummyRootP->nextInT = newNode;
  newNode->enteringTArc = stopArcsP--;
  newNode->potential = 0;
  #if(QUADRATICCOST)
   newNode->sumQuadratic = 0;
  #endif
  }
 else {
  #if( QUADRATICCOST == 0 )
   nodeDType *newNode = nodesD + n - 1;
   stopArcsD->tail = newNode;
   stopArcsD->head = dummyRootD;
   stopArcsD->upper = 0;
   stopArcsD->flow = 0;
   stopArcsD->cost = Inf<CNumber>();
   stopArcsD->ident = BASIC;
   newNode->balance = aDfct;
   newNode->prevInT = dummyRootD;
   newNode->nextInT = dummyRootD->nextInT;
   (dummyRootD->nextInT)->prevInT = newNode;
   dummyRootD->nextInT = newNode;
   newNode->enteringTArc = stopArcsD;
   newNode->potential = 0;
   newNode->firstFs = stopArcsD;
   newNode->firstBs = NULL;
   stopArcsD->nextFs = NULL;
   stopArcsD->nextBs = dummyRootD->firstBs;
   dummyRootD->firstBs = stopArcsD;
   stopArcsD++;
  #endif
  }

 return( n );

 }  // end( MCFSimplex::AddNode )

/*--------------------------------------------------------------------------*/

void MCFSimplex::ChangeArc( cIndex name , cIndex nSN , cIndex nEN )
{
 if( name >= m )
  return;

 CloseArc( name );

 if( usePrimalSimplex ) {
  if( nSN <= n )
   (arcsP + name)->tail = (nodesP + nSN + USENAME0 - 1);
  if( nEN <= n )
   (arcsP + name)->head = (nodesP + nEN + USENAME0 - 1);
  }
 else {
  if( nSN <= n )
   (arcsD + name)->tail = (nodesD + nSN + USENAME0 - 1);
  if( nEN <= n )
   (arcsD + name)->head = (nodesD + nEN + USENAME0 - 1);
  }

 OpenArc( name );

 }  // end( MCFSimplex::ChangeArc )

/*--------------------------------------------------------------------------*/

void MCFSimplex::DelArc( cIndex name )
{
 if( name >= m )
  return;

 if( usePrimalSimplex ) {
  arcPType *arc = arcsP + name;
  #if( QUADRATICCOST )
   if( arc->upper == -Inf<FNumber>() )
    return;

   if( arc->cost < Inf<CNumber>() )
  #else
   if( arc->cost == DELETED )
    return;

   if( arc->cost >= BASIC )
  #endif
    CloseArc( name );

  #if( QUADRATICCOST )
   arc->upper = -Inf<FNumber>();
  #else
   arc->ident = DELETED;
  #endif
  }
 else {
  #if( QUADRATICCOST == 0 )
   arcDType *arc = arcsD + name;
   if( arc->cost == DELETED )
    return;

   if( arc->cost >= BASIC )
    CloseArc( name );

   arc->ident = DELETED;
  #endif
  }
 }  // end( MCFSimplex::DelArc )

/*--------------------------------------------------------------------------*/

MCFSimplex::Index MCFSimplex::AddArc( cIndex Start , cIndex End ,
				      cFNumber aU , cCNumber aC ) 
{
 if( usePrimalSimplex ) {
  arcPType *arc = arcsP;
  #if( QUADRATICCOST )
   while( ( arc->upper != -Inf<FNumber>() ) && ( arc != stopArcsP ) )
  #else
   while( ( arc->ident > DELETED ) && ( arc != stopArcsP ) )
  #endif
    arc++;

  if( arc == stopArcsP ) {
   if( m >= mmax )
    return( Inf<Index>() );

   m++;
   stopArcsP++;
   }

  Index pos = ( arc - arcsP ) + 1;
  arc->tail = nodesP + Start + USENAME0 - 1;
  arc->head = nodesP + End + USENAME0 - 1;
  arc->upper = aU;
  arc->cost = aC;
  arc->flow = 0;
  #if( QUADRATICCOST )
   arc->quadraticCost = 0;
  #else
   arc->ident = AT_LOWER;
  #endif
  ComputePotential( dummyRootP );
  return( pos );
  }
 else {
  #if( QUADRATICCOST == 0 )
   arcDType *arc = arcsD;
   while( ( arc->ident > DELETED ) && ( arc != stopArcsD ) )
    arc++;

   if( arc == stopArcsD ) {
    if( m >= mmax )
     return( Inf<Index>() );

    m++;
    stopArcsD++;
    }

   Index pos = ( arc - arcsD ) + 1;
   arc->tail = nodesD + Start + USENAME0 - 1;
   arc->head = nodesD + End + USENAME0 - 1;
   arc->upper = aU;
   arc->cost = aC;
   if( GEZ( ReductCost( arc ) , EpsCst ) ) {
    arc->flow = 0;
    arc->ident = AT_LOWER;
    if( Senstv && ( status != kUnSolved ) ) {
     CreateInitialDModifiedBalanceVector();
     PostDVisit( dummyRootD );
     ComputePotential( dummyRootD );
     }
    else
     status = kUnSolved;
    }
   else {
    arc->flow = arc->upper;
    arc->ident = AT_UPPER;
    if( Senstv && ( status != kUnSolved ) ) {
     CreateInitialDModifiedBalanceVector();
     PostDVisit( dummyRootD );
     ComputePotential( dummyRootD );
     }
    else
     status = kUnSolved;
    }

   ComputePotential( dummyRootD );
   return( pos );
  #endif
  }
 }  // end( MCFSimplex::AddArc )

/*--------------------------------------------------------------------------*/

bool MCFSimplex::IsDeletedArc( cIndex name )
{
 if( name >= m )
  return( false );

 #if( QUADRATICCOST )
  return( ( ( arcsP + name )->upper == -Inf<FNumber>() ) );
 #else
  if( usePrimalSimplex )
   return( ( arcsP + name )->ident == DELETED );
  else
   return( ( arcsD + name )->ident == DELETED );
 #endif
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

MCFSimplex::~MCFSimplex()
{
 MemDeAllocCandidateList();
 MemDeAlloc(true);
 MemDeAlloc(false);
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void MCFSimplex::MemAlloc( void )
{
 if( usePrimalSimplex )        {
  nodesP = new nodePType[ nmax + 1 ];   // array of nodes
  arcsP = new arcPType[ mmax + nmax ];  // array of arcs
  dummyArcsP = arcsP + mmax;            // artificial arcs are in the last
                                        // nmax positions of the array arcs[]
  }
 else {
  nodesD = new nodeDType[ nmax + 1 ];   // array of nodes
  arcsD = new arcDType[ mmax + nmax ];  // array of arcs
  dummyArcsD = arcsD + mmax;            // artificial arcs are in the last nmax
                                        // positions of the array arcs[]
  }
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::MemDeAlloc( bool whatDeAlloc )
{
 if( whatDeAlloc ) {
  delete[] nodesP;
  delete[] arcsP;
  nodesP = NULL;
  arcsP = NULL;
  }
 else {
  delete[] nodesD;
  delete[] arcsD;
  nodesD = NULL;
  arcsD = NULL;
 }
 MemDeAllocCandidateList( );
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::MemAllocCandidateList( void )
{
 if( usePrimalSimplex ) {
  if( m < 10000 ) {
   numCandidateList = PRIMAL_LOW_NUM_CANDIDATE_LIST; 
   hotListSize = PRIMAL_LOW_HOT_LIST_SIZE;
   }
  else
   if( m > 100000 ) {
    numCandidateList = PRIMAL_HIGH_NUM_CANDIDATE_LIST; 
    hotListSize = PRIMAL_HIGH_HOT_LIST_SIZE ;
    }
   else {
    numCandidateList = PRIMAL_MEDIUM_NUM_CANDIDATE_LIST; 
    hotListSize = PRIMAL_MEDIUM_HOT_LIST_SIZE;
    }

  #if( QUADRATICCOST )
   int coef = 1;
   // If the number of the arcs is more than 10000, numCandidateList and hotListSize 
   // are increased to improve the performance of the Quadratic Primal Simplex
   if( m > 10000 )
    coef = 10;

   numCandidateList = numCandidateList * coef;
   hotListSize = hotListSize * coef;
  #endif

  if( forcedNumCandidateList > 0 )
   numCandidateList = forcedNumCandidateList;

  if( forcedHotListSize > 0 )
   hotListSize = forcedHotListSize;

  candP = new primalCandidType[ hotListSize + numCandidateList + 1 ];
  }
 else {
  if( n < 10000 ) {
   numCandidateList = DUAL_LOW_NUM_CANDIDATE_LIST;
   hotListSize = DUAL_LOW_HOT_LIST_SIZE;
   }
  else {
   numCandidateList = DUAL_HIGH_NUM_CANDIDATE_LIST;
   hotListSize = DUAL_HIGH_HOT_LIST_SIZE;
   }

  if( forcedNumCandidateList > 0 )
   numCandidateList = forcedNumCandidateList;

  if( forcedHotListSize > 0 )
   hotListSize = forcedHotListSize;

  candD = new dualCandidType[ hotListSize + numCandidateList + 1 ];
  }
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::MemDeAllocCandidateList( void )
{
 delete[] candP;
 candP = NULL;
 delete[] candD; 
 candD = NULL;
}

/*--------------------------------------------------------------------------*/

void MCFSimplex::CreateInitialPrimalBase( void )
{
 arcPType *arc;
 nodePType *node;
 for( arc = arcsP ; arc != stopArcsP ; arc++ ) {  // initialize real arcs
  arc->flow = 0;
  #if( QUADRATICCOST == 0 )
   arc->ident = AT_LOWER;
  #endif
  }

 for( arc = dummyArcsP ; arc != stopDummyP ; arc++ ) {  // initialize dummy arcs
  node = nodesP + ( arc - dummyArcsP );
  if( node->balance > 0 ) {  // sink nodes 
   arc->tail = dummyRootP;
   arc->head = node;
   arc->flow = node->balance;
   }
  else {  // source nodes or transit node
   arc->tail = node;
   arc->head = dummyRootP;
   arc->flow = -node->balance;
   }

  arc->cost = MAX_ART_COST;
  #if( QUADRATICCOST )
   arc->quadraticCost = 0; 
  #else
   arc->ident = BASIC;
  #endif
  arc->upper = Inf<FNumber>();
  }

 dummyRootP->balance = 0;
 dummyRootP->prevInT = NULL;
 dummyRootP->nextInT = nodesP;
 dummyRootP->enteringTArc = NULL;
 #if( QUADRATICCOST )
  dummyRootP->sumQuadratic = 0;
 #endif
 dummyRootP->potential = MAX_ART_COST;
 dummyRootP->subTreeLevel = 0;
 for( node = nodesP ; node != stopNodesP ; node++) {  // initialize nodes
  node->prevInT = node - 1;
  node->nextInT = node + 1;
  node->enteringTArc = dummyArcsP + (node - nodesP);
  #if( QUADRATICCOST )
   node->sumQuadratic = (node->enteringTArc)->quadraticCost;
  #endif
  if( node->balance > 0 )  // sink nodes
   node->potential = 2 * MAX_ART_COST;
  else                     // source nodes or transit node
   node->potential = 0;

  node->subTreeLevel = 1;
  }

 nodesP->prevInT = dummyRootP;
 ( nodesP + n - 1 )->nextInT = NULL;
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::CreateInitialDualBase( void )
{
 arcDType *arc;
 nodeDType *node;
 for( arc = dummyArcsD ; arc != stopDummyD ; arc++ ) {  // initialize dummy arcs
  node = nodesD + ( arc - dummyArcsD );
  arc->tail = node;
  arc->head = dummyRootD;
  arc->flow = -node->balance;
  arc->cost = MAX_ART_COST;
  #if( QUADRATICCOST )
   arc->quadraticCost = 0;
  #else
   arc->ident = BASIC;
  #endif
  arc->upper = 0;
  }

 for( arc = arcsD ; arc != stopArcsD ; arc++ ) {  // initialize real arcs
  if( GTZ( arc->cost , EpsCst ) ) {
   arc->flow = 0;
   #if( QUADRATICCOST == 0 )
    arc->ident = AT_LOWER;
   #endif
   }
  else {
   #if( QUADRATICCOST == 0 )
    arc->ident = AT_UPPER;
   #endif
   arc->flow = arc->upper;
   ( dummyArcsD + ( ( arc->tail ) - nodesD ) )->flow =
     ( dummyArcsD + ( ( arc->tail ) - nodesD ) )->flow - arc->upper;

   ( dummyArcsD + ( ( arc->head ) - nodesD ) )->flow =
    ( dummyArcsD + ( ( arc->head ) - nodesD ) )->flow + arc->upper;
   }
  }

 dummyRootD->balance = 0;
 dummyRootD->prevInT = NULL;
 dummyRootD->nextInT = nodesD;
 dummyRootD->enteringTArc = NULL;
 #if( QUADRATICCOST )
  dummyRootD->sumQuadratic = 0;
 #endif
 dummyRootD->potential = MAX_ART_COST;
 dummyRootD->subTreeLevel = 0;
 for( node = nodesD ; node != stopNodesD ; node++ ) {  // initialize nodes
  node->prevInT = node - 1;
  node->nextInT = node + 1;
  node->enteringTArc = dummyArcsD + ( node - nodesD );
  #if( QUADRATICCOST )
   node->sumQuadratic = ( node->enteringTArc )->quadraticCost;
  #endif
  node->potential = 0;
  node->subTreeLevel = 1;
  node->whenInT2 = 0;
  }

 nodesD->prevInT = dummyRootD;
 ( nodesD + n - 1 )->nextInT = NULL;
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::CreateAdditionalDualStructures( void )
{
 // this method creates, in a Dual context, the Backward Star and the
 // Forward Star of every node

 for( nodeDType *node = nodesD ; node != stopNodesD ; node++) {
  // initialize nodes
  node->firstBs = NULL;
  node->firstFs = NULL;
  node->numArcs = 0;
  }

 dummyRootD->firstBs = NULL;
 dummyRootD->firstFs = NULL;
 dummyRootD->numArcs = 0;
 for( arcDType *arc = arcsD ; arc != stopArcsD ; arc++ ) {
  // initialize real arcs
  arc->nextFs = ( arc->tail )->firstFs;
  ( arc->tail )->firstFs = arc;
  arc->nextBs = ( arc->head )->firstBs;
  ( arc->head )->firstBs = arc;
  ( arc->tail )->numArcs++;
  ( arc->head )->numArcs++;
  }

 ResetWhenInT2();
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::PrimalSimplex( void )
{
 #if( UNIPI_PRIMAL_INITIAL_SHOW == 0 )
  #if( UNIPI_PRIMAL_ITER_SHOW == 0 )
   #if( UNIPI_PRIMAL_FINAL_SHOW == 0 ) 
    cout << endl;
   #endif
  #endif
 #endif
 #if( UNIPI_PRIMAL_INITIAL_SHOW )
  cout << endl;
  for( int t = 0; t < 3; t++ )
   cout << "\t";
  cout << "PRIMALE MCFSimplex: ARCHI E NODI ALL' INIZIO" << endl;
  ShowSituation( 3 );
 #endif
 #if( QUADRATICCOST )
  #if( LIMITATEPRECISION )
   foValue = GetFO();
   int cont = 0;
  #endif
 #endif

 status = kUnSolved;
 if( pricingRule != kCandidateListPivot )
  arcToStartP = arcsP;

 iterator = 0;  // setting the initial arc for the Dantzig or First Elibigle Rule

 arcPType *enteringArc;
 arcPType *leavingArc;
 if( pricingRule == kCandidateListPivot )
  InitializePrimalCandidateList();

 while( status == kUnSolved ) {
  iterator++;
  switch( pricingRule ) {
  case( kDantzig ):
   enteringArc = RuleDantzig();
   break;
  case( kFirstEligibleArc ):
   enteringArc = PRuleFirstEligibleArc();
   break;
  case( kCandidateListPivot ):
   enteringArc = RulePrimalCandidateListPivot();
   break;
  }

 #if( QUADRATICCOST )
  #if( LIMITATEPRECISION )
   /* In the quadratic case with LIMITATEPRECISION == 1, the entering arcs
      are selected according to a thresold.
      This thresold is definited according to the old f.o. value.
      If Primal Simplex doesn't find an entering arc, it calculates again
      the f.o. value, and try again. */
   if( enteringArc == NULL ) {
    foValue = GetFO();
    switch( pricingRule ) {
    case( kDantzig ):
     enteringArc = RuleDantzig();
     break;
    case( kFirstEligibleArc ):
     enteringArc = PRuleFirstEligibleArc();
     break;
    case( kCandidateListPivot ):
     enteringArc = RulePrimalCandidateListPivot();
     break;
     }
    }
  #endif
 #endif

 if( pricingRule != kCandidateListPivot ) {
  // in every iteration the algorithm changes the initial arc for
  // Dantzig and First Eligible Rule.
  arcToStartP++;
  if( arcToStartP == stopArcsP )
   arcToStartP = arcsP;
  }

 if( enteringArc ) {
  arcPType *arc;
  nodePType *k1;
  nodePType *k2;
  /* If the reduced cost of the entering arc is > 0, the Primal Simplex
     pushes flow in the cycle determinated by T and entering arc for decreases 
     flow in the entering arc: in the linear case entering arc's flow goes to 0,
     in the quadratic case it decreases while it's possibile.
     If the reduced cost of the entering arc is < 0, the Primal Simplex pushes
     flow in the cycle  determinated by T and entering arc for increases flow
     in the entering arc: in the linear case entering arc's flow goes to upper
     bound, in the quadratic case it increases while it's possibile.  */

  #if( QUADRATICCOST )
   FONumber t;
   FONumber theta; 
   FONumber deltaFO;
   FNumber theta2;
   CNumber Q = ( enteringArc->tail )->sumQuadratic +
               ( enteringArc->head )->sumQuadratic + enteringArc->quadraticCost;
   // Q is the sum of the quadratic coefficient in the cycle determinated by T
   // and entering arc.
   FONumber rc = ReductCost( enteringArc );
   if( ETZ( Q, EpsCst ) )
    theta = Inf<FNumber>();  // This value will be certainly decreased 
   else
    theta = ABS( rc / Q );
    // This is the best theta value (with best f.o. value decrease) 

   leavingArc = enteringArc;
   nodePType *cycleRoot; // The radix of the cycle determinated by T and entering arc.
   if( GTZ( rc , EpsCst ) ) {
  #else
   FNumber t;
   FNumber theta;
   if( enteringArc->ident == AT_UPPER ) {
  #endif
    /*  Primal Simplex increases or decreases entering arc's flow.
	"theta" is a positive value.
	For this reason the algorithm uses two nodes ("k1" and "k2") to push
	flow ("theta") from "k1" to "k2".
	According to entering arc's reduct cost, the algorithm determinates
	"k1" and "k2" */

    k1 = enteringArc->head;
    k2 = enteringArc->tail;
    #if( QUADRATICCOST )
     theta = min( theta , enteringArc->flow );
     // The best value for theta is compared with the entering arc's bound
     theta2 = - theta;
    #else
     theta = enteringArc->flow;
    #endif
    }
   else {
    k1 = enteringArc->tail;
    k2 = enteringArc->head;        
    #if( QUADRATICCOST )
     theta = min( theta , enteringArc->upper - enteringArc->flow );
     // The best value for theta is compared with the entering arc's bound
     theta2 = theta;
    #else
     theta = enteringArc->upper - enteringArc->flow;
    #endif
    }

   nodePType *memK1 = k1;
   nodePType *memK2 = k2;
   leavingArc = NULL;
   #if( QUADRATICCOST )
    #if( LIMITATEPRECISION )
     deltaFO = rc * theta2 + Q * theta2 * theta2 / 2;
    #endif
    bool leavingReducesFlow = GTZ( rc , EpsCst );
   #else
    bool leavingReducesFlow = GTZ( ReductCost( enteringArc ) , EpsCst );
   #endif
   // Compute "theta", find outgoing arc and "root" of the cycle
   bool leave;
   // Actual "theta" is compared with the bounds of the other cycle's arcs
   while( k1 != k2 ) {
    if( k1->subTreeLevel > k2->subTreeLevel ) {
     arc = k1->enteringTArc;
     if( arc->tail != k1 ) {
      t = arc->upper - arc->flow;
      leave = false;
      }
     else {
      t = arc->flow;
      leave = true;
      }

     if( t < theta ) {
      // The algorithm has found a possible leaving arc
      theta = t;
      leavingArc = arc;
      leavingReducesFlow = leave;
      // If "leavingReducesFlow" == true, if this arc is selected to exit the base, 
      // it decreases its flow
      }

     k1 = Father( k1 , arc );
     }
    else {
     arc = k2->enteringTArc;
     if( arc->tail == k2 ) {
      t = arc->upper - arc->flow;
      leave = false;
      }
     else {
      t = arc->flow;
      leave = true;
      }

     if( t <= theta ) {
      // The algorithm has found a possible leaving arc
      theta = t;
      leavingArc = arc;
      leavingReducesFlow = leave;
      // If "leavingReducesFlow" == true, if this arc is selected to exit the base, 
      // it decreases its flow
      }

     k2 = Father(k2, arc);
     }
    }

   #if( QUADRATICCOST )
    cycleRoot = k1;
   #endif

   if( leavingArc == NULL )
    leavingArc = enteringArc;

   if( theta >= Inf<FNumber>() ) {
    status = kUnbounded;
    break;
    }

   // Update flow with "theta"
   k1 = memK1;
   k2 = memK2;
   #if( QUADRATICCOST )
    if( enteringArc->tail == k1 )
     theta2 = theta;
    else
     theta2 = -theta;

    // "theta" is a positive value in every case.
    // "theta2" is the real theta value according to the real
    // direction of the entering arc
    #if( LIMITATEPRECISION )
     deltaFO = rc * theta2 + Q * theta2 * theta2 / 2;
     // The decrease of the f.o. value in the quadratic case
    #endif
   #endif

     if( ! ETZ(theta , EpsFlw ) ) {
      if( enteringArc->tail == k1 )
       enteringArc->flow = enteringArc->flow + theta;
      else
       enteringArc->flow = enteringArc->flow - theta;

      while( k1 != k2 ) {
       if( k1->subTreeLevel > k2->subTreeLevel ) {
	arc = k1->enteringTArc;
	if( arc->tail != k1 )
	 arc->flow = arc->flow + theta;
	else
	 arc->flow = arc->flow - theta;

	k1 = Father(k1, k1->enteringTArc);
        }
       else {
	arc = k2->enteringTArc;
	if( arc->tail == k2 )
	 arc->flow = arc->flow + theta;
	else
	 arc->flow = arc->flow - theta;

	k2 = Father(k2, k2->enteringTArc);
        }
       }
      }

     if( enteringArc != leavingArc ) {
      bool leavingBringFlowInT2 = ( leavingReducesFlow == 
	( ( leavingArc->tail )->subTreeLevel > ( leavingArc->head )->subTreeLevel ) );
      // "leavingBringFlowInT2" == true if leaving arc brings flow to the subtree T2
      if( leavingBringFlowInT2 == ( memK1 == enteringArc->tail ) ) {
       k2 = enteringArc->tail;
       k1 = enteringArc->head;
       }
      else {
       k2 = enteringArc->head; 
       k1 = enteringArc->tail;
       }
      }

     #if( QUADRATICCOST == 0 )
      if( leavingReducesFlow )
       leavingArc->ident = AT_LOWER;
      else
       leavingArc->ident = AT_UPPER;

      if( leavingArc != enteringArc ) {
       enteringArc->ident = BASIC;
       nodePType *h1;
       nodePType *h2;
       // "h1" is the node in the leaving arc with smaller tree's level 
       if( ( leavingArc->tail )->subTreeLevel < ( leavingArc->head )->subTreeLevel ) {
	h1 = leavingArc->tail;
	h2 = leavingArc->head;
        }
       else {
	h1 = leavingArc->head;
	h2 = leavingArc->tail;
        }

       UpdateT(leavingArc, enteringArc, h1, h2, k1, k2);
       // Update potential of the subtree T2
       k2 = enteringArc->head;
       CNumber delta = ReductCost(enteringArc);
       if( ( enteringArc->tail )->subTreeLevel > ( enteringArc->head )->subTreeLevel ) {
	delta = -delta;
	k2 = enteringArc->tail;
        }

       AddPotential( k2 , delta );
       // In the linear case Primal Simplex only updates the potential of the nodes of
       // subtree T2
       }
     #else
      if( leavingArc != enteringArc ) {
       nodePType *h1;
       nodePType *h2;
       // "h1" is the node in the leaving arc with smaller tree's level 
       if( ( leavingArc->tail )->subTreeLevel <
	   ( leavingArc->head )->subTreeLevel ) {
	h1 = leavingArc->tail;
	h2 = leavingArc->head;
        }
       else {
	h1 = leavingArc->head;
	h2 = leavingArc->tail;
        } 

       // Update the basic tree
       UpdateT( leavingArc , enteringArc , h1 , h2 , k1 , k2 );
       }

      #if( OPTQUADRATIC )
       nodePType *h1;
       nodePType *h2;
       if( ( leavingArc->tail )->subTreeLevel <
	   ( leavingArc->head )->subTreeLevel ) {
	h1 = leavingArc->tail;
	h2 = leavingArc->head;
        }
       else {
	h1 = leavingArc->head;
	h2 = leavingArc->tail;
        }

       nodePType *node = h1;
       nodePType *updateNode = h1;
       if( h1 == cycleRoot )
	ComputePotential( cycleRoot );
       else {
	while( node != cycleRoot ) {
	 arcPType *entArc = node->enteringTArc;
	 if( ! ETZ( entArc->quadraticCost , EpsCst ) )
	  updateNode = node;

	 node = Father( node , entArc );
	 }

	ComputePotential( updateNode );
	node = h2;
	updateNode = h2;
	while( node != cycleRoot ) {
	 arcPType *entArc = node->enteringTArc;
	 if( ! ETZ( entArc->quadraticCost , EpsCst ) )
	  updateNode = node;

	 node = Father( node , entArc );
	 }

	ComputePotential( updateNode );
        }
      #else
       // Update the potential of the node "manually"
       ComputePotential( cycleRoot );
      #endif

      #if( LIMITATEPRECISION )
       cont = cont + 1;
       if( cont == recomputeFOLimits ) {
	cont = 0;
	foValue = GetFO();  // Calculate f.o. value manually
        }
       else
	foValue = foValue + deltaFO;
        // Calculate the f.o. value with the estimated decrease
      #endif
     #endif
    }
   else {
    status = kOK;
    // If one of dummy arcs has flow bigger than 0, the solution is unfeasible.
    for( arcPType *arc = dummyArcsP ; arc != stopDummyP ; arc++ )
     if( GTZ( arc->flow , EpsFlw ) ) 
      status = kUnfeasible;
    }

   if( ( status == kUnSolved ) && MaxTime && MCFt ) {
    double tu, ts;
    TimeMCF( tu , ts );
    if( MaxTime < tu + ts )
     status = kStopped;
    }

   if( ( status == kUnSolved ) && MaxIter)
    if( MaxIter < (int) iterator )
     status = kStopped;

   #if( UNIPI_PRIMAL_ITER_SHOW )
    int it = (int) iterator;
    if( it % UNIPI_PRIMAL_ITER_SHOW == 0 ) {        
     cout << endl;
     for( int t = 0; t < 3; t++ )
      cout << "\t";
     cout << "PRIMALE MCFSimplex: ARCHI E NODI ALLA " << it << " ITERAZIONE" << endl;
     ShowSituation( 3 );
     #if( FOSHOW )
      if( (int) iterator % FOSHOW == 0 ) 
       clog << "Iteration = " << iterator << " of = "
        #if( LIMITATEPRECISION && QUADRATICCOST )
	    << foValue
        #else
	    << GetFO()
        #endif
            << endl;
     #endif
     }
   #endif
  }

 #if( UNIPI_PRIMAL_FINAL_SHOW )
  cout << endl;
  for( int t = 0; t < 3; t++ )
   cout << "\t";
  cout << "PRIMALE UniPi: ARCHI E NODI ALLA FINE" << endl;
  ShowSituation( 3 );
 #endif

 }  // end( PrimalSimplex )

/*--------------------------------------------------------------------------*/

void MCFSimplex::DualSimplex( void )
{
 #if( UNIPI_PRIMAL_INITIAL_SHOW == 0 )
  #if( UNIPI_PRIMAL_ITER_SHOW == 0 )
   #if( UNIPI_PRIMAL_FINAL_SHOW == 0 ) 
    cout << endl;
   #endif
  #endif
 #endif
 #if( UNIPI_DUAL_INITIAL_SHOW )
  cout << endl;
  for( int t = 0; t < 3; t++ )
   cout << "\t";
  cout << "DUALE MCFSimplex: ARCHI E NODI ALL' INIZIO" << endl;
  ShowSituation( 3 );
 #endif

 if( pricingRule != kCandidateListPivot )
  arcToStartD = arcsD;

 iterator = 0;
 arcDType *enteringArc;
 arcDType *leavingArc;
 if( pricingRule == kCandidateListPivot )
  InitializeDualCandidateList();

 status = kUnSolved;
 while( status == kUnSolved ) {
  iterator++;
  if( iterator == Inf<iteratorType>() ) {
   ResetWhenInT2(); // Restore to 0 every nodes' field "whenInT2"
   iterator = 1;
   }

  switch( pricingRule ) {
  case( kDantzig ):
   leavingArc = DRuleFirstEligibleArc();
   break; // si esegue cmq FEA
  case( kFirstEligibleArc ):
   leavingArc = DRuleFirstEligibleArc();
   break;
  case( kCandidateListPivot ):
   leavingArc = RuleDualCandidateListPivot();
   break;
   }

  if( pricingRule != kCandidateListPivot ) {
   arcToStartD++;
   if( arcToStartD == stopArcsD )
    arcToStartD = dummyArcsD;

   if( arcToStartD == stopDummyD )
    arcToStartD = arcsD;
   // Setting the initial arc for the Dantzig or First Elibigle Rule
   }

  if( leavingArc ) {
   bool leavingArcInL = false;
   bool leavingArcFromT1toT2;
   if( LTZ( leavingArc->flow , EpsFlw ) )
    leavingArcInL = true;

   nodeDType *h1;
   nodeDType *h2;
   if( ( leavingArc->tail )->subTreeLevel <
       ( leavingArc->head )->subTreeLevel ) {
    h1 = leavingArc->tail;
    h2 = leavingArc->head;
    leavingArcFromT1toT2 = true;
    }
   else {
    h1 = leavingArc->head;
    h2 = leavingArc->tail;
    leavingArcFromT1toT2 = false;
    }

   Index numOfT2Arcs = 0;
   int level = h2->subTreeLevel;
   nodeDType *node = h2;
   node->whenInT2 = iterator;
   nodeDType *lastNodeOfT2 = h2;
   numOfT2Arcs = node->numArcs;
   while( node->nextInT && ( ( node->nextInT )->subTreeLevel > level ) ) {
    node = node->nextInT;
    lastNodeOfT2 = node;
    numOfT2Arcs = numOfT2Arcs + node->numArcs;
    node->whenInT2 = iterator;
    }

   /* The Dual Simplex has determinated the leaving arc, and so the
      subtrees T1 and T2. Dual Simplex scans T2 to fix the fields "whenInT2"
      of T2's nodes to the iteration value, and counts the entering and
      outgoing arcs from these nodes. According to this number, it decides
      to scan the Backward and Forward of the subtree (T1 or T2) with the
      minor number of entering/outgoing arcs from its nodes. */
   enteringArc = NULL;
   bool lv = ( leavingArcFromT1toT2 == leavingArcInL );
   CNumber maxRc = Inf<CNumber>();
   //Search arc in the Forward Star and Backward Star of nodes of T1
   if( numOfT2Arcs > m ) {
    // Dual Simplex starts from the node which follows the dummy root.
    node = dummyRootD->nextInT;
    bool fine = false;
    while( fine == false ) {
     /* If node is the root of subtree T2, Dual Simplex jumps to the node
	(if exists) which follows the last node of T2 */
     if( node == h2 )
      if( lastNodeOfT2->nextInT )
       node = lastNodeOfT2->nextInT;
      else
       break;

     // Search arc in the Backward Star of nodes of T1
     arcDType *arc = node->firstBs;
     while( arc ) {
      if( ( arc->tail )->whenInT2 == iterator ) {
       // Evaluate each arc from T2 to T1 which isn't in T
       if( arc->ident == AT_LOWER ) {
	if( lv ) {
	 CNumber rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
             search is ended: this is the arc which enters in T */
	  if( ETZ( rc , EpsCst) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }

       if( arc->ident == AT_UPPER ) {
	if( ! lv ) {
	 CNumber rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
             search is ended: this is the arc which enters in T */

	  if( ETZ( rc , EpsCst ) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }
       }

      arc = arc->nextBs;
      }

     // Search arc in the Forward Star of nodes of T1
     arc = node->firstFs;
     while( arc ) {
      if( ( arc->head )->whenInT2 == iterator ) {
       // Evaluate each arc from T1 to T2 which isn't in T
       if( arc->ident == AT_LOWER ) {
	if( ! lv ) {
	 CNumber rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
	     search is ended: this is the arc which enters in T */
	  if( ETZ( rc , EpsCst ) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }

       if( arc->ident == AT_UPPER ) {
	if( lv ) {
	 CNumber rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
	     search is ended: this is the arc which enters in T */
	  if( ETZ( rc , EpsCst ) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }
       }

      arc = arc->nextFs;
      }

     node = node->nextInT;
     if( node == NULL )
      fine = true;
     }
    }
   // Search arc in the Forward Star and Backward Star of nodes of T2
   else {
    node = h2;
    bool fine = false;
    while( fine == false ) {
     // Search arc in the Backward Star of nodes of T2
     arcDType *arc = node->firstBs;
     CNumber rc;
     while( arc ) {
      if( ( arc->tail )->whenInT2 != iterator ) {
       // Evaluate each arc from T1 to T2 which isn't in T
       if( arc->ident == AT_LOWER ) {
	if( ! lv ) {
	 rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
	     search is ended: this is the arc which enters in T */
	  if( ETZ( rc , EpsCst ) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }

       if( arc->ident == AT_UPPER ) {
	if( lv ) {
	 rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
	     search is ended: this is the arc which enters in T */
	  if( ETZ( rc , EpsCst ) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }
       }

      arc = arc->nextBs;
      }

     // Search arc in the Forward Star of nodes of T2
     arc = node->firstFs;
     while( arc ) {
      if( ( arc->head )->whenInT2 != iterator ) {
       // Evaluate each arc from T2 to T1 which isn't in T
       if( arc->ident == AT_LOWER ) {
	if( lv ) {
	 rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
	     search is ended: this is the arc which enters in T */
	  if( ETZ( rc , EpsCst ) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }

       if( arc->ident == AT_UPPER ) {
	if( ! lv ) {
	 rc = ABS( ReductCost( arc ) );
	 if( LT( rc , maxRc , EpsCst ) ) {
	  enteringArc = arc;
	  maxRc = rc;
	  /* If arc is appropriate to enter in T and its reduct cost is 0, 
	     search is ended: this is the arc which enters in T */
	  if( ETZ( rc , EpsCst) ) {
	   fine = true;
	   break;
	   }
	  }
	 }
        }
       }

      arc = arc->nextFs;
      }

     if( node == lastNodeOfT2 )
      fine = true;
     else
      node = node->nextInT;

     }
    }

   if( enteringArc ) {
    FNumber theta = -leavingArc->flow;
    if( GTZ( leavingArc->flow , EpsFlw ) )
     theta = leavingArc->flow - leavingArc->upper;
    // Initial value of theta is the infeasibility of the leaving arc

    FNumber t;
    nodeDType *k1;
    nodeDType *k2;
    /* if entering arc is in U, Dual Simplex pushs flow in the cycle 
       determinated by T and entering arc for decreases flow in the entering arc:
       if entering arc is in L, Dual Simplex pushs flow in the cycle 
       determinated by T and entering arc for increases flow in the entering arc:
       Dual Simplex increases or decreases entering arc's flow.
       theta is a positive value.
       For this reason the algorithm uses two nodes (k1 and k2) to push flow
       (theta) from k1 to k2. According to entering arc's reduct cost, the
       algorithm determinates k1 and k2 */

    if( enteringArc->ident == AT_UPPER ) {
     k1 = enteringArc->head;
     k2 = enteringArc->tail;
     }
    else {
     k1 = enteringArc->tail;
     k2 = enteringArc->head;
     }

    nodeDType *memK1 = k1;
    nodeDType *memK2 = k2;
    arcDType *arc;                                
    k1 = memK1;
    k2 = memK2;
    // Update the flow
    while( k1 != k2 ) {
     if( k1->subTreeLevel > k2->subTreeLevel ) {
      arc = k1->enteringTArc;
      if( arc->tail != k1 )
       arc->flow = arc->flow + theta;
      else
       arc->flow = arc->flow - theta;

      k1 = Father(k1, k1->enteringTArc);
      }
     else {
      arc = k2->enteringTArc;
      if( arc->tail == k2 )
       arc->flow = arc->flow + theta;
      else
       arc->flow = arc->flow - theta;

      k2 = Father( k2 , k2->enteringTArc );
      }
     }

    if(leavingArcInL )
     leavingArc->ident = AT_LOWER;
    else
     leavingArc->ident = AT_UPPER;

    bool leavingBringFlowInT2 = ( leavingArcInL == 
	 ( ( leavingArc->tail )->subTreeLevel >
	   ( leavingArc->head )->subTreeLevel ) );
    // leavingBringFlowInT2 == true if leaving arc brings flow to the subtree T2
    if( leavingBringFlowInT2 != ( memK1 == enteringArc->tail ) ) {
     k2 = enteringArc->tail;
     k1 = enteringArc->head;
     }
    else {
     k2 = enteringArc->head; 
     k1 = enteringArc->tail;
     }

    if( enteringArc->ident == AT_LOWER )
     enteringArc->flow = enteringArc->flow + theta;
    else
     enteringArc->flow = enteringArc->flow - theta;

    enteringArc->ident = BASIC;
    UpdateT( leavingArc , enteringArc , h1 , h2 , k1 , k2 );
    // update potential of the subtree T2
    k2 = enteringArc->head;
    CNumber delta = ReductCost( enteringArc );
    if( ( enteringArc->tail) ->subTreeLevel >
	( enteringArc->head )->subTreeLevel ) {
     delta = -delta;
     k2 = enteringArc->tail;
     }

    // Dual Simplex only updates the potential of the T2's nodes
    AddPotential( k2 , delta );
    }
   else
    status = kUnfeasible;
    /* If Dual Simplex finds a leaving arc but it doesn't find an entering arc,
       the algorithm stops. At this point Dual Simplex has an unfeasible primal
       solution. */
   }
  else {
   status = kOK;
   // If one of dummy arcs has flow different than 0, the solution is unfeasible.
   for( arcDType *arc = dummyArcsD ; arc != stopDummyD ; arc++ )
    if( ! ETZ( arc->flow , EpsFlw ) ) {
     status = kUnfeasible;
     break;
     }
   }

  if( ( status == kUnSolved ) && MaxTime ) {
   double tu, ts;
   TimeMCF( tu , ts );
   if( MaxTime < tu + ts )
    status = kStopped;
   }

  if( ( status == kUnSolved ) && MaxIter && MCFt )
   if( MaxIter < (int) iterator )
    status = kStopped;

  #if( UNIPI_DUAL_ITER_SHOW )
   if( (int) iterator % UNIPI_DUAL_ITER_SHOW == 0 ) {
    cout << endl;
    for( int t = 0; t < 3; t++ )
     cout << "\t";
    cout << "DUALE MCFSimplex: ARCHI E NODI ALLA " << iterator << " ITERAZIONE"
	 << endl;
    ShowSituation( 3 );
    #if( FOSHOW )
     cout << "of = " << GetFO() << endl;
    #endif
    }
  #endif
  }

 #if( UNIPI_DUAL_ITER_SHOW )
  int it = (int) iterator;
  if( it % UNIPI_DUAL_ITER_SHOW == 0 ) {        
   cout << endl;
   for( int t = 0; t < 3; t++ )
    cout << "\t";
   cout << "DUALE MCFSimplex: ARCHI E NODI ALLA " << iterator << " ITERAZIONE"
	<< endl;
   Showsituation( 3 );
   }
 #endif

 }  // end( DualSimplex )

/*--------------------------------------------------------------------------*/

template<class N, class A>
void MCFSimplex::UpdateT( A *h , A *k , N *h1 , N *h2 , N *k1 , N *k2 )
{
 /* In subtree T2 there is a path from node h2 (deepest node of the leaving
    arc h and root of T2) to node k2 (deepest node of the leaving arc h and
    coming root of T2). With the update of T, this path will be overturned:
    node k2 will become the root of T2...
    The subtree T2 must be reordered and the field "subTreeLevel", which
    represents the depth in T of every node, of every T2's nodes is changed.
    Variable delta represents the increase of "subTreeLevel" value for node
    k2 and its descendants: probably this value is a negative value. */

 int delta = (k1->subTreeLevel) + 1 - (k2->subTreeLevel);
 N *root = k2;
 N *dad;
 
 /*To reorder T2, the method analyses every nodes of the path h2->k2, starting
   from k2. For every node, it moves the node's subtree from its original
   position to a new appropriate position. In particular k2 and its subtree
   (k2 is the new root of T2, so the first nodes of the new T2) will be moved
   next to node k1 (new father of k2), the next subtree will be moved beside
   the last node of k2's subtree....
   "previousNode" represents the node where the new subtree will be moved
   beside in this iterative action. At the start "previousNode" is the node
   k1 (T2 will be at the right of k1). */

 N *previousNode = k1;
 N *lastNode;
 /* "arc1" is the entering arc in T (passed by parameters).
    For every analysed node of path h2->k2, the method changes
    "enteringTArc" but it must remember the old arc, which will be the
    "enteringTArc" of the next analysed node.
    At the start "arc1" is k (the new "enteringTArc" of k2). */

 A *arc1 = k;
 A *arc2;
 bool fine = false;
 while( fine == false ) {
  // If analysed node is h2, this is the last iteration
  if( root == h2 )
   fine = true;

  dad = Father( root , root->enteringTArc );
  // Cut the root's subtree from T and update the "subLevelTree" of its nodes
  lastNode = CutAndUpdateSubtree( root , delta );
  // Paste the root's subtree in the right position;
  PasteSubtree( root , lastNode , previousNode );
  // In the next iteration the subtree will be moved beside the last node of
  // the actual analysed subtree.
  previousNode = lastNode;
  /* The increase of the subtree in the next iteration is different from
     the actual increase: in particual the increase increases itself (+2 
     at every iteration). */
  delta = delta + 2; 
  /* A this point "enteringTArc" of actual root is stored in "arc2" and
     changed; then "arc1" and "root" are changed. */
  arc2 = root->enteringTArc;
  root->enteringTArc = arc1;
  arc1 = arc2;
  root = dad;
  } 
 }

/*--------------------------------------------------------------------------*/

template<class N>
N* MCFSimplex::CutAndUpdateSubtree( N *root , int delta )
{
 int level = root->subTreeLevel;
 N *node = root;
 // The root of this subtree is passed by parameters, the last node is searched.
 while ( ( node->nextInT ) && ( ( node->nextInT )->subTreeLevel > level ) ) {
  node = node->nextInT;
  // The "subTreeLevel" of every nodes of subtree is updated
  node->subTreeLevel = node->subTreeLevel + delta;
  }

 root->subTreeLevel = root->subTreeLevel + delta;
 /* The 2 neighbouring nodes of the subtree (the node at the left of the root
    and the node at the right of the last node) is joined. */

 if( root->prevInT )
  ( root->prevInT )->nextInT = node->nextInT;
 if( node->nextInT )
  ( node->nextInT )->prevInT = root->prevInT;

 return( node );  // the method returns the last node of the subtree
 }

/*--------------------------------------------------------------------------*/

template<class N>
void MCFSimplex::PasteSubtree( N *root , N *lastNode , N *previousNode )
{
 /* The method inserts subtree ("root" and "lastNode" are the extremity of the
    subtree) after "previousNode". The method starts to identify the next node
    of "previousNode" ("nextNode"), so it joins "root" with "previousNode" and
    "lastNode" with "nextNode" (if exists). */

 N *nextNode = previousNode->nextInT;
 root->prevInT = previousNode;
 previousNode->nextInT = root;
 lastNode->nextInT = nextNode;
 if( nextNode )
  nextNode->prevInT = lastNode;
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::arcPType* MCFSimplex::RuleDantzig( void )
{
 arcPType *arc = arcToStartP;
 arcPType *enteringArc = NULL;
 #if( QUADRATICCOST )
  /* In the quadratic case used type for reduct cost is FONumber.
     Value "lim" is the fixed thresold for the decrease of the f.o. value */
  FONumber lim = EpsOpt * foValue / n;
  FONumber RC;
  FONumber maxValue = 0;
 #else
  CNumber RC;
  CNumber maxValue = 0;
 #endif

 do {
  // The method analyses every arc 
  #if( QUADRATICCOST )
   RC = ReductCost( arc );
   FNumber theta;
   /* If reduct cost of arc is lower than 0, the flow of the arc must increase.
      If reduct cost of arc is bigger than 0, the flow of the arc must decrease.
      "theta" is the difference between lower (upper) bound and the actual flow.
      */

   if( LTZ( RC , EpsCst ) )
    theta = arc->upper - arc->flow;

   if( GTZ( RC , EpsCst ) ) 
    theta = -arc->flow;

   // If it's possible to increase (or decrease) the flow in this arc
   if( ! ETZ( theta , EpsFlw ) ) {
    /* "Q" is the sum of the quadratic coefficient of the arc belonging the
       T path from tail's arc to head's arc
       "Q" is always bigger than 0 or equals to 0.
       If "Q" > 0, the value - RC / Q is the increase (decrease) of the flow
       with the best decrease of f.o. value.
       - RC/ Q must be compare with "theta" to avoid that the best increase
       (decrease) of the flow violates the bounds of the arc.
       This confront determines "theta". */

    CNumber Q = ( arc->tail )->sumQuadratic + ( arc->head )->sumQuadratic +
                arc->quadraticCost;

    if( GTZ( Q , EpsCst ) )
     if( GTZ( theta , EpsFlw ) )
      theta = min( theta , - RC / Q );        
     else
      theta = max( theta , - RC / Q );        

    /* Calculate the estimate decrease of the f.o. value if this arc is
       selected and flow is increased (decreased) by "theta" */

    CNumber deltaFO = RC * theta + Q * theta * theta / 2;
    // if deltaFO < 0 this arc is appropriate; if deltaFO is lower than
    // old decrease value, arc is the best arc.

    if( deltaFO < maxValue ) {
     maxValue = deltaFO;
     enteringArc = arc;
     }
    }
  #else
   if( arc->ident > BASIC ) {
    RC = ReductCost( arc );

    if( ( LTZ( RC , EpsCst ) && ( arc->ident == AT_LOWER ) ) || 
	( GTZ( RC , EpsCst ) && ( arc->ident == AT_UPPER ) ) ) {

     if( RC < 0 )
      RC = -RC;

     if( RC > maxValue ) {
      maxValue = RC;
      enteringArc = arc;
      }
     }
    }
  #endif

  arc++;
  if( arc == stopArcsP )
   arc = arcsP;

  } while( arc != arcToStartP );

 #if( ( LIMITATEPRECISION ) && ( QUADRATICCOST ) )
  if( -maxValue <= lim )
   enteringArc = NULL;
 #endif

 return( enteringArc );
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::arcPType* MCFSimplex::PRuleFirstEligibleArc( void )
{
 arcPType *arc = arcToStartP;
 arcPType *enteringArc = NULL;
 #if( QUADRATICCOST )
  FONumber RC;
 #else
  CNumber RC;
 #endif

 do {
  #if( QUADRATICCOST )
   // In this method the "decrease f.o. value" criterion is not used.
   RC = ReductCost( arc );
   if( ( LTZ( RC , EpsCst ) && LT( arc->flow , arc->upper , EpsFlw ) ) || 
       ( GTZ( RC , EpsCst ) && GTZ( arc->flow , EpsFlw ) ) )
    enteringArc = arc;
  #else
   if( arc->ident > BASIC ) {
    RC = ReductCost( arc );
    if( ( LTZ( RC , EpsCst ) && ( arc->ident == AT_LOWER ) ) ||
	( GTZ( RC , EpsCst ) && ( arc->ident == AT_UPPER ) ) )
     enteringArc = arc;
    }
  #endif

   arc++;
   if( arc == stopArcsP )
    arc = dummyArcsP;
   if( arc == stopDummyP )
    arc = arcsP;

  } while( ( enteringArc == NULL ) && ( arc != arcToStartP ) );

 return( enteringArc );
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::arcDType* MCFSimplex::DRuleFirstEligibleArc( void )
{
 arcDType *arc = arcToStartD;
 arcDType *leavingArc = NULL;
 do {
  if( GT( arc->flow , arc->upper , EpsFlw ) || LTZ( arc->flow , EpsFlw ) )
   leavingArc = arc;

  arc++;
  if( arc == stopArcsD )
   arc = dummyArcsD;
  if( arc == stopDummyD )
   arc = arcsD;

  } while( ( leavingArc == NULL ) && ( arc != arcToStartD ) );

 return(leavingArc);
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::arcPType* MCFSimplex::RulePrimalCandidateListPivot( void )
{
 Index next = 0;
 Index i;
 Index minimeValue;
 if( hotListSize < tempCandidateListSize )
  minimeValue = hotListSize;
 else
  minimeValue = tempCandidateListSize;

 #if( QUADRATICCOST )
  // Check if the left arcs in the list continue to violate the dual condition
  for( i = 2 ; i <= minimeValue ; i++ ) {
   arcPType *arc = candP[ i ].arc;
   FONumber red_cost = ReductCost( arc );
   FNumber theta = 0;
   /* If reduct cost of arc is lower than 0, the flow of the arc must increase.
      If reduct cost of arc is bigger than 0, the flow of the arc must decrease.
      "theta" is the difference between lower (upper) bound and the actual flow.
      */

   if( LTZ( red_cost , EpsCst ) )
    theta = arc->upper - arc->flow;
   else
    if( GTZ( red_cost , EpsCst ) )
     theta = - arc->flow;

   // if it's possible to increase (or decrease) the flow in this arc
   if( theta != 0 ) {
    /* "Q" is the sum of the quadratic coefficient of the arc belonging the T
       path from tail's arc to head's arc
       "Q" is always bigger than 0 or equals to 0.
       If "Q" > 0, the value - RC / Q is the increase (decrease) of the flow
       with the best decrease of f.o. value.
       - RC/ Q must be compare with "theta" to avoid that the best increase
       (decrease) of the flow violates the bounds of the arc.
       This confront determines "theta". */

    CNumber Q = ( arc->tail )->sumQuadratic + ( arc->head )->sumQuadratic +
                arc->quadraticCost;

    if( GTZ( Q , EpsCst  ) )
     if( GTZ( theta , EpsFlw ) )
      theta = min( theta , - red_cost / Q );
     else
      theta = max( theta , - red_cost / Q );

    /* Calculate the estimate decrease of the f.o. value if this arc is
       selected and flow is increased (decreased) by "theta" */

    CNumber deltaFO = red_cost * theta + Q * theta * theta / 2;
    #if( LIMITATEPRECISION )
     // if deltaFO < 0 this arc is appropriate; if deltaFO is lower than
     // old decrease value, arc is the best arc.
     if( - deltaFO > ( EpsOpt * foValue / n ) ) {
    #else
     if( LTZ( deltaFO , EpsCst ) ) {
    #endif
      next++;
      candP[ next ].arc = arc;
      candP[ next ].absRC = -deltaFO;
      }
     }
    }

   tempCandidateListSize = next;
   Index oldGroupPos = groupPos;
   // Search other arcs to fill the list
   do {
    arcPType *arc;
    for( arc = arcsP + groupPos ; arc < stopArcsP ; arc += numGroup ) {
     FONumber red_cost = ReductCost( arc );
     FNumber theta = 0;
     /* If reduced cost is lower than 0, the flow of the arc must increase.
	If reduced cost is larger than 0, the flow of the arc must decrease.
	"theta" is the difference between lower (upper) bound and the actual
	flow. */

     if( LTZ( red_cost , EpsCst ) ) 
      theta = arc->upper - arc->flow;
     else
      if( GTZ( red_cost , EpsCst ) ) 
       theta = - arc->flow;

     // if it's possible to increase (or decrease) the flow in this arc
     if( theta != 0 ) {
      /* "Q" is the sum of the quadratic coefficient of the arc belonging the T
	 path from tail's arc to head's arc
	 "Q" is always bigger than 0 or equals to 0.
	 If "Q" > 0, the value - RC / Q is the increase (decrease) of the flow
	 with the best decrease of f.o. value.
	 - RC/ Q must be compare with "theta" to avoid that the best increase
	 (decrease) of the flow violates the bounds of the arc.
	 This confront determines "theta". */
      CNumber Q = ( arc->tail )->sumQuadratic + ( arc->head )->sumQuadratic +
                  arc->quadraticCost;

      if( GTZ( Q , EpsCst  ) )
       if( GTZ( theta , EpsFlw ) )
	theta = min( theta , - red_cost / Q );        
       else
	theta = max( theta , - red_cost / Q );        

      /* Calculate the estimate decrease of the f.o. value if this arc is
	 selected and flow is increased (decreased) by "theta" */

      CNumber deltaFO = red_cost * theta + Q * theta * theta / 2;
      #if( LIMITATEPRECISION )
       // if deltaFO < 0 this arc is appropriate; if deltaFO is lower than
       // old decrease value, arc is the best arc.
       if( -deltaFO > ( EpsOpt * foValue / n ) ) {
      #else
       if( LTZ( deltaFO , EpsCst ) ) {
      #endif
	tempCandidateListSize++;
	candP[ tempCandidateListSize ].arc = arc;
	candP[ tempCandidateListSize ].absRC = -deltaFO;
        }
       }
      }

     groupPos++;
     if( groupPos == numGroup )
      groupPos = 0;

     } while( ( tempCandidateListSize < hotListSize ) &&
	      ( groupPos != oldGroupPos ) );
 #else
  // Check if the left arcs in the list continue to violate the dual condition
  for( i = 2 ; i <= minimeValue ; i++ ) {
   arcPType *arc = candP[i].arc;
   CNumber red_cost = ReductCost( arc );

   if( ( LTZ( red_cost , EpsCst ) && ( arc->ident == AT_LOWER ) ) ||
       ( GTZ( red_cost , EpsCst ) && ( arc->ident == AT_UPPER ) ) ) {
    next++;
    candP[ next ].arc = arc;
    candP[ next ].absRC = ABS( red_cost );
    }
   }

  tempCandidateListSize = next;
  Index oldGroupPos = groupPos;
  // Search other arcs to fill the list
  do {
   arcPType *arc;
   for( arc = arcsP + groupPos ; arc < stopArcsP ; arc += numGroup ) {
    if( arc->ident == AT_LOWER ) {
     CNumber red_cost = ReductCost( arc );
     if( LTZ( red_cost , EpsCst ) ) {
      tempCandidateListSize++;
      candP[ tempCandidateListSize ].arc = arc;
      candP[ tempCandidateListSize ].absRC = ABS( red_cost );
      }
     }
    else
     if( arc->ident == AT_UPPER ) {
      CNumber red_cost = ReductCost( arc );
      if( GTZ( red_cost , EpsCst ) ) {
       tempCandidateListSize++;
       candP[ tempCandidateListSize ].arc = arc;
       candP[ tempCandidateListSize ].absRC = ABS( red_cost );
       }
      }
    }

   groupPos++;
   if( groupPos == numGroup )
    groupPos = 0;

   } while( ( tempCandidateListSize < hotListSize ) && ( groupPos != oldGroupPos ) );
 #endif

 if( tempCandidateListSize ) {
  SortPrimalCandidateList( 1 , tempCandidateListSize );
  return( candP[ 1 ].arc );
  }
 else
  return( NULL );
 }

/*--------------------------------------------------------------------------*/

inline void MCFSimplex::InitializePrimalCandidateList( void )
{
 numGroup = ( ( m - 1 ) / numCandidateList ) + 1;
 groupPos = 0;
 tempCandidateListSize = 0;
 }

/*--------------------------------------------------------------------------*/

inline void MCFSimplex::SortPrimalCandidateList( Index min , Index max )
{
 Index left = min;
 Index right = max;
 #if( QUADRATICCOST )
  FONumber cut = candP[ ( left + right ) / 2 ].absRC;
 #else
  CNumber cut = candP[ ( left + right ) / 2 ].absRC;
 #endif
 do {
  while( candP[ left ].absRC > cut) 
   left++;
  while( cut > candP[ right ].absRC) 
   right--;

  if( left < right )
   Swap( candP[ left ] , candP[ right ] );

  if(left <= right) {
   left++;
   right--;
   }
  } while( left <= right );

 if( min < right ) 
  SortPrimalCandidateList( min , right );
 if( ( left < max ) && ( left <= hotListSize ) ) 
  SortPrimalCandidateList( left , max );
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::arcDType* MCFSimplex::RuleDualCandidateListPivot( void )
{
 Index next = 0;
 // Check if the left arcs in the list continue to violate the primal condition
 for( Index i = 2 ; ( i <= hotListSize ) && ( i <= tempCandidateListSize ) ;
      i++ ) {
  nodeDType *node = candD[ i ].node;
  arcDType *arc = node->enteringTArc;
  cFNumber flow = arc->flow;
  if( LTZ( flow , EpsFlw ) ) {
   next++;
   candD[ next ].node = node;
   candD[ next ].absInfeas = ABS( flow );
   }

  if( GT( flow , arc->upper , EpsFlw ) ) {
   next++;
   candD[ next ].node = node;
   candD[ next ].absInfeas = flow - arc->upper;
   }
  }

 tempCandidateListSize = next;
 Index oldGroupPos = groupPos;
 // Search other arcs to fill the list
 do {
  nodeDType *node = nodesD + groupPos;
  for( node ; node < stopNodesD ; node += numGroup ) {
   arcDType *arc = node->enteringTArc;
   cFNumber flow = arc->flow;
   if( LTZ( flow , EpsFlw ) ) {
    tempCandidateListSize++;
    candD[ tempCandidateListSize ].node = node;
    candD[ tempCandidateListSize ].absInfeas = ABS( flow );
    }

   if( GT( flow , arc->upper , EpsFlw) ) {
    tempCandidateListSize++;
    candD[ tempCandidateListSize ].node = node;
    candD[ tempCandidateListSize ].absInfeas = flow - arc->upper;
    }
   }

  groupPos++;
  if( groupPos == numGroup )
   groupPos = 0;

  } while( ( tempCandidateListSize < hotListSize ) &&
	   ( groupPos != oldGroupPos ) );

 if( tempCandidateListSize ) {
  SortDualCandidateList( 1 , tempCandidateListSize );
  return( (candD[ 1 ].node)->enteringTArc );
  }
 else
  return( NULL );
 }

/*--------------------------------------------------------------------------*/

inline void MCFSimplex::InitializeDualCandidateList( void )
{
 numGroup = ( ( n - 1 ) / numCandidateList ) + 1;
 groupPos = 0;
 tempCandidateListSize = 0;
 }

/*--------------------------------------------------------------------------*/

inline void MCFSimplex::SortDualCandidateList(Index min, Index max)
{
 Index left = min;
 Index right = max;
 FNumber cut = candD[ ( left + right ) / 2 ].absInfeas;
 do {
  while( candD[ left ].absInfeas > cut )
   left++;
  while( cut > candD[ right ].absInfeas )
   right--;
  if( left < right ) 
   Swap( candD[left ] , candD[ right ] );
  if( left <= right) {
   left++;
   right--;
   }
  } while( left <= right );

 if( min < right )
  SortDualCandidateList( min , right );
 if( (left < max) && ( left <= hotListSize ) ) 
  SortDualCandidateList( left , max );
 }

/*--------------------------------------------------------------------------*/

template<class N, class RCT>
inline void MCFSimplex::AddPotential( N *r , RCT delta )
{
 int level = r->subTreeLevel;
 N *n = r;
 
 do {
  n->potential = n->potential + delta;
  n = n->nextInT;
  } while ( ( n ) && ( n->subTreeLevel > level ) );
 }

/*--------------------------------------------------------------------------*/

template<class N>
inline void MCFSimplex::ComputePotential( N *r )
{
 N *n = r;
 int level = r->subTreeLevel;
 FONumber cost;
 // If "n" is not the dummy root, the potential of "r" is computed.
 // If "n" is the dummy root, the potential of dummy root is a constant.
 
 do {
  if( n->enteringTArc ) {
   cost = ( n->enteringTArc )->cost;
   #if (QUADRATICCOST)
    // Also field "sumQuadratic" is updated
    n->sumQuadratic = ( Father( n , n->enteringTArc ) )->sumQuadratic +
                      ( n->enteringTArc )->quadraticCost;

    if( ! ETZ( ( n->enteringTArc )->flow , EpsFlw ) )
     cost = cost + ( ( n->enteringTArc )->quadraticCost * ( n->enteringTArc )->flow );
   #endif

   if( n == ( n->enteringTArc )->head ) 
    n->potential = ( Father( n , n->enteringTArc ) )->potential + cost;
   else
    n->potential = ( Father( n , n->enteringTArc ) )->potential - cost;
   }
  n = n->nextInT;
  } while( ( n ) && ( n->subTreeLevel > level ) );
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::CreateInitialPModifiedBalanceVector( void )
{
 int i = 0;
 delete[] modifiedBalance;
 modifiedBalance = new FNumber[ n ];
 // Initialited every node's modifiedBalance to his balance
 for ( nodePType *node = nodesP ; node != stopNodesP ; node++ ) {
  modifiedBalance[i] = node->balance;
  i++;
  }

 // Modify the vector according to the arcs out of base with flow non zero
 // Scan the real arcs

 for( arcPType *arc = arcsP ; arc != stopArcsP ; arc++ ) {
  #if( QUADRATICCOST )
   if( ( ! ETZ( arc->flow , EpsFlw ) ) &&
       ( ( arc->tail )->enteringTArc != arc ) && 
       ( ( arc->head )->enteringTArc != arc ) ) {
    i = (arc->tail) - nodesP;
    modifiedBalance[ i ] += arc->flow;
    i = (arc->head) - nodesP;
    modifiedBalance[ i ] -= arc->flow;
    }
  #else
   if( arc->ident == AT_UPPER ) {
    i = (arc->tail) - nodesP;
    modifiedBalance[ i ] += arc->upper;
    i = (arc->head) - nodesP;
    modifiedBalance[ i ] -= arc->upper;
    }
  #endif
  }

 // Scan the dummy arcs
 for( arcPType *arc = dummyArcsP ; arc != stopDummyP ; arc++ ) {
  #if( QUADRATICCOST )
   if ( ( ! ETZ( arc->flow , EpsFlw ) ) &&
	( ( arc->tail )->enteringTArc != arc ) && 
	( ( arc->head )->enteringTArc != arc ) ) {
    i = (arc->tail) - nodesP;
    modifiedBalance[ i ] += arc->flow;
    i = (arc->head) - nodesP;
    modifiedBalance[ i ] -= arc->flow;
    }
  #else
   if (arc->ident == AT_UPPER) {
    i = (arc->tail) - nodesP;
    modifiedBalance[ i ] += arc->upper;
    i = (arc->head) - nodesP;
    modifiedBalance[ i ] -= arc->upper;
   }
  #endif
  }
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::PostPVisit( nodePType *r )
{
 // The method controls if "r" is a leaf in T
 bool rLeaf = false;
 int i = r - nodesP;
 if( r->nextInT ) 
  if( ( r->nextInT )->subTreeLevel <= r->subTreeLevel )
   rLeaf = true;
  else
   rLeaf = true;

 if( rLeaf )  // If "r" is a leaf
  if( ( r->enteringTArc)->head == r ) // If enteringTArc of "r" goes in "r"
   ( r->enteringTArc )->flow = modifiedBalance[ i ];
  else // If enteringTArc of "r" goes out "r"
   ( r->enteringTArc )->flow = - modifiedBalance[ i ];
 else { // If "r" isn't a leaf
  nodePType *desc = r->nextInT;
  // Call PostPVisit for every child of "r"
  while( ( desc ) && ( desc->subTreeLevel > r->subTreeLevel ) ) {
   if( desc->subTreeLevel - 1 == r->subTreeLevel ) { // desc is a son of r
    PostPVisit( desc );

    if( ( desc->enteringTArc )->head == r ) // enteringTArc of desc goes in r
     modifiedBalance[ i ] -= ( desc->enteringTArc )->flow;
    else // If enteringTArc of "desc" goes out "r"
     modifiedBalance[ i ] += ( desc->enteringTArc )->flow;
    }
   desc = desc->nextInT;
   }

  if( r != dummyRootP )
   if( ( r->enteringTArc )->head == r ) // If enteringTArc of "r" goes in "r"
    ( r->enteringTArc )->flow = modifiedBalance[ i ];
   else // If enteringTArc of "r" goes out "r"
    ( r->enteringTArc )->flow = - modifiedBalance[ i ];
  }
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::BalanceFlow( nodePType *r )
{
 // used only by Primal Simplex to restore a primal feasible solution.
 if( r == dummyRootP ) {
  nodePType *node = dummyRootP->nextInT;
  while( node ) {
   // call this function recursively for every son of dummy root
   if( node->subTreeLevel == 1 ) 
    BalanceFlow( node );

   node = node->nextInT;
   }
  }
 else {
  // The method controls if "r" is a leaf in T
  bool rLeaf = false;
  if( r->nextInT )
   if( ( r->nextInT )->subTreeLevel <= r->subTreeLevel )
    rLeaf = true;
   else
    rLeaf = true;

  if( rLeaf )  // If "r" is a leaf
   AdjustFlow( r );  // The method controls if entering basic arc in "r" is
                     // not feasible; in case adjust its flow
  else { // If "r" isn't a leaf
   nodePType *node = r->nextInT;
   // Balance the flow of every child of "r"
   while ( ( node ) && ( node->subTreeLevel > r->subTreeLevel ) ) {
    if( node->subTreeLevel == r->subTreeLevel + 1 ) 
     BalanceFlow( node );
    node = node->nextInT;
    }

   // The method controls if entering basic arc in "r" is not feasible;
   //in case adjust its flow
   AdjustFlow( r );
   }
  }
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::AdjustFlow( nodePType *r )
{
 arcPType *arc = r->enteringTArc;
 if( arc >= dummyArcsP ) // If entering arc of "r" is a dummy arc
  if( LTZ( arc->flow , EpsFlw ) ) {
   // If this dummy arc has flow < 0, the algorithm overturns the arc
   nodePType *temp = arc->tail;
   arc->tail = arc->head;
   arc->head = temp;
   arc->flow = -arc->flow;
   }
 else {  // If entering arc of "r" is not a dummy arc
  bool orientationDown = ( arc->head == r );
  FNumber delta = 0;
  if( LTZ( arc->flow , EpsFlw ) ) { // If flow is < 0
   delta = -arc->flow;
   arc->flow = 0;
   #if( ! QUADRATICCOST )
    arc->ident = AT_LOWER;
   #endif
   }

  if( GT( arc->flow , arc->upper , EpsFlw ) ) {
   // If flow goes over the capacity of the arc
   delta = arc->upper - arc->flow;
   arc->flow = arc->upper;
   #if( ! QUADRATICCOST )
    arc->ident = AT_UPPER;
   #endif
   }

  /* This arc goes out from the basis, and the relative dummy arc goes in T.
     Then the algorithm push flow in the cycle made by the arc and some arcs
     of T to balance the flow. */

  if( ! ETZ( delta , EpsFlw ) ) {
   nodePType *node = Father( r , arc );
   while( node != dummyRootP ) {
    arc = node->enteringTArc;
    if( ( arc->head == node ) == orientationDown )
     arc->flow += delta;
    else
     arc->flow -= delta;

    node = Father( node , arc );
    }

   arcPType *dummy = dummyArcsP + ( r - nodesP );
   #if( ! QUADRATICCOST )
    dummy->ident = BASIC;
   #endif

   /* Update the structure of the tree. If entering basic arc of "r" is
      changed, subtree of "r"is moved next dummy root. */
   r->enteringTArc = dummy;
   int deltaLevel = 1 - r->subTreeLevel;
   nodePType *lastNode = CutAndUpdateSubtree( r , deltaLevel ); 
   PasteSubtree( r , lastNode , dummyRootP );
   if( ( dummy->head == r ) != orientationDown )
    dummy->flow += delta;
   else
    dummy->flow -= delta;

   if( LTZ( dummy->flow , EpsFlw ) ) {
    nodePType *temp = dummy->tail;
    dummy->tail = dummy->head;
    dummy->head = temp;
    dummy->flow = -dummy->flow;
    }
   }
  }
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::CreateInitialDModifiedBalanceVector( void )
{
 #if( ! QUADRATICCOST )
  int i = 0;
  modifiedBalance = new FNumber[ n ];
  // Initialited every node's modifiedBalance to his balance
  for( nodeDType *node = nodesD ; node != stopNodesD ; node++ ) {
   modifiedBalance[ i ] = node->balance;
   i++;
   }

  // Modify the vector according to the arcs out of base with flow non zero
  // Scan the real arcs
  for( arcDType *arc = arcsD ; arc != stopArcsD ; arc++ )
   if( arc->ident == AT_UPPER ) {
    i = (arc->tail) - nodesD;
    modifiedBalance[ i ] += arc->upper;
    i = (arc->head) - nodesD;
    modifiedBalance[ i ] -= arc->upper;
    }

  // Scan the dummy arcs
  for( arcDType *arc = dummyArcsD ; arc != stopDummyD ; arc++ )
   if( arc->ident == AT_UPPER ) {
    i = (arc->tail) - nodesD;
    modifiedBalance[ i ] += arc->upper;
    i = (arc->head) - nodesD;
    modifiedBalance[ i ] -= arc->upper;
    }
 #endif
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::PostDVisit( nodeDType *r )
{
 #if( ! QUADRATICCOST )
  // The method controls if "r" is a leaf in T
  bool rLeaf = false;
  int i = r - nodesD;
  if( r->nextInT )
   if( ( r->nextInT )->subTreeLevel <= r->subTreeLevel )
    rLeaf = true;
   else
    rLeaf = true;

  if( rLeaf ) // If "r" is a leaf
   if( ( r->enteringTArc)->head == r ) // If enteringTArc of "r" goes in "r"
    ( r->enteringTArc )->flow = modifiedBalance[ i ];
   else // If enteringTArc of "r" goes out "r"
    ( r->enteringTArc )->flow = - modifiedBalance[ i ];
  else { // If "r" isn't a leaf
   nodeDType *desc = r->nextInT;
   // Call PostDVisit for every child of "r"
   while( ( desc ) && ( desc->subTreeLevel > r->subTreeLevel ) ) {
    if( desc->subTreeLevel -1 == r->subTreeLevel ) { // desc is a son of r
     PostDVisit( desc );

     if( ( desc->enteringTArc )->head == r ) // enteringTArc of desc goes in r
      modifiedBalance[ i ] -= ( desc->enteringTArc )->flow;
     else // If enteringTArc of "desc" goes out "r"
      modifiedBalance[ i ] += ( desc->enteringTArc )->flow;
     }

    desc = desc->nextInT;
    }

   if( r != dummyRootD )
    if( ( r->enteringTArc )->head == r ) // If enteringTArc of "r" goes in "r"
     ( r->enteringTArc )->flow = modifiedBalance[ i ];
    else // If enteringTArc of "r" goes out "r"
     ( r->enteringTArc )->flow = - modifiedBalance[ i ];
   }
 #endif
 }

/*--------------------------------------------------------------------------*/

inline void MCFSimplex::ResetWhenInT2( void )
{
 for( nodeDType *n = nodesD ; n != stopNodesD ; n++)
  n->whenInT2 = 0;
 }

/*--------------------------------------------------------------------------*/

template<class N, class A>
inline N* MCFSimplex::Father( N *n , A *a )
{
 if( a == NULL )
  return NULL;

 if( a->tail == n )
  return( a->head );
 else
  return( a->tail );
 }

/*-------------------------------------------------------------------------*/

inline MCFSimplex::FONumber MCFSimplex::GetFO( void )
{
 FONumber fo = 0;
 if( usePrimalSimplex ) {
  arcPType *arco;
  for( arco = arcsP ; arco != stopArcsP ; arco++ ) {
   #if( QUADRATICCOST ) 
    if( ! ETZ( arco->flow , EpsFlw ) ) 
     fo += arco->flow * ( arco->cost + arco->flow * arco->quadraticCost / 2 );
   #else
    if( ( arco->ident == BASIC ) || ( arco->ident == AT_UPPER ) )
     fo += arco->cost * arco->flow;
   #endif
   }

  for( arco = dummyArcsP ; arco != stopDummyP ; arco++ ) {
   #if( QUADRATICCOST ) 
    if( ! ETZ( arco->flow , EpsFlw ) )
     fo += arco->flow * ( arco->cost + arco->flow * arco->quadraticCost / 2 );
   #else
    if( ( arco->ident == BASIC ) || ( arco->ident == AT_UPPER ) ) 
     fo += arco->cost * arco->flow;
   #endif
   }
  }
 else {
  arcDType *a;
  for( a = arcsD ; a != stopArcsD ; a++ ) {
   #if (QUADRATICCOST) 
    fo += ( a->cost * a->flow ) + ( a->quadraticCost * a->flow * a->flow ) / 2;
   #else
    if( ( a->ident == BASIC ) || (a->ident == AT_UPPER ) ) 
     fo += a->cost * a->flow;
   #endif
   }

  for( a = dummyArcsD ; a != stopDummyD ; a++) {
   #if (QUADRATICCOST) 
    fo += ( a->cost * a->flow ) + ( a->quadraticCost * a->flow * a->flow ) / 2;
   #else
    if( ( a->ident == BASIC ) || ( a->ident == AT_UPPER ) ) 
     fo += a->cost * a->flow;
   #endif
   }
  }

 return( fo );
 }

/*-------------------------------------------------------------------------*/

void MCFSimplex::PrintPNode( nodePType *nodo )
{
 if( nodo )
  if( nodo != dummyRootP )
   cout << ( nodo - nodesP + 1 );
  else
   cout << "r";
 else
  cout << "..";
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::PrintPArc( arcPType *arc )
{
 if( arc ) {
  cout << "(";
  PrintPNode( arc->tail );
  cout << ", ";
  PrintPNode( arc->head );
  cout << ")";
  }
 else
  cout << "..";
}

/*--------------------------------------------------------------------------*/

void MCFSimplex::PrintDNode( nodeDType *nodo )
{
 if( nodo )
  if( nodo != dummyRootD )
   cout << ( nodo - nodesD + 1 );
  else
   cout << "r";
 else
  cout << "..";
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::PrintDArc( arcDType *arc )
{
 if( arc ) {
  cout << "(";
  PrintDNode( arc->tail );
  cout << ", ";
  PrintDNode( arc->head );
  cout << ")";
  }
 else
  cout << "..";
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::nodePType* MCFSimplex::RecoverPNode( Index ind ) 
{
 if( ( ind < 0 ) || ( ind > n ) )
  return( NULL );
 if( ind )
  return( nodesP + ind - 1 );
 else
  return( dummyRootP );
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::arcPType* MCFSimplex::RecoverPArc( nodePType *tail ,
					       nodePType *head )
{
 if( ( tail == NULL ) || ( head == NULL ) )
  return( NULL );

 arcPType *arc = arcsP;
 while( ( arc->tail != tail ) || ( arc->head != head ) ) {
  arc++;
  if( arc == stopArcsP )
   arc = dummyArcsP;
   if( arc == stopDummyP )
    return( NULL );
  }

 return( arc );
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::nodeDType* MCFSimplex::RecoverDNode( Index ind )
{
 if( ( ind < 0 ) || ( ind > n ) ) 
  return( NULL );

 if( ind )
  return( nodesD + ind - 1 );
 else
  return( dummyRootD );
 }

/*--------------------------------------------------------------------------*/

MCFSimplex::arcDType* MCFSimplex::RecoverDArc( nodeDType *tail ,
					       nodeDType *head )
{
 if( ( tail == NULL ) || ( head == NULL ) ) 
  return( NULL );

 arcDType *arc = arcsD;
 while( ( arc->tail != tail ) || ( arc->head != head ) ) {
  arc++;
  if( arc == stopArcsD )
   arc = dummyArcsD;
  if( arc == stopDummyD )
   return( NULL );
  }

 return( arc );
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::infoPNode( nodePType *node , int tab )
{
 for( int t = 0 ; t < tab ; t++ )
  cout << "\t";
 cout << "Nodo ";
 PrintPNode( node );
 cout << ": b = " << node->balance << " y = " << node->potential << endl;
 #if( UNIPI_VIS_NODE_BASIC_ARC )
  cout << ": TArc=";
  PrintPArc( node->enteringTArc );
  cout << endl;
 #endif
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::infoPArc( arcPType *arc , int ind , int tab )
{
 for( int t = 0 ; t < tab ; t++ )
  cout << "\t";
 cout << "Arco ";
 PrintPArc( arc );
 cout << ": x = " << arc->flow;
 #if( UNIPI_VIS_ARC_UPPER )
  cout << " u = " << arc->upper;
 #endif
 #if( UNIPI_VIS_ARC_COST )
  cout << " c = " << arc->cost;
 #endif
 #if( QUADRATICCOST )
  #if( UNIPI_VIS_ARC_Q_COST )
   cout << " q = " << arc->quadraticCost;
  #endif
  cout << endl;
  for( int t = 0 ; t < tab ; t++ )
   cout << "\t";
  #if( UNIPI_VIS_ARC_REDUCT_COST )
   cout << " rc = " << MCFGetRC( ind );
  #endif
 #else
  cout << endl;
  for( int t = 0 ; t < tab ; t++ )
   cout << "\t";
  #if( UNIPI_VIS_ARC_REDUCT_COST )
   cout << " rc = " << MCFGetRC( ind );
  #endif
  #if( UNIPI_VIS_ARC_STATE )
   switch( arc->ident ) {
   case( BASIC ):    cout << " in T"; break;
   case( AT_LOWER ): cout << " in L"; break;
   case( AT_UPPER ): cout << " in U"; break;
   case( DELETED ):  cout << " canceled"; break;
   case( CLOSED ):   cout << " closed";
   }
  #endif
 #endif
 cout << endl;
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::infoDNode( nodeDType *node , int tab )
{
 for( int t = 0 ; t < tab; t++ )
  cout << "\t";
 cout << "Nodo ";
 PrintDNode( node );
 cout << ": b = " << node->balance << " y = " << node->potential;
 #if( UNIPI_VIS_NODE_BASIC_ARC )    
  cout << ": TArc=";
  PrintDArc( node->enteringTArc );
  cout << endl;
 #endif
 }

/*--------------------------------------------------------------------------*/

void MCFSimplex::infoDArc( arcDType *arc , int ind , int tab )
{
 for( int t = 0 ; t < tab ; t++ )
  cout << "\t";
 cout << "Arco ";
 PrintDArc( arc );
 cout << " x = " << arc->flow;
 #if( UNIPI_VIS_ARC_UPPER )
  cout << " u = " << arc->upper;
 #endif
 #if( UNIPI_VIS_ARC_COST )
  cout << " c = " << arc->cost;
 #endif
 #if( QUADRATICCOST )
  #if( UNIPI_VIS_ARC_Q_COST )
   cout << " q = " << arc->quadraticCost;
  #endif
  cout << endl;
  for( int t = 0 ; t < tab ; t++ )
   cout << "\t";
  #if( UNIPI_VIS_ARC_REDUCT_COST )
   cout << " rc = " << MCFGetRC( ind );
  #endif
 #else
  cout << endl;
  for( int t = 0 ; t < tab ; t++ )
   cout << "\t";
  #if( UNIPI_VIS_ARC_REDUCT_COST )
   cout << " rc = " << MCFGetRC( ind );
  #endif
  #if (UNIPI_VIS_ARC_STATE)
  switch( arc->ident ) {
  case( BASIC ):    cout << " in T"; break;
  case( AT_LOWER ): cout << " in L"; break;
  case( AT_UPPER ): cout << " in U"; break;
  case( DELETED ):  cout << " canceled"; break;
  case( CLOSED ):   cout << " closed";
  }
  #endif
 #endif
 cout << endl;
 } 

/*--------------------------------------------------------------------------*/

void MCFSimplex::ShowSituation( int tab )
{
 if( usePrimalSimplex ) {
  arcPType *arc;
  nodePType *node;
  int i = 0;
  for( arc = arcsP ; arc != stopArcsP ; arc++ ) {
   infoPArc( arc , i , tab );
   i++;
   }
  cout << endl;
  #if( UNIPI_VIS_DUMMY_ARCS )
   i = 0;
   for( arc = dummyArcsP ; arc != stopDummyP ; arc++ ) {
    infoPArc( arc , i , tab );                
    i++;
    }
   cout << endl;
  #endif
  infoPNode( dummyRootP , tab );
  for( node = nodesP ; node != stopNodesP ; node++ )
   infoPNode( node , tab );
  }
 else {
  arcDType *arc;
  nodeDType *node;
  int i = 0;
  for( arc = arcsD ; arc != stopArcsD ; arc++ ) {
   infoDArc( arc , i , tab );
   i++;
   }
  cout << endl;
  #if( UNIPI_VIS_DUMMY_ARCS )
   i = 0;
   for( arc = dummyArcsD ; arc != stopDummyD ; arc++) {
    infoDArc( arc , i , tab );
    i++;
    }
   cout << endl;
  #endif
  infoDNode( dummyRootD , tab );
  for( node = nodesD ; node != stopNodesD ; node++ )
   infoDNode( node , tab );
  }
 }

/*-------------------------------------------------------------------------*/
/*---------------------- End File MCFSimplex.C ----------------------------*/
/*-------------------------------------------------------------------------*/
