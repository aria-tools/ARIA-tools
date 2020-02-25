/*--------------------------------------------------------------------------*/
/*------------------------- File MCFCplex.C --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  Quadratic Min Cost Flow problems solver, based on calls to the Cplex--*/
/*--  Callable Libraries. Conforms to the standard (MCF) interface        --*/
/*--  defined in MCFClass.h.                                              --*/
/*--                                                                      --*/
/*--                            VERSION 1.40                              --*/
/*--                           22 - 09 - 2015                             --*/
/*--                                                                      --*/
/*--                  Original Idea and Implementation by:                --*/
/*--                                                                      --*/
/*--                          Antonio Frangioni                           --*/
/*--                            Antonio Manca                             --*/
/*--                           Matteo Sammartino                          --*/
/*--                                                                      --*/
/*--                       Operations Research Group                      --*/
/*--                      Dipartimento di Informatica                     --*/
/*--                         Universita' di Pisa                          --*/
/*--                                                                      --*/
/*--             Copyright (C) 1997 - 2015 by Antonio Frangioni.          --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFCplex.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( CPX_VERSION >= 800 )
 #define CPX_OPTIMAL CPX_STAT_OPTIMAL
 #define CPX_UNBOUNDED CPX_STAT_UNBOUNDED
 #define CPX_INFEASIBLE CPX_STAT_INFEASIBLE

 #define CPX_INForUNBD CPX_STAT_INForUNBD
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace MCFClass_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  These functions are not implemented as methods of the class, since  --*/
/*--  they don't use directly its data structures.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

#if( ( ! FNUMBER_IS_DOUBLE ) || ( ! CNUMBER_IS_DOUBLE ) )

template<class T>
static inline Index_Set SparseAssign( T* g1 , const double *g2 ,
				      Index_Set nms , Index n , Index Bs = 0 )
{
 // takes the dense n-vector g2 and assigns its nonzeroes to the sparse vector
 // g1, meanwhile constructiong its set of nonzeroes in nms, with names from
 // Bs onwards; returns a pointer to the first element in nms after the last
 // index vritten: this can be used for computing the number of nonzeroes in
 // the "sparsified" vector and/or for InINF-terminating nms[]

 for( ; n-- ; Bs++ ) {
  const double tg2 = *(g2++);
  if( tg2 ) {
   *(nms++) = Bs;
   *(g1++) = T( tg2 );
   }
  }

 return( nms );
 }

/*--------------------------------------------------------------------------*/

template<class T>
static inline Index_Set SparseAssign( T* g1 , const double *g2 ,
				      Index_Set nms , Index n ,
				      const double eps , Index Bs = 0 )
{
 // as SparseAssign(), but elements are considered nonzero only if they are
 // >= eps (the idea is that all elements are >= 0)

 for( ; n-- ; Bs++ ) {
  const double tg2 = *(g2++);
  if( tg2 >= eps ) {
   *(nms++) = Bs;
   *(g1++) = T( tg2 );
   }
  }

 return( nms );
 }

#endif

/*--------------------------------------------------------------------------*/

static inline MCFCplex::Index VectLength( const MCFCplex::cIndex_Set nms ,
					  MCFCplex::cIndex stp = 
					             Inf<MCFCplex::Index>() )
{
 // count the number of elements in the vector nms that are < stp; stops
 // as soon as the first element >= stp is found

 MCFCplex::cIndex_Set tnms = nms;
 while( *tnms < stp )
  tnms++;

 return( tnms - nms );
 }

/*--------------------------------------------------------------------------*/

static inline void VectFill( int *nms , int strt ,  MCFCplex::Index n )
{
 // fills nms with the indices strt, strt + 1, ..., strt + n - 1

 for( ; n-- ; )
  *(nms++) = strt++;
 }

/*--------------------------------------------------------------------------*/

template<class T>
static inline void VectAssign( T *const g , const T x , MCFCplex::cIndex n )
{
 // g[ i ] = x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) = x;
 }

/*--------------------------------------------------------------------------*/

template<class T>
static inline void VectAssign( T *const g1 , const T *const g2 ,
			       MCFCplex::cIndex n )
{
 // g1[ i ] = g2[ i ] for each i = 0 .. n - 1

 const T *tg2 = g2 + n;
 for( T *tg1 = g1 + n ; tg1 > g1 ; )
  *(--tg1) = *(--tg2);
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
static inline void VectMAssign( T1 *g1 , const T2 *g2 , MCFCplex::Index n )
{
 // g1 := - g2

 for( ; n-- ; )
  *(g1++) = - *(g2++);
 }

/*--------------------------------------------------------------------------*/

template<class T>
static inline MCFCplex::Index_Set Sparsify( T* g , MCFCplex::Index_Set B ,
				  MCFCplex::Index n , MCFCplex::Index Bs = 0 )
{
 // turns g from a "dense" n-vector to a "sparse" one, eliminating all items
 // that are exactly == 0; writes the set of nonzero items in B, with names
 // from Bs onwards, ordered in increasing sense; returns a pointer to the
 // first element in B after the last index vritten: this can be used for
 // computing the number of nonzeroes in the "sparsified" vector and/or for
 // InINF-terminating the set

 for( ; n ; n-- , g++ )
  if( *g )
   *(B++) = Bs++;
  else
   break;

 if( n ) {
  T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g ) {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( Sparsify )

/*--------------------------------------------------------------------------*/

template<class T>
static inline MCFCplex::Index_Set SparsifyT( T* g , MCFCplex::Index_Set B ,
					     MCFCplex::Index n , const T eps ,
					     MCFCplex::Index Bs = 0 )
{
 // as Sparsify(), but elements are considered nonzero only if they are >=
 // eps (the idea is that all elements are >= 0)

 for( ; n ; n-- , g++ )
  if( *g >= eps )
   *(B++) = Bs++;
  else
   break;

 if( n ) {
  T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g >= eps ) {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( SparsifyT )

/*--------------------------------------------------------------------------*/

template<class T>
static inline void VectSum( T *const g , const T x , MCFCplex::cIndex n )
{
 // g[ i ] += x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) += x;
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( USENAME0 )
 #define MINUSONE
 #define PLUSONE
#else
 #define MINUSONE -1
 #define PLUSONE +1
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFCplex --------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------PUBLIC METHODS------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------CONSTRUCTOR--------------------------------------*/
/*--------------------------------------------------------------------------*/

MCFCplex::MCFCplex( cIndex nmx , cIndex mmx , CPXENVptr extenv )
          :
          MCFClass( nmx , mmx )
{
 // setup environment, if necessary - - - - - - - - - - - - - - - - - - - - -

 if( extenv )
  env = extenv;
 else {
  if( ! EnvICntr++ ) {
   int ts;
   #if( CPX_VERSION >= 700 )
    genv = CPXopenCPLEX( &ts );
   #else
    genv = CPXopenCPLEXdevelop( &ts );
   #endif

   if( ( ! genv ) || ts )
    throw( MCFException( "Problem opening Cplex environment" ) );
   }

  env = genv;
  }

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( nmax && mmax )
  MemAlloc();
 else
  nmax = mmax = 0;

 // other initializations - - - - - - - - - - - - - - - - - - - - - - - - - -

 InstCntr++;

 }  // end( MCFCplex )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::LoadNet( cIndex nmx , cIndex mmx , cIndex pn , cIndex pm ,
			cFRow pU , cCRow pC , cFRow pDfct , cIndex_Set pSn ,
			cIndex_Set pEn )
{
 // allocating and deallocating memory- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( mmax && nmax ) && ( ( nmx != nmax ) || ( mmx != mmax ) ) ) {
  MemDeAlloc();
  nmax = mmax = 0;
  }

 if( ( mmx && nmx ) && ( ( nmx != nmax ) || ( mmx != mmax ) ) ) {
  nmax = nmx;
  mmax = mmx;
  MemAlloc();
  }

 if( ( ! nmax ) || ( ! mmax ) )  // just sit down in the corner and wait
  return;

 // now setting up data - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 n = pn;
 m = pm;

 // setup data structures for arc creation/deletion- - - - - - - - - - - - -

 #if( DYNMC_MCF_CPX )
  VectAssign( ArcPos , FNumber( 0 ) , FreePos = m );
  VectAssign( ArcPos + m , FNumber( Inf<FNumber>() ), mmax - m );
 #endif

 // create and set up temporary data structures - - - - - - - - - - - - - - -

 #if( ( ! INDEX_IS_UINT ) || ( ! USENAME0 ) )
  int* stn = new int[ m ];
  int* enn = new int[ m ];

  for( Index i = m ; i-- ; ) {
   stn[ i ] = pSn[ i ] MINUSONE;
   enn[ i ] = pEn[ i ] MINUSONE;
   }
 #else
  #define stn (int *) pSn
  #define enn (int *) pEn
 #endif

 double* sup = new double[ n ];

 if( pDfct )
  VectMAssign( sup , pDfct , n );  // invert the sign of deficits
 else
  VectAssign( sup , double( 0 ) , n );

 double* upc = new double[ m ];

 if( pU )
  for( Index i = m ; i-- ; )
   if( pU[ i ] == Inf<FNumber>() )
    upc[ i ] = CPX_INFBOUND;
   else
    upc[ i ] = double( pU[ i ] );
 else
  VectAssign( upc , CPX_INFBOUND , m );

 double* obj = new double[ m ];

 if( pC )
  for( Index i = m ; i-- ; )
   if( pC[ i ] == Inf<CNumber>() ) {
    #if( DYNMC_MCF_CPX )
     ArcPos[ i ] = pU[ i ];
    #endif
    obj[ i ] = upc[ i ] = 0;
    }
   else
    obj[ i ] = double( pC[ i ] );
 else
  VectAssign( obj , double( 0 ) , m );

 #if( DYNMC_MCF_CPX )
  while( FreePos && ( ArcPos[ FreePos - 1 ] == Inf<FNumber>() ) )
   FreePos--;
 #endif

 // load internal structure of Cplex- - - - - - - - - - - - - - - - - - - - -

 status = CPXNETcopynet( env , net , CPX_MIN , n , sup , NULL , m , stn ,
			 enn , NULL , upc , obj , NULL );
 if( status )
  throw( MCFException( "Problem loading data" ) );

 // setup QP data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 qp = NULL;             // problem is initially Network
 QPMthd = qpNSimplex;   // set default QP solving method to Network Simplex
 

 // delete temporaries- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] obj;
 delete[] upc;
 delete[] sup;

 #if( ( ! INDEX_IS_UINT ) || ( ! USENAME0 ) )
  delete[] enn;
  delete[] stn;
 #endif

 status = MCFClass::kUnSolved;

 }  // end( MCFCplex::LoadNet )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::SolveMCF( void )
{
 if( MCFt )
  MCFt->Start();
 
 if( net ) {
  CPXNETprimopt( env , net );  // call the network simplex- - - - - - - - - - 
  status = CPXNETgetstat( env , net );
  }
 else {
  CPXqpopt( env , qp );       // call the QP solver - - - - - - - - - - - - -
  status = CPXgetstat( env , qp );
  }

 switch( status ) {
   case( CPX_OPTIMAL ):       status = MCFClass::kOK;
                              break;
   case( CPX_INForUNBD ):
   case( CPX_INFEASIBLE ):    status = MCFClass::kUnfeasible;
                              break;
   case( CPX_UNBOUNDED ):     status = MCFClass::kUnbounded;
                              break;
   default:                   status = MCFClass::kStopped;
   }

 if( MCFt )
  MCFt->Stop();

 }  // end( MCFCplex::SolveMCF )

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR READING RESULTS  -------------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::MCFGetX( FRow F , Index_Set nms , cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 #if( FNUMBER_IS_DOUBLE )
  if( net ) 
   CPXNETgetx( env , net , F , strt , stp - 1 );
  else
   CPXgetx( env , qp , F , strt , stp - 1 );
   
  if( nms )
   #if( EPS_FLOW )
    *SparsifyT( F , nms , stp - strt , EpsFlw , strt ) = Inf<Index>();
   #else
    *Sparsify( F , nms , stp - strt , strt ) = Inf<Index>();
   #endif
 #else
  if( net ) 
   CPXNETgetx( env , net , F , strt , stp - 1 );
  else
   CPXgetx( env , qp , F , strt , stp - 1 );

  if( nms )
   #if( EPS_FLOW )
    *SparseAssign( F , val , nms , stp - strt , EpsFlw , strt ) = Inf<Index>();
   #else
    *SparseAssign( F , val , nms , stp - strt , strt ) = Inf<Index>();
   #endif
  else
   VectAssign( F , val , stp - strt );
 #endif

 }  // end( MCFCplex::MCFGetX( F ) )

/*--------------------------------------------------------------------------*/

void MCFCplex::MCFGetRC( CRow CR , cIndex_Set nms , cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( nms ) {
  if( net )
   CPXNETgetdj( env , net , val , strt , stp - 1 );
  else
   CPXgetdj( env , qp , val , strt , stp - 1 );

  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(CR++) = val[ h - strt ];
  }
 else {
  #if( CNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetdj( env , net , CR , strt , stp - 1 );
   else
    CPXgetdj( env , qp , CR , strt , stp - 1);
  #else
   if( net )
    CPXNETgetdj( env , net , val , strt , stp - 1 );
   else
    CPXgetdj( env , qp , val , strt , stp - 1);

   VectAssign( CR , val , stp - strt );
  #endif
  }
 }  // end( MCFCplex::MCFGetRC( CR ) )

/*--------------------------------------------------------------------------*/

MCFCplex::CNumber MCFCplex::MCFGetRC( cIndex i )
{
 double temp;
 if( net )
  CPXNETgetdj( env , net , &temp , int( i ) , int( i ) );
 else
  CPXgetdj( env , qp , &temp , int( i ) , int( i ) );

 return( CNumber( temp ) );

 }  // end( MCFCplex::MCFGetRC( i ) )

/*--------------------------------------------------------------------------*/

void MCFCplex::MCFGetPi( CRow P , cIndex_Set nms , cIndex strt , Index stp )
{
 if( stp > n )
  stp = n;

 if( nms ) {
  if( net )
   CPXNETgetpi( env , net , val , strt , stp - 1 ); 
  else
   CPXgetpi( env , qp , val , strt , stp - 1 ); 

  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(P++) = - val[ h - strt ];
  }
 else {
  #if( CNUMBER_IS_DOUBLE ) 
   if( net )
    CPXNETgetpi( env , net , P , strt , stp - 1 );          
   else
    CPXgetpi( env , qp , P , strt , stp - 1 );          

   for( Index i = 0 ; i < stp - strt ; i++ )
    P[ i ] = - P[ i ];
  #else
   if( net )	   
    CPXNETgetpi( env , net , val , strt , stp - 1 );
   else
    CPXgetpi( env , qp , val , strt , stp - 1 );

   VectMAssign( P , val , stp - strt );
  #endif
  }    
 }  // end( MCFGetPi( P ) )

/*--------------------------------------------------------------------------*/

MCFCplex::FONumber MCFCplex::MCFGetFO( void )
{
 if( status == MCFClass::kOK ) {
  double objval;
  if( net )
   CPXNETgetobjval( env , net , &objval );
  else
   CPXgetobjval(env, qp, &objval);

  return( FONumber( objval ) );  
  }
 else
  if( status == MCFClass::kUnbounded )
   return( - Inf<FONumber>() );
  else
   return( Inf<FONumber>() );

 }  // end( MCFCplex::MCFGetFO )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::MCFArcs( Index_Set Startv , Index_Set Endv ,
			cIndex_Set nms , cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; ) {
   int st , en;
   if( net ) 
    CPXNETgetarcnodes( env , net , &st , &en , int( i ) , int( i ) );
   else {
    st = Startn[ i ];
    en = Endn[ i ];
    }

   if( Startv )
    *(Startv++) = Index( st PLUSONE );
   if( Endv )
    *(Endv++) = Index( en PLUSONE );
   }
  }
 else {
  if( stp > m )
   stp = m;

  if( net ) {
   #if( INDEX_IS_UINT )
     CPXNETgetarcnodes( env , net , (int *) Startv , (int *) Endv , strt ,
	                stp - 1 );
   #else
    if( Startv )
     CPXNETgetarcnodes( env , net , ind , NULL , strt , stp - 1 );
    
    if( Endv )
     CPXNETgetarcnodes( env , net , NULL , ind , strt , stp - 1 );
   #endif
   }
  else {
   if( Startv )
    VectAssign( Startv , Index_Set( Startn + strt ) , stp - strt );
  
   if( Endv )
    VectAssign( Endv , Index_Set( Endn + strt ) , stp - strt );
   }

  #if( ! USENAME0 )
   if( Startv )
    VectSum( Startv , Index( 1 ) , stp - strt );

   if( Endv )
    VectSum( Endv , Index( 1 ) , stp - strt );
  #endif
  }
 }  // end( MCFArcs )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFCosts( CRow Costv , cIndex_Set nms ,
			 cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; ) {
   double cst;
   if( net )
    CPXNETgetobj( env , net , &cst , int( i ) , int( i ) );
   else
    CPXgetobj( env , qp , &cst , int( i ) , int( i ) );

   *(Costv++) = cst;
   }
  }
 else {
  if( stp > m )
   stp = m;

  #if( CNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetobj( env , net , Costv , strt , stp - 1 );
   else
    CPXgetobj( env , qp , Costv , strt , stp - 1 );
  #else
   if( net )
    CPXNETgetobj( env , net , val , strt , stp - 1 );
   else
    CPXgetobj( env , qp , val , strt , stp - 1 );

   VectAssign( Costv , val , stp - strt );
  #endif
  }
 }  // end( MCFCplex::MCFCosts )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFQCoef( CRow Qv , cIndex_Set nms , cIndex strt , Index stp ) 
{
 double qcoef = 0;
 if( nms ) { 
  for( Index i = 0 , arc ; i < m && ( arc = nms[ i ] ) < stp ; i++ )
   if( arc >= strt ) {
    if( qp )
     CPXgetqpcoef( env , qp , int( arc ) , int( arc ) , &qcoef );

    Qv[ i ] = CNumber( qcoef );
    }    
  }
 else {
  if( stp > m)
   stp = m;

  for( Index arc = strt , i = 0 ; arc < stp ; arc++ , i++ ) {
   if( qp )
    CPXgetqpcoef( env , qp , int( arc ) , int( arc ) , &qcoef );

   Qv[ i ] = CNumber( qcoef );
   }
  }
 } // end( MCFCplex::MCFQCoef )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFUCaps( FRow UCapv , cIndex_Set nms ,
			 cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; )
   #if( DYNMC_MCF_CPX )
    if( ArcPos[ i ] && ( ArcPos[ i ] < Inf<FNumber>() ) )
     *(UCapv++) = ArcPos[ i ];
    else
   #endif
    {
     double cap;
     if( net )
      CPXNETgetub( env , net , &cap , int( i ) , int( i ) );
     else
      CPXgetub( env , qp , &cap , int( i ) , int( i ) );

     *(UCapv++) = cap;
     }
  }
 else {
  if( stp > m )
   stp = m;

  #if( FNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetub( env , net , UCapv , strt , stp - 1 );  
   else
    CPXgetub( env , qp , UCapv , strt , stp - 1 );
  #else
   if( net )
    CPXNETgetub( env , net , val , strt , stp - 1 );
   else
    CPXgetub( env , qp , val , strt , stp - 1 );

   VectAssign( UCapv , val , stp - strt );
  #endif

  #if( DYNMC_MCF_CPX )
   for( Index i = strt ; i < stp ; i++ )
    if( ArcPos[ i ] && ( ArcPos[ i ] < Inf<FNumber>() ) )
     UCapv[ i ] = ArcPos[ i ];
  #endif
  }
 }  // end( MCFCplex::MCFUCaps )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFDfcts( FRow Dfctv , cIndex_Set nms ,
			 cIndex strt , Index stp )
{
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; ) {
   double dfct;
   if( net )
    CPXNETgetsupply( env , net , &dfct , int( i ) , int( i ) );
   else
    CPXgetrhs( env , qp , &dfct , int( i ) , int( i ) );
   
   *(Dfctv++) = - dfct;
   }
  }
 else {
  if( stp > n )
   stp = n;

  #if( FNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetsupply( env , net , Dfctv , strt , stp - 1 ); 
   else
    CPXgetrhs( env , qp  , Dfctv , strt , stp - 1 ); 

   for( FRow tDf = Dfctv + stp - strt ; tDf-- > Dfctv ; )
    *tDf = - *tDf;
  #else
   if( net )
    CPXNETgetsupply( env , net , val , strt , stp - 1 );
   else
    CPXgetrhs( env , qp , val , strt , stp - 1 );

   VectMAssign( Dfctv , val , stp - strt );
  #endif
  }
 }  // end( MCFCplex::MCFDfcts )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::ChgCosts( cCRow NCost , cIndex_Set nms ,
			 cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 Index cnt;
 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NCost++;
   }

  cnt = VectLength( nms , stp );
  }
 else
  cnt = stp - strt;

 #undef VALUEPASSED
 #undef INDEXPASSED

 #if( CNUMBER_IS_DOUBLE )
  #define VALUEPASSED CRow( NCost )
 #else
  VectAssign( val , NCost , cnt );
  #define VALUEPASSED val
 #endif

 if( nms ) {
  #if( INDEX_IS_UINT )
   #define INDEXPASSED (int *) nms
  #else
   VectAssign( ind , nms , cnt );
   #define INDEXPASSED ind
  #endif

  if( net )
   CPXNETchgobj( env , net , cnt , INDEXPASSED , VALUEPASSED );
  else
   CPXchgobj( env , qp , cnt , INDEXPASSED , VALUEPASSED );
  }
 else {
  VectFill( ind , strt , cnt );
  
  if( net )
   CPXNETchgobj( env , net , cnt , ind , VALUEPASSED );
  else
   CPXchgobj( env , qp , cnt , ind , VALUEPASSED );
  }
 }  // end( MCFCplex::ChgCosts )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgCost( Index arc , cCNumber NCost )
{
 double cost = double( NCost );
 int which = int( arc );

 if( net )
  CPXNETchgobj( env , net , 1 , &which , &cost );
 else
  CPXchgobj( env , qp  , 1 , &which , &cost);

 }  // end( MCFCplex::ChgCost )


/*--------------------------------------------------------------------------*/

void MCFCplex::ChgQCoef( cCRow NQCoef, cIndex_Set nms,
			 cIndex strt , Index stp ) 
{
 if( ! qp )   // the problem in not QP:
  TurnToQP(); // turn it into QP

 double qcoef = 0;
 if( nms ) {
  for( Index i = 0 , arc ; i < m && ( arc = nms[ i ] ) < stp ; i++ )
   if( arc >= strt ) 
   #if( DYNMC_MCF_CPX )
    if( ! ArcPos[ arc ] )
   #endif
   {
    if( NQCoef )
     qcoef = double( NQCoef[ i ] );

    CPXchgqpcoef( env , qp , int( arc ) , int( arc ) , qcoef );
    }
  }
 else {
  if( stp > m )
   stp = m;
 
  for( Index arc = strt , i = 0 ; arc < stp ; arc++ , i++ )
   #if( DYNMC_MCF_CPX )
    if( ! ArcPos[ arc ] )
   #endif
   {
    if( NQCoef )
     qcoef = double( NQCoef[ i ] );

    CPXchgqpcoef( env , qp , int( arc ) , int( arc ) , qcoef );
    }
  }
 }  // end( MCFCplex::ChgQCoef )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgQCoef( Index arc , cCNumber NQCoef ) 
{
 if( ! qp )
  TurnToQP();

 #if( DYNMC_MCF_CPX )
  if( ! ArcPos[ arc ] )
 #endif
   CPXchgqpcoef( env , qp , int( arc ) , int( arc ) , NQCoef );

} // end( MCFCplex::ChgQCoef )

/*--------------------------------------------------------------------------*/
   
void MCFCplex::ChgDfcts( cFRow NDfct , cIndex_Set nms ,
			 cIndex strt , Index stp ) 
{
 if( stp > n )
  stp = n;

 Index cnt;
 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NDfct++;
   }

  cnt = VectLength( nms , stp );
  }
 else
  cnt = stp - strt;

 #undef INDEXPASSED

 VectMAssign( val , NDfct , cnt );

 if( nms ) {
  #if( INDEX_IS_UINT )
   #define INDEXPASSED (int *) nms
  #else
   VectAssign( ind , nms , cnt );
   #define INDEXPASSED ind
  #endif

  if( net )
   CPXNETchgsupply( env , net , cnt , INDEXPASSED , val );
  else
   CPXchgrhs( env , qp , cnt , INDEXPASSED , val );
  }
 else {
  VectFill( ind , strt , cnt );

  if( net )
   CPXNETchgsupply( env , net , cnt , ind , val );
  else
   CPXchgrhs( env , qp , cnt , ind , val );
  }
 }  // end( MCFCplex::ChgDfcts )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgDfct( Index node , cFNumber NDfct )
{
 double dfct = - double( NDfct );
 int which = int( node );

 if( net ) 
  CPXNETchgsupply( env , net , 1 , &which , &dfct );
 else
  CPXchgrhs( env , qp , 1 , &which , &dfct );

 }  // end( MCFCplex::ChgDfct )

/*--------------------------------------------------------------------------*/
 
void MCFCplex::ChgUCaps( cFRow NCap ,  cIndex_Set nms ,
			 cIndex strt , Index stp )
{
 if( stp > m )
  stp = m;

 Index cnt;
 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NCap++;
   }

  cnt = VectLength( nms , stp );
  }
 else
  cnt = stp - strt;

 #undef VALUEPASSED
 #undef INDEXPASSED

 #if( FNUMBER_IS_DOUBLE )
  #define VALUEPASSED FRow( NCap )
 #else 
  VectAssign( val , NCap , cnt );
  #define VALUEPASSED val
 #endif

 if( nms ) {
  #if( INDEX_IS_UINT )
   #define INDEXPASSED (int *) nms
  #else 
   VectAssign( ind , nms , cnt );
   #define INDEXPASSED ind
  #endif

  if( net )
   CPXNETchgbds( env , net , cnt , INDEXPASSED , ChangeUB , VALUEPASSED );
  else
   CPXchgbds( env , qp , cnt , INDEXPASSED , ChangeUB , VALUEPASSED );
  }
 else {
  VectFill( ind , strt , cnt );

  if( net )
   CPXNETchgbds( env , net , cnt , ind , ChangeUB , VALUEPASSED );
  else
   CPXchgbds( env , qp , cnt , ind , ChangeUB , VALUEPASSED );
  }
 }  // end( MCFCplex::ChgUCaps )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgUCap( Index arc , cFNumber NCap )
{  
 double cap = double( NCap );
 int which = int( arc );

 if( net )
  CPXNETchgbds( env , net , 1 , &which , "U", &cap );
 else
  CPXchgbds( env , qp , 1 , &which , "U" , &cap );

 }  // end( MCFCplex::ChgUCap )

/*--------------------------------------------------------------------------*/
/*----------------- MODIFYING THE STRUCTURE OF THE GRAPH -------------------*/
/*--------------------------------------------------------------------------*/ 

void MCFCplex::CloseArc( cIndex name )
{
 double ub;
 int temp = int( name );

 #if( DYNMC_MCF_CPX )
  if( net )
   CPXNETgetub( env , net , &ub , temp , temp );
  else
   CPXgetub( env , qp , &ub , temp , temp );

  ArcPos[ name ] = FNumber( ub );  // save current upper bound
 #endif

 // set upper bound to 0
 ub = 0;
 if( net )
  CPXNETchgbds( env , net , 1 , &temp , "U" , &ub );
 else
  CPXchgbds( env , qp , 1 , &temp , "U" , &ub );

 }  // end( MCFCplex::CloseArc )

/*--------------------------------------------------------------------------*/

void MCFCplex::DelNode( cIndex name )
{
 int index = int( name MINUSONE );
 double newdfct = 0;
 int arccount_p = 0;

 if( net ) { 
  int arcbeg = 0 , surplus = 0;
  CPXNETchgsupply( env , net , 1 , &index , &newdfct );  // set supply to 0
  CPXNETgetnodearcs( env , net , &arccount_p , &arcbeg , ind , m , &surplus ,
		    index , index );  
  }
 else {
  for( Index arc = 0 ; arc < m ; arc++ )     // get incident arcs
   if( Startn[ arc ] == index || Endn[ arc ] == index ) 
    ind[ arccount_p++ ] = arc;
    
  CPXchgrhs( env , qp , 1 , &index , &newdfct );  // set rhs to 0
  }

 for( int i = 0 ; i < arccount_p ; )
  CloseArc( Index( ind[ i++ ] ) );  // incident arcs are "closed"
  
 }  // end( MCFCplex::DelNode )

/*--------------------------------------------------------------------------*/

void MCFCplex::OpenArc( cIndex name )
{
 #if( DYNMC_MCF_CPX )
  int temp = int( name );
  double ub = double( ArcPos[ name ] );
  ArcPos[ name ] = 0;

  // restore upper bound
  if( net )
   CPXNETchgbds( env , net , 1 , &temp , "U" , &ub ); 
  else
   CPXchgbds( env , qp , 1 , &temp , "U" , &ub ); 

 #else
  throw(
   MCFException( "MCFCplex::ChangeArc() not implemented if DYNMC_MCF_CPX == 0"
		 ) );
 #endif

 }  // end( MCFCplex::OpenArc )
 
/*--------------------------------------------------------------------------*/

MCFCplex::Index MCFCplex::AddNode( cFNumber aDfct )
{
 if( n < nmax ) {
  double dfct = double( aDfct );
 
  if( net ) 
   CPXNETaddnodes( env , net , 1 , &dfct , NULL );
  else
   CPXnewrows( env, qp , 1 , &dfct , "E" , NULL , NULL ); // add new row in
                                                          // constraint matrix 
  n++;
  }

 return( n );

 }  // end( MCFCplex::AddNode )

/*--------------------------------------------------------------------------*/

void MCFCplex::DelArc( cIndex name ) 
{
 #if( DYNMC_MCF_CPX )
  ArcPos[ name ] = Inf<FNumber>();  // position is available for a new arc
  if( name < FreePos )
   FreePos = name; 

  int which = int( name );

  if( name == m - 1 ) {  // deleting the last arc(s) - - - - - - - - - - - -
   do {
    m--;               // decrement number of arcs
 
    if( net )
     CPXNETdelarcs( env , net , which , which );
    else
     CPXdelcols( env , qp , which , which );  // delete matching column in
                                              // constraint matrix

    } while( ArcPos[ m - 1 ] == Inf<FNumber>() );

   if( FreePos > m ) 
    FreePos = m; 
   }
  else {               // just set the bound to 0- - - - - - - - - - - - - -
   double cpct = 0; 
 
   if( net )
    CPXNETchgbds( env , net , 1 , &which , "U" , &cpct );
   else
    CPXchgbds( env , qp , 1 , &which , "U" , &cpct );
   }
 #else
  throw(
   MCFException( "MCFCplex::DelArc() not implemented if DYNMC_MCF_CPX == 0"
	 	 ) );
 #endif

 }  // end( MCFCplex::DelArc )

/*-------------------------------------------------------------------------*/

MCFCplex::Index MCFCplex::AddArc( cIndex Start , cIndex End , cFNumber aU ,
				  cCNumber aC )
{
 #if( DYNMC_MCF_CPX )
  if( m >= mmax )
   return( Inf<Index>() );

  int strn = int( Start MINUSONE );
  int endn = int( End MINUSONE );
 
  double cst = double( aC );
  double cpct = double( aU );

  if( FreePos >= m ) {  // there is no arc marked as deleted
   ArcPos[ FreePos ] = 0;
   FreePos = ++m;     // increase number of arcs

   int NewArc = int( m - 1 );
   
   if( net ) 
    CPXNETaddarcs( env , net , 1 , &strn , &endn , NULL , &cpct , &cst ,
		   NULL );
   else {
    // add new column in objective function and constraint matrix
    CPXnewcols( env , qp , 1 , &cst , NULL, &cpct , NULL , NULL );
    
    // modify arc-node incidence matrix
    CPXchgcoef( env , qp , strn , NewArc , 1 );   
    CPXchgcoef( env , qp , endn , NewArc , -1 );

    Startn[ NewArc ] = strn;
    Endn[ NewArc ]   = endn;
    }

   return( Index( NewArc ) );
   }
  else {
   ArcPos[ FreePos ] = 0;
   int temp = int( FreePos );

   if( net ) {
    CPXNETchgarcnodes( env , net , 1 , &temp , &strn , &endn );
    CPXNETchgobj( env , net , 1 , &temp , &cst );
    CPXNETchgbds( env , net , 1 , &temp , "U" , &cpct );
    }
   else {
    QPchgarcnode( temp, strn, endn );
    CPXchgobj( env , qp , 1 , &temp , &cst );
    CPXchgbds( env , qp , 1 , &temp , "U" , &cpct );

    Startn[ temp ] = strn;
    Endn[ temp ] = endn;
    }

   while( ( FreePos < mmax ) && ( ArcPos[ FreePos ] < Inf<FNumber>() ) )
    FreePos++;

   return( temp );
   }
 #else
  throw(
   MCFException( "MCFCplex::AddArc() not implemented if DYNMC_MCF_CPX == 0"
		 ) );

  return( Inf<Index>() );
 #endif

 }  // end( MCFCplex::AddArc )

/*-------------------------------------------------------------------------*/

void MCFCplex::ChangeArc( cIndex name , cIndex nSN , cIndex nEN ) 
{
 int temp = int( name );
 int sn , en;
 if( ( nSN == Inf<Index>() ) || ( nEN == Inf<Index>() ) ) {
  if( net ) {
   CPXNETgetarcnodes( env, net, &sn, &en, temp, temp );
  } else {
   sn = Startn[ name ];
   en = Endn[ name ];
  }
 }

 if( nSN < Inf<Index>() )
  sn = int( nSN MINUSONE );

 if( nEN < Inf<Index>() )
  en = int( nEN MINUSONE );

 if( net ) 
  CPXNETchgarcnodes( env , net , 1 , &temp , &sn , &en ) ;
 else
  QPchgarcnode( temp , sn , en );
       
 }  // end( MCFCplex::ChangeArc )

/*--------------------------------------------------------------------------*/
/*---------------------------- DESTRUCTOR ----------------------------------*/
/*-------------------------------------------------------------------------*/

MCFCplex::~MCFCplex()
{
 if( mmax && nmax )
  MemDeAlloc();

 if( ! --InstCntr ) {
  delete[] ChangeUB;
  delete[] val;
  delete[] ind;
  nmaxarc = 0;
  }

 if( env == genv )
  if( ! --EnvICntr ) {
   CPXcloseCPLEX( &genv );
   genv = NULL;
   }

 }  // end( MCFCplex::~MCFCplex )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::MemAlloc( void )
{
 // create a new Cplex problem- - - - - - - - - - - - - - - - - - - - - - - -

 int ts;
 net = CPXNETcreateprob( env , &ts , "NET" );
 if( ( ! net ) || ts )
  throw( MCFException( "Problem creating Cplex problem" ) );

 // create data structures for arc creation/deletion- - - - - - - - - - - - -

 #if( DYNMC_MCF_CPX )
  ArcPos = new FNumber[ mmax ];
 #endif

 // if necessary, resize static members - - - - - - - - - - - - - - - - - - -
 
 if( mmax > nmaxarc ) {
  if( nmaxarc ) {
   delete[] ChangeUB;
   delete[] val;
   delete[] ind;
   }

  ind = new int[ nmaxarc = mmax ];
  val = new double[ nmaxarc ]; 
  ChangeUB = new char[ nmaxarc ];

  VectAssign( ChangeUB , 'U' , nmaxarc );
  }
 }  // end( MemAlloc )

/*--------------------------------------------------------------------------*/

void MCFCplex::MemDeAlloc( void )
{
 #if( DYNMC_MCF_CPX )
  delete[] ArcPos;
 #endif

 if( net ) 
  CPXNETfreeprob( env , &net ); 
 else {
  CPXfreeprob( env , &qp );
  delete[] Startn;
  delete[] Endn;
  }
 }  // end( MemDeAlloc )

/*--------------------------------------------------------------------------*/

void MCFCplex::TurnToQP( void )
{
 int status_p;
 qp = CPXcreateprob( env , &status_p , "QP" );
 if( ( ! qp ) || status_p )
  throw( MCFException( "Problem creating Cplex LP" ) );

 Startn = new int[ mmax ];
 Endn   = new int[ mmax ];

 CPXNETgetarcnodes( env , net , Startn , Endn , 0 , m - 1);

 status_p = CPXcopynettolp( env , qp , net );
 if( ( ! qp ) || status_p )
  throw( MCFException( "Problem creating copying NET to LP" ) );

 CPXNETfreeprob( env , &net );
 net = NULL;

 status_p = CPXchgprobtype( env , qp , CPXPROB_QP );
 if( status_p )
  throw( MCFException( "Problem changing problem type to QP" ) );
 
 status_p = CPXsetintparam( env , CPX_PARAM_QPMETHOD , QPMthd );
 if( status_p )
  throw( MCFException( "Problem setting QP solving method" ) );
 }

/*--------------------------------------------------------------------------*/

void MCFCplex::QPchgarcnode( int name , int sn , int en ) 
{
 // Modify node-arc incidence matrix - - - - - - - - - - - - - - - - - - - - -

 CPXchgcoef( env , qp , Startn[ name ] , name , 0 );
 CPXchgcoef( env , qp , Endn[ name ] , name , 0 );
 // set to 0 components of the node-arc incidence matrix corresponding to
 // arc "name"

 CPXchgcoef( env , qp , sn , name , 1 );  // setup new
 CPXchgcoef( env , qp , en , name , -1 ); // components

 Startn[ name ] = sn; 
 Endn[ name ]   = en;
 }

/*-------------------------------------------------------------------------*/
/*------------------initialize static members------------------------------*/
/*-------------------------------------------------------------------------*/

CPXENVptr MCFCplex::genv = 0;

int* MCFCplex::ind = NULL;
double* MCFCplex::val = NULL;
char* MCFCplex::ChangeUB = NULL;

MCFCplex::Index MCFCplex::InstCntr = 0;
MCFCplex::Index MCFCplex::EnvICntr = 0;
MCFCplex::Index MCFCplex::nmaxarc = 0;

/*--------------------------------------------------------------------------*/
/*-------------------- End File MCFCplex.C ---------------------------------*/
/*--------------------------------------------------------------------------*/
