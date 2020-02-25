/*--------------------------------------------------------------------------*/
/*------------------------- File MCFCplex.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Quadratic Min Cost Flow problems solver, based on calls to the Cplex
 * Callable Libraries. Conforms to the standard MCF interface
 * defined in MCFClass.h
 *
 * \version 1.42
 *
 * \date 16 - 10 - 2018
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Manca \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Matteo Sammartino \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 1997 - 2018 by Antonio Frangioni.
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFCplex
 #define __MCFCplex  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*----------------------------- INCLUDES -----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

#include <ilcplex/cplex.h>

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFCPLEX_MACROS Compile-time switches in MCFCplex.h
    These macros control some important details of the implementation.
    Although using macros for activating features of the implementation is
    not very C++, switching off some unused features may make the code
    more efficient in running time or memory.
    @{ */

#define CNUMBER_IS_DOUBLE 1

/**< Tells if CNumber is in fact double.
   Although the MCFClass interface is designed to work seamlessly for every
   possibile choice of the basic types FNumber and CNumber, Cplex only works
   with doubles. Thus, when CNumber != double some conversions have to be done
   between the internal Cplex format and the CNumbers that are received in
   input/expected in output; when CNumber == double this can be saved, and
   setting CNUMBER_IS_DOUBLE == 1 heqps the code to perform some operations
   faster and using less memory. */

#define FNUMBER_IS_DOUBLE 1

/**< Tells if FNumber is in fact double.
   Although the MCFClass interface is designed to work seamlessly for every
   possibile choice of the basic types FNumber and CNumber, Cplex only works
   with doubles. Thus, when FNumber != double some conversions have to be done
   between the internal Cplex format and the FNumbers that are received in
   input/expected in output; when FNumber == double this can be saved, and
   setting FNUMBER_IS_DOUBLE == 1 heqps the code to perform some operations
   faster and using less memory. */

#define INDEX_IS_UINT 1

/**< Tells if Index is in fact [unsigned] int.
   Indices in Cplex are ints; for Index == [unsigned] int, some conversions
   can be avoided [note: this assumes that unsigned ints and ints are
   phisically the same type], and setting INDEX_IS_UINT == 1 does that. */

/*----------------------------- DYNMC_MCF_CPX ------------------------------*/

#define DYNMC_MCF_CPX 1

/**< Decides if the graph topology (arcs, nodes) can be changed.
   If DYNMC_MCF_CPX > 0, some of the methods of the public interface of
   class that allow to change the topology of the underlying network are
   actually implemented. Possible values of this macro are:

   - 0 => arcs cannot be added or deleted, closed arcs cannot be reopened;
     all the other operations are possible;

   - 1 => all the methods that change the topology of the graph are
          implemented. */

/*@} -----------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFCplex --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFCPLEX_CLASSES Classes in MCFCplex.h
    @{ */

/** The MCFCplex class derives from the abstract base class MCFClass, thus
    sharing its (standard) interface, and solves (Linear) Min Cost Flow
    problems via calls to Cplex Callable Library functions. */

class MCFCplex: public MCFClass {

/*--------------------------------------------------------------------------*/
/*-------------------PUBLIC PART OF THE CLASS-------------------------------*/
/*--                                                                     ---*/
/*--  The following methods and data are the actual interface of the     ---*/
/*--  class: the standard user should use these methods and data only    ---*/
/*--                                                                     ---*/
/*--------------------------------------------------------------------------*/
   
 public:

/*--------------------------------------------------------------------------*/
/*------------------------- PUBLIC DATA STRUCTURES -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** Public enum describing the possible parameters of the MCF solver,
    "extended" from MCFClass::MCFParam, to be used with the methods
    SetPar() and GetPar(). */

  enum MCFCParam { kQPMethod = kLastParam     ///< solution method
                   };

/*--------------------------------------------------------------------------*/
/** enum describing possible ways for solving QP problem: see the Cplex
    manual for details */

   enum QPMethod { qpAutomatic = 0,
		   qpPSimplex  = 1,
		   qpDSimplex  = 2,
		   qpNSimplex  = 3,
		   qpBarrier   = 4
                   }; 

/*--------------------------------------------------------------------------*/
/*------------------------PUBLIC METHODS------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------CONSTRUCTOR--------------------------------------*/
/*--------------------------------------------------------------------------*/

   MCFCplex( cIndex nmx = 0 , cIndex mmx = 0 , CPXENVptr extenv = NULL );

/**< Constructor of the class.

   For the meaning of nmx and mmx see MCFClass::MCFClass().

   If extenv != NULL, it is taken as a pointer to a valid Cplex environment
   [see GetCplexEnv() below] that will be used for all the lifetime of the
   object. Otherwise, a Cplex environment is possibly initialized in the
   constructor. The environment is shared among all the active instances of
   MCFCplex objects; this is done in order to save on the number of (costly)
   Cplex licenses required to have multiple MCF solvers active at the same
   time, since any environment consumes a licence. Thus, the environment is
   actually initialized only when the first instance is constructed, and it
   is released when the last destructor (of an instance using it) is invoked.

   The possibility of passing an external environment is let because having
   all the instances to share the same environment has the drawback that
   all changes made to optimization parameters of the environment [see
   SetPar() below] are "seen" be all active instances, even though the
   changes are invoked for one specific object. If for some reason this is
   unacceptable, the user should provide each "sensitive" instance with its
   own private environment, letting all the others to share the "static" one.
   Of course, it is then user's responsibility to initialize and free the
   environment. */

/*--------------------------------------------------------------------------*/
/*---------------------- OTHER INITIALIZATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
		 cIndex pm = 0 , cFRow pU = NULL , cCRow pC = NULL ,
		 cFRow pDfct = NULL , cIndex_Set pSn = NULL ,
		 cIndex_Set pEn = NULL ) override;

/**< Inputs a new network, as in MCFClass::LoadNet().

   Passing pC[ i ] == C_INF means that the arc `i' does not exist in the
   problem. These arcs are just "closed" and their cost is set to 0: this is
   done for being (if DYNMC_MCF_CPX > 0) subsequently capable of "opening"
   them back with OpenArc(). If the corresponding pU[ i ] is == F_INF then
   the arc is just "deleted". */

/*--------------------------------------------------------------------------*/

   virtual inline void SetPar( int par , int val ) override;

/**< Set integer parameters of the algorithm.

   @param par   is the parameter to be set;

   @param value is the value to assign to the parameter.  

   Apart from the parameters of the base class, this method handles:

   - kQPMethod:   the alorithm used to solve the QP, possible values are
                  defined in the enum QPMethod.

   - <any other>: any unrecognized value is taken to be one of the the many
                  "int" algorithmic parameters of Cplex and passed right
                  away via CPXsetintparam() [see the documentation in the
                  Cplex manual for details. */

/*--------------------------------------------------------------------------*/

   virtual inline void SetPar( int par , double val ) override;

/**< Set float parameters of the algorithm.

   @param par   is the parameter to be set;

   @param value is the value to assign to the parameter.  

   Apart from the parameters of the base class, this method handles:

   - <any other>: any unrecognized value is taken to be one of the the many
                  "int" algorithmic parameters of Cplex and passed right
                  away via CPXsetintparam() [see the documentation in the
                  Cplex manual for details. */

/*--------------------------------------------------------------------------*/

   virtual inline void GetPar( int par , int &val ) const override
   {
    if( par == kQPMethod )
     val = (int) QPMthd;
    else
     MCFClass::GetPar( par , val );
    }

/**< This method returns one of the integer parameter of the algorithm.

   @param par  is the parameter to return [see SetPar( int ) for comments];

   @param val  upon return, it will contain the value of the parameter.

   Apart from the parameters of the base class, this method handles
   kQPMethod. */

/*--------------------------------------------------------------------------*/

   inline CPXENVptr GetCplexEnv( void ) const { return( env ); }

/**< Returns a pointer to the internal Cplex environment.

   This method is provided as an alternative to the two forms of SetPar()
   [see above]; by getting the pointer to the internal Cplex environment with
   GetCplexEnv(), CPXset***param() can be called directly. This also allows to
   perform any other operation with the environment, such as reading the value
   of the parameters with CPXgetintparam() and CPXgetdbqparam(), so care must
   be taken.

   The returned pointer is the same passed to the constructor [see above],
   if any; otherwise it is the "static" environment shared by all the active
   MCFCplex instances. In the latter case, any change in the environment
   simultaneously affect *all* the existing (and future) MCFCplex instances
   which have not been given a "private" environment. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void ) override;

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR READING RESULTS  -------------------------*/
/*--------------------------------------------------------------------------*/

  using MCFClass::MCFGetX;  // the ( void ) method, which is otherwise hidden

  void MCFGetX( FRow F , Index_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

  using MCFClass::MCFGetRC;  // the ( void ) method, which is otherwise hidden

  void MCFGetRC( CRow CR,  cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  CNumber MCFGetRC( cIndex i ) override;

/*--------------------------------------------------------------------------*/

  using MCFClass::MCFGetPi;  // the ( void ) method, which is otherwise hidden

  void MCFGetPi( CRow P , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

  FONumber MCFGetFO( void ) override;

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

  void MCFArcs( Index_Set Startv , Index_Set Endv , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  inline Index MCFSNde( cIndex i ) override;
     
  inline Index MCFENde( cIndex i ) override;

/*--------------------------------------------------------------------------*/
 
  void MCFCosts( CRow Costv , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;
    
  inline CNumber MCFCost( cIndex i ) override;

/*--------------------------------------------------------------------------*/

  void MCFQCoef( CRow Qv , cIndex_Set nms = NULL ,    
                 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  inline CNumber MCFQCoef( cIndex i ) override;

/*--------------------------------------------------------------------------*/

  void MCFUCaps( FRow UCapv , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  inline FNumber MCFUCap( cIndex i ) override;
     
/*--------------------------------------------------------------------------*/
     
  void MCFDfcts( FRow Dfctv , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  inline FNumber MCFDfct( cIndex i ) override;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

  void ChgCosts( cCRow NCost , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  void ChgCost( Index arc , cCNumber NCost ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  void ChgQCoef( cCRow NQCoef, cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  void ChgQCoef( Index arc , cCNumber NQCoef ) override;

/*--------------------------------------------------------------------------*/
  
  void ChgDfcts( cFRow NDfct , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  void ChgDfct( Index node , cFNumber NDfct ) override;

/*--------------------------------------------------------------------------*/

  void ChgUCaps( cFRow NCap , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

  void ChgUCap( Index arc , cFNumber NCap  ) override;


/*--------------------------------------------------------------------------*/
/*----------------- MODIFYING THE STRUCTURE OF THE GRAPH -------------------*/
/*--------------------------------------------------------------------------*/ 

  void CloseArc( cIndex name ) override;
     
  inline bool IsClosedArc( cIndex name ) override;

  void DelNode( cIndex name ) override;

  void OpenArc( cIndex name ) override;

  Index AddNode( cFNumber aDfct ) override;

  void ChangeArc( cIndex name ,
		  cIndex nSN = Inf<Index>() , cIndex nEN = Inf<Index>() )
   override;

  void DelArc( cIndex name ) override;

  inline bool IsDeletedArc( cIndex name ) override;

  Index AddArc( cIndex Start , cIndex End , cFNumber aU , cCNumber aC )
   override; 

  inline void printQ();

/*--------------------------------------------------------------------------*/
/*---------------------------- DESTRUCTOR ----------------------------------*/
/*--------------------------------------------------------------------------*/

  ~MCFCplex();  

/*--------------------------------------------------------------------------*/ 
/*----------------- PRIVATE PART OF THE CLASS ------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 void MemAlloc( void );

 void MemDeAlloc( void );

 void TurnToQP( void );

 void QPchgarcnode( int name , int sn , int en );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  CPXENVptr env;      // Cplex environment pointer 
  CPXNETptr net;      // network pointer
  CPXLPptr qp;        // QP pointer 

  int* Startn;        // arcs' start nodes
  int* Endn;          // arcs' end nodes
 
  QPMethod QPMthd;    // QP solving method

  #if( DYNMC_MCF_CPX )
   FRow ArcPos;       // ArcPos[ i ] == 0 means that arc i exists, == F_INF
                      // means that the position is available for creating
                      // a new arc, anything in between means that the arc
                      // is closed and that is its original capacity
   Index FreePos;     // first position available for creating a new arc, i.e.
                      // smallest index i s.t. ArcPos[ i ] == F_INF
  #endif

  // static members - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static CPXENVptr genv;      // "global" Cplex environment pointer 
  static Index InstCntr;      // active instances counter
  static Index EnvICntr;      // counter of active instances that use genv

  static int* ind;            // contains index of arcs for cplex methods
  static double* val;         // contains new values for cplex methods
  static char* ChangeUB;      // vector used in McfSolve::ChangeUCap(..)
  static Index nmaxarc;

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

  };  // end( class MCFCplex )

/* @} end( group( MCFCPLEX_CLASSES ) ) */
/*-------------------------------------------------------------------------*/
/*-------------------inline methods implementation-------------------------*/
/*-------------------------------------------------------------------------*/

inline void MCFCplex::SetPar( int par , int val )
{
 try {
  MCFClass::SetPar( par , val );  // is it handled by the base class?

  if( par == kMaxIter )           // let MaxIter to be handled by Cplex
   CPXsetintparam( env , CPX_PARAM_NETITLIM ,  val > 0 ? val : 2100000000 );
  }
 catch( MCFException &e ) {      // it is *not* handled by the base class
  if( par == kQPMethod )
   QPMthd = (QPMethod) val;
  else
   CPXsetintparam( env , par , val );
  }
 }

/*--------------------------------------------------------------------------*/

inline void MCFCplex::SetPar( int par , double val )
{
 try {
  MCFClass::SetPar( par , val );  // is it handled by the base class?

  if( par == kMaxTime )           // let MaxTime to be handled by Cplex
   CPXsetdblparam( env , CPX_PARAM_TILIM ,  val > 0 ? val : 1e+75 );
  }
 catch( MCFException &e ) {      // it is *not* handled by the base class
   CPXsetdblparam( env , par , val );
  }
 }

/*-------------------------------------------------------------------------*/

inline MCFCplex::Index MCFCplex::MCFSNde( MCFCplex::cIndex i )
{
 int strtn;
 if( net )
  CPXNETgetarcnodes( env , net , (int *) &strtn , NULL , int( i ) , int( i ) );
 else
  strtn = Startn[ i ];

 #if( ! USENAME0 )
  strtn++;
 #endif

 return( Index( strtn ) );
 }

/*-------------------------------------------------------------------------*/

inline MCFCplex::Index MCFCplex::MCFENde( MCFCplex::cIndex i )
{
 int endn;
 if( net )
  CPXNETgetarcnodes( env , net , NULL , (int *) &endn , int( i ) , int( i ) );
 else
  endn = Endn[ i ];

 #if( ! USENAME0 )
  endn++;
 #endif

 return( Index( endn ) );
 }

/*-------------------------------------------------------------------------*/

inline MCFCplex::CNumber MCFCplex::MCFCost( MCFCplex::cIndex i )
{
 double cst;
 if( net )  
  CPXNETgetobj( env , net , &cst , int( i ) , int( i ) );
 else
  CPXgetobj( env , qp , &cst , int( i ) , int( i ) );

 return( CNumber( cst ) );
 }

/*-------------------------------------------------------------------------*/

inline MCFCplex::CNumber MCFCplex::MCFQCoef( MCFCplex::cIndex i )
{
 double qcoef = 0;
 if( qp )
  CPXgetqpcoef( env , qp , int( i ) , int( i ) , &qcoef );

 return( qcoef );
 }

/*-------------------------------------------------------------------------*/

inline MCFCplex::FNumber MCFCplex::MCFUCap( MCFCplex::cIndex i )
{
 #if( DYNMC_MCF_CPX )
  if( ArcPos[ i ] && ( ArcPos[ i ] < Inf<FNumber>() ) )
   return( ArcPos[ i ] );
  else
 #endif
  {
   double ucap;
   if( net )
    CPXNETgetub( env , net , &ucap , int( i ) , int( i ) );
   else
    CPXgetub( env , qp , &ucap , int( i ) , int ( i ) );

   return( FNumber( ucap ) );
   }
 } 

/*-------------------------------------------------------------------------*/

inline MCFCplex::FNumber MCFCplex::MCFDfct( MCFCplex::cIndex i )
{
 double dfct;
 if( net )
  CPXNETgetsupply( env , net , &dfct , int( i ) , int( i ) );
 else
  CPXgetrhs( env , qp , &dfct , int( i ) , int( i ) );

 return( - FNumber( dfct ) ); 
 }

/*--------------------------------------------------------------------------*/

inline bool MCFCplex::IsClosedArc( MCFCplex::cIndex name )
{
 #if( DYNMC_MCF_CPX )
  return( ArcPos[ name ] && ( ArcPos[ name ] < Inf<Index>() ) );
 #else
  return( MCFCplex::MCFUCap( name ) == 0 );
 #endif
 }

/*--------------------------------------------------------------------------*/

inline bool MCFCplex::IsDeletedArc( MCFCplex::cIndex name )
{
 #if( DYNMC_MCF_CPX )
  return( ArcPos[ name ] == Inf<Index>() );
 #else
  return( false );
 #endif
 }

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MCFCplex.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File MCFCplex.h ---------------------------------*/
/*--------------------------------------------------------------------------*/
